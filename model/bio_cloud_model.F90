!======================================================================
! bio_cloud_model.F90
! Venus cloud biosphere model (Seager-cycle carrying-capacity diagnostic)
! ENVIRONMENT MODULE: the static cloud column + its physics.  The organism
! store, 3-state lifecycle and transport live in module bio_venus (which uses
! this one); the run driver is bio_venus_driver.F90.
!
! Provides:
!   cloud_read          : ingest the VPCM column (dynamics + master grid,
!                         and the VPCM mode-1/2 droplet fields)
!   cloud_assemble      : fill each mode's column (N, radius, width, dist
!                         type) from a chosen SOURCE per mode
!   cloud_build_spectrum: bin every active mode onto a shared radius grid
!   cloud_diagnostics / cloud_summary / cloud_spectrum_check : logging
!   host_capacity / host_supply / draw_host / host_capacity_sweep : packing
!   viscosity_co2 / mean_free_path / cunningham / settling_velocity(_z) :
!                         Stokes+slip settling (matched to VPCM new_cloud_sedim)
!   col_interp / dKzz_dz : continuous-z interpolation of column fields
!
! Each cloud particle mode m is described per level by:
!   N_m(z)   number density          [m^-3]
!   r_m(z)   radius parameter         [m]   (lognormal median, or Gaussian mean)
!   s_m(z)   width parameter                (lognormal sigma_g [-], or Gaussian sigma [m])
!   dist_m   distribution type (lognormal | gaussian)
! and is sourced INDEPENDENTLY from one of:
!   VPCM  : modes 1,2 from the cloudvenus two-moment column (this run's input)
!   HAUS  : Haus et al. 2016 analytic standard model (modes 1,2,3; lognormal)
!   KH80  : Knollenberg & Hunten 1980 LCPS Table 4 (modes 1,3 lognormal,
!           mode 2 Gaussian), read from inputs/kh80_table4.txt
! So one can run e.g. VPCM(1,2)+HAUS(3), all-HAUS, all-KH80, etc.
!
! Dynamics (T,P,rho,w,Kzz) and the altitude grid ALWAYS come from the VPCM
! column.  Acidity (WSA) is read but is NOT a habitability criterion in this
! project; temperature is the sole AHZ axis.
!
! Units: SI throughout (m, kg, s, K).  Altitude is km in files -> m internally.
!======================================================================
module bio_cloud_model
  implicit none
  private

  integer, parameter :: dp = 8

  ! ---- Mode and source/distribution identifiers ----
  integer, parameter, public :: NMODES = 3
  integer, parameter, public :: SRC_OFF = 0, SRC_VPCM = 1, SRC_HAUS = 2, SRC_KH80 = 3
  integer, parameter, public :: DIST_LOGNORMAL = 1, DIST_GAUSSIAN = 2

  ! ---- Cloud column grid + dynamics (always from the VPCM column) ----
  integer,  save              :: nz = 0
  real(dp), allocatable, save :: z(:)        ! altitude        [m]
  real(dp), allocatable, save :: T(:)        ! temperature     [K]
  real(dp), allocatable, save :: P(:)        ! pressure        [Pa]
  real(dp), allocatable, save :: rho_air(:)  ! gas density     [kg/m^3]
  real(dp), allocatable, save :: w(:)        ! vertical wind   [m/s]
  real(dp), allocatable, save :: Kzz(:)      ! eddy diffusivity[m^2/s]
  real(dp), allocatable, save :: WSA(:)      ! H2SO4 mass fraction   [-] (diagnostic only)
  real(dp), allocatable, save :: rho_host(:) ! droplet density       [kg/m^3]

  ! ---- Raw VPCM mode-1/2 droplet fields (used by the VPCM source) ----
  real(dp), allocatable, save :: vN1(:), vr1(:), vN2(:), vr2(:)
  logical,  allocatable, save :: vvalid1(:), vvalid2(:)

  ! ---- Assembled per-mode columns (filled by cloud_assemble) ----
  real(dp), allocatable, save :: Nm(:,:)     ! number density   [m^-3]  (nz, NMODES)
  real(dp), allocatable, save :: rm(:,:)     ! radius param     [m]     (nz, NMODES)
  real(dp), allocatable, save :: sm(:,:)     ! width param (sigma_g [-] or sigma [m]) (nz, NMODES)
  logical,  allocatable, save :: validm(:,:) ! mode active at level     (nz, NMODES)
  integer,  save              :: distm(NMODES) ! distribution type per mode
  integer,  save              :: srcm(NMODES)  ! source per mode

  ! ---- Reconstructed binned host spectrum ----
  integer,  parameter, public :: NBIN = 80
  real(dp), parameter         :: R_GRID_MIN = 1.0e-8_dp   ! grid lower edge [m] (0.01 um)
  real(dp), parameter         :: R_GRID_MAX = 5.0e-5_dp   ! grid upper edge [m] (50 um; covers mode 3)
  real(dp), allocatable, save :: r_edge(:)   ! bin edges   [m]  (NBIN+1)
  real(dp), allocatable, save :: r_bin(:)    ! bin centres [m]  (NBIN)
  real(dp), allocatable, save :: specm(:,:,:)! per-mode number/bin [m^-3] (NBIN, nz, NMODES)
  real(dp), allocatable, save :: spec(:,:)   ! total host spectrum [m^-3] (NBIN, nz)

  ! ---- VPCM fixed log-normal widths (conf_phys.F90 defaults) ----
  real(dp), parameter, public :: SIGMA1 = 1.56_dp   ! mode 1
  real(dp), parameter, public :: SIGMA2 = 1.29_dp   ! mode 2

  ! ---- VPCM host data-hygiene masking thresholds ----
  real(dp), parameter, public :: N_DROP_MIN = 1.0e3_dp   ! min trustworthy droplet number [m^-3]
  real(dp), parameter, public :: R1_MAX = 1.0e-6_dp      ! mode-1 max trusted radius [m] (~3 x ri1)
  real(dp), parameter, public :: R2_MAX = 3.0e-6_dp      ! mode-2 max trusted radius [m] (~3 x ri2)
  ! Generic minimum number density for a non-VPCM mode to count as "active".
  real(dp), parameter, public :: N_HOST_MIN = 1.0e3_dp   ! [m^-3]

  ! ---- Host-capacity engine ----
  ! A host droplet of radius r_host holds at most C = floor(phi*(r_host/r_cell)^3)
  ! cells (volume packing).  phi is the packing fraction (~0.6 random close packing;
  ! phi=1 reads C as available liquid volume).  Settable via the namelist.
  real(dp), parameter, public :: PHI_PACK_DEFAULT = 0.6_dp
  real(dp), save,      public :: phi_pack = PHI_PACK_DEFAULT

  ! ---- Settling physics (matched to VPCM cloudvenus/new_cloud_sedim.F) ----
  ! Stokes drag with Cunningham slip:
  !   v_settle = (2/9)*(rho_p - rho_air)*g*r^2/mu * Cc(r)
  !   Cc = 1 + (lambda/r)*[A + B*exp(-C*r/lambda)]   (A,B,C = 1.246,0.42,0.87)
  ! At um sizes in thin Venus air Kn ~ O(1), so slip (Cc>1) is non-negligible;
  ! pure Stokes overestimates drag and underestimates the fall speed.
  real(dp), parameter, public :: PI         = 3.141592653589793_dp
  real(dp), parameter, public :: GRAV_VENUS = 8.87_dp          ! g [m/s^2] (VPCM RG)
  real(dp), parameter :: R_UNIV     = 8.314462618_dp           ! universal gas const [J/mol/K]
  real(dp), parameter :: N_AVO      = 6.02214076e23_dp         ! Avogadro [1/mol]
  real(dp), parameter :: CO2_MOLRAD = 2.2e-10_dp               ! CO2 molecular radius [m] (VPCM)
  real(dp), parameter :: SLIP_A = 1.246_dp, SLIP_B = 0.42_dp, SLIP_C = 0.87_dp
  ! (The organism store, 3-state machine, spore density and altitude bands now
  !  live in module bio_venus, which uses this environment module.)

  ! ---- Haus et al. 2016 standard-model analytic parameters (per mode) ----
  ! N(z): piecewise top-hat with exponential tails (Haus Table 1 + Eq. T1).
  ! Modal radius (= lognormal median r0) and sigma_g from Pollack et al. 1993.
  real(dp), parameter :: HAUS_ZB(NMODES)   = [49.0_dp, 65.0_dp, 49.0_dp]   ! base [km]
  real(dp), parameter :: HAUS_ZC(NMODES)   = [16.0_dp,  1.0_dp,  8.0_dp]   ! const-N thickness [km]
  real(dp), parameter :: HAUS_HUP(NMODES)  = [ 3.5_dp,  3.5_dp,  1.0_dp]   ! upper scale height [km]
  real(dp), parameter :: HAUS_HLO(NMODES)  = [ 1.0_dp,  3.0_dp,  0.5_dp]   ! lower scale height [km]
  real(dp), parameter :: HAUS_N0(NMODES)   = [193.5_dp,100.0_dp, 14.0_dp]  ! peak number [cm^-3]
  real(dp), parameter :: HAUS_R0(NMODES)   = [0.30_dp,  1.00_dp, 3.65_dp]  ! median radius [um]
  real(dp), parameter :: HAUS_SIG(NMODES)  = [1.56_dp,  1.29_dp, 1.28_dp]  ! sigma_g [-]

  public :: cloud_read, cloud_assemble, cloud_build_spectrum
  public :: cloud_diagnostics, cloud_summary, cloud_spectrum_check, cloud_cleanup
  public :: src_name, set_source
  public :: set_phi
  public :: host_capacity, host_supply, draw_host, host_capacity_sweep
  public :: viscosity_co2, mean_free_path, cunningham
  public :: settling_velocity, settling_velocity_z, col_interp, dKzz_dz
  public :: nz, z, T, P, rho_air, w, Kzz, WSA, rho_host
  public :: Nm, rm, sm, validm, distm, srcm
  public :: r_bin, r_edge, specm, spec

contains

  !====================================================================
  ! Low-level column file I/O
  !====================================================================
  integer function count_rows(path) result(n)
    character(*), intent(in) :: path
    integer :: u, ios
    character(len=512) :: line
    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(a)') 'ERROR: cannot open '//trim(path); stop 1
    end if
    n = 0
    do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#')    cycle
      n = n + 1
    end do
    close(u)
  end function count_rows

  subroutine read_column(path, zc, vc)
    character(*), intent(in)  :: path
    real(dp),     intent(out) :: zc(:), vc(:)
    integer :: u, ios, k
    character(len=512) :: line
    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(a)') 'ERROR: cannot open '//trim(path); stop 1
    end if
    k = 0
    do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#')    cycle
      k = k + 1
      if (k > size(zc)) exit
      read(line,*,iostat=ios) zc(k), vc(k)
      if (ios /= 0) then
        write(*,'(a,i0)') 'ERROR: parse failure in '//trim(path)//' at data row ', k; stop 1
      end if
    end do
    close(u)
    if (k /= size(zc)) then
      write(*,'(a,i0,a,i0)') 'ERROR: '//trim(path)//' row count ', k, ' /= expected ', size(zc); stop 1
    end if
  end subroutine read_column

  subroutine read_field(indir, fname, dest, zkm_ref)
    character(*), intent(in)  :: indir, fname
    real(dp),     intent(out) :: dest(:)
    real(dp),     intent(in)  :: zkm_ref(:)
    real(dp) :: zkm(size(dest)), dzmax
    call read_column(trim(indir)//'/'//trim(fname), zkm, dest)
    dzmax = maxval(abs(zkm - zkm_ref))
    if (dzmax > 1.0e-3_dp) then
      write(*,'(a,es10.3)') 'ERROR: '//trim(fname)//' altitude grid mismatch, max |dz| [km] = ', dzmax
      stop 1
    end if
  end subroutine read_field

  !====================================================================
  ! Read the VPCM column: dynamics + master grid + raw mode-1/2 fields.
  !====================================================================
  subroutine cloud_read(indir)
    character(*), intent(in) :: indir
    real(dp), allocatable :: zkm(:)
    integer :: k

    nz = count_rows(trim(indir)//'/vpcm_temp.txt')
    if (nz < 2) then
      write(*,'(a)') 'ERROR: fewer than 2 data rows in vpcm_temp.txt'; stop 1
    end if

    call cloud_cleanup()
    allocate(z(nz), T(nz), P(nz), rho_air(nz), w(nz), Kzz(nz), WSA(nz), rho_host(nz))
    allocate(vN1(nz), vr1(nz), vN2(nz), vr2(nz), vvalid1(nz), vvalid2(nz))
    allocate(Nm(nz,NMODES), rm(nz,NMODES), sm(nz,NMODES), validm(nz,NMODES))
    allocate(zkm(nz))

    call read_column(trim(indir)//'/vpcm_temp.txt', zkm, T)
    z = zkm * 1000.0_dp

    call read_field(indir, 'vpcm_pres.txt',         P,        zkm)
    call read_field(indir, 'vpcm_rho.txt',          rho_air,  zkm)
    call read_field(indir, 'vpcm_vitwz.txt',        w,        zkm)
    call read_field(indir, 'vpcm_kz.txt',           Kzz,      zkm)
    call read_field(indir, 'vpcm_m0_mode1drop.txt', vN1,      zkm)
    call read_field(indir, 'vpcm_r1.txt',           vr1,      zkm)
    call read_field(indir, 'vpcm_m0_mode2drop.txt', vN2,      zkm)
    call read_field(indir, 'vpcm_r2.txt',           vr2,      zkm)
    call read_field(indir, 'vpcm_WSAVtab.txt',      WSA,      zkm)
    call read_field(indir, 'vpcm_rhotab.txt',       rho_host, zkm)

    ! VPCM host data-hygiene mask (number threshold + per-mode radius guard).
    do k = 1, nz
      vvalid1(k) = (vN1(k) > N_DROP_MIN) .and. (vr1(k) > 0.0_dp) .and. (vr1(k) < R1_MAX)
      vvalid2(k) = (vN2(k) > N_DROP_MIN) .and. (vr2(k) > 0.0_dp) .and. (vr2(k) < R2_MAX)
    end do

    ! Default sources: VPCM for modes 1&2, mode 3 off (reproduces prior behaviour).
    srcm = [SRC_VPCM, SRC_VPCM, SRC_OFF]

    deallocate(zkm)
  end subroutine cloud_read

  !====================================================================
  ! Source selection helpers
  !====================================================================
  subroutine set_source(imode, isrc)
    integer, intent(in) :: imode, isrc
    if (imode < 1 .or. imode > NMODES) then
      write(*,'(a,i0)') 'ERROR: set_source bad mode ', imode; stop 1
    end if
    srcm(imode) = isrc
  end subroutine set_source

  pure function src_name(isrc) result(nm)
    integer, intent(in) :: isrc
    character(len=4) :: nm
    select case (isrc)
    case (SRC_OFF);  nm = 'off '
    case (SRC_VPCM); nm = 'vpcm'
    case (SRC_HAUS); nm = 'haus'
    case (SRC_KH80); nm = 'kh80'
    case default;    nm = '??? '
    end select
  end function src_name

  !====================================================================
  ! Assemble each mode's column from its chosen source.
  !====================================================================
  subroutine cloud_assemble(indir)
    character(*), intent(in) :: indir
    integer :: m
    if (nz <= 0) then
      write(*,'(a)') 'ERROR: cloud_assemble called before cloud_read'; stop 1
    end if
    Nm = 0.0_dp; rm = 0.0_dp; sm = 0.0_dp; validm = .false.
    distm = DIST_LOGNORMAL
    do m = 1, NMODES
      select case (srcm(m))
      case (SRC_OFF)
        ! leave empty
      case (SRC_VPCM)
        call adapt_vpcm(m)
      case (SRC_HAUS)
        call adapt_haus(m)
      case (SRC_KH80)
        call adapt_kh80(m, indir)
      case default
        write(*,'(a,i0)') 'ERROR: unknown source for mode ', m; stop 1
      end select
    end do
  end subroutine cloud_assemble

  !--------------------------------------------------------------------
  ! VPCM source: modes 1 and 2 only (lognormal, fixed sigma).
  !--------------------------------------------------------------------
  subroutine adapt_vpcm(m)
    integer, intent(in) :: m
    select case (m)
    case (1)
      Nm(:,1) = vN1; rm(:,1) = vr1; sm(:,1) = SIGMA1
      validm(:,1) = vvalid1; distm(1) = DIST_LOGNORMAL
    case (2)
      Nm(:,2) = vN2; rm(:,2) = vr2; sm(:,2) = SIGMA2
      validm(:,2) = vvalid2; distm(2) = DIST_LOGNORMAL
    case default
      write(*,'(a)') 'ERROR: VPCM source has no mode 3 (use HAUS or KH80)'; stop 1
    end select
  end subroutine adapt_vpcm

  !--------------------------------------------------------------------
  ! Haus 2016 source: analytic lognormal profile (modes 1,2,3).
  !--------------------------------------------------------------------
  subroutine adapt_haus(m)
    integer, intent(in) :: m
    integer  :: k
    real(dp) :: zkm, zb, zc, hup, hlo, n0, Nloc
    zb = HAUS_ZB(m); zc = HAUS_ZC(m); hup = HAUS_HUP(m); hlo = HAUS_HLO(m); n0 = HAUS_N0(m)
    distm(m) = DIST_LOGNORMAL
    do k = 1, nz
      zkm = z(k) / 1000.0_dp
      if (zkm > zb + zc) then
        Nloc = n0 * exp(-(zkm - (zb + zc)) / hup)
      else if (zkm >= zb) then
        Nloc = n0
      else
        Nloc = n0 * exp(-(zb - zkm) / hlo)
      end if
      Nm(k,m)     = Nloc * 1.0e6_dp              ! cm^-3 -> m^-3
      rm(k,m)     = HAUS_R0(m) * 1.0e-6_dp       ! um -> m (lognormal median)
      sm(k,m)     = HAUS_SIG(m)                  ! sigma_g [-]
      validm(k,m) = Nm(k,m) > N_HOST_MIN
    end do
  end subroutine adapt_haus

  !--------------------------------------------------------------------
  ! KH80 source: tabulated LCPS Table 4, sampled piecewise-constant onto the z grid.
  ! Modes 1,3 lognormal; mode 2 Gaussian.  Size columns are DIAMETERS.
  ! File columns: z_km m1_N m1_sigg m1_Dum  m2_N m2_sigD m2_Dmean  m3_N m3_sigg m3_Dum
  !--------------------------------------------------------------------
  subroutine adapt_kh80(m, indir)
    integer,      intent(in) :: m
    character(*), intent(in) :: indir
    character(len=256) :: path
    integer :: u, ios, nt, k, j
    real(dp), allocatable :: zt(:), col(:,:)
    real(dp) :: tmp(10), Nval, wval, rval
    character(len=512) :: line

    path = trim(indir)//'/kh80_table4.txt'
    nt = count_rows(path)
    allocate(zt(nt), col(nt,9))
    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(a)') 'ERROR: cannot open '//trim(path); stop 1
    end if
    j = 0
    do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      j = j + 1
      read(line,*,iostat=ios) tmp(1:10)
      if (ios /= 0) then
        write(*,'(a,i0)') 'ERROR: KH80 parse failure at row ', j; stop 1
      end if
      zt(j)    = tmp(1)
      col(j,:) = tmp(2:10)
    end do
    close(u)
    ! File is ascending in z (generator sorted it); assume so.

    select case (m)
    case (1); distm(1) = DIST_LOGNORMAL   ! cols 1-3: N, sigma_g, D[um]
    case (2); distm(2) = DIST_GAUSSIAN    ! cols 4-6: N, sigma_D[um], D_mean[um]
    case (3); distm(3) = DIST_LOGNORMAL   ! cols 7-9: N, sigma_g, D[um]
    end select

    do k = 1, nz
      call kh80_interp(zt, col, nt, z(k)/1000.0_dp, m, Nval, wval, rval)
      Nm(k,m) = Nval * 1.0e6_dp                  ! cm^-3 -> m^-3
      select case (m)
      case (1, 3)
        rm(k,m) = 0.5_dp * rval * 1.0e-6_dp      ! diameter um -> radius m
        sm(k,m) = wval                           ! sigma_g [-]
      case (2)
        rm(k,m) = 0.5_dp * rval * 1.0e-6_dp      ! mean diameter um -> mean radius m
        sm(k,m) = 0.5_dp * wval * 1.0e-6_dp      ! sigma diameter um -> sigma radius m
      end select
      validm(k,m) = Nm(k,m) > N_HOST_MIN .and. rm(k,m) > 0.0_dp
    end do
    deallocate(zt, col)
  end subroutine adapt_kh80

  ! Piecewise-constant sampling of KH80 columns for mode m at altitude zq [km].
  ! The LCPS values are bin AVERAGES over altitude intervals, so we treat each
  ! row as a constant value owned by its bin and assign zq to the nearest bin
  ! centre (nearest-neighbour in altitude).  Centres are used rather than the
  ! reported edges because the digitized edges have OCR typos in the 50-53 km
  ! region whereas the centres are stable.  Returns 0 outside the measured span
  ! or where the nearest row has the mode absent (N <= 0).
  subroutine kh80_interp(zt, col, nt, zq, m, Nval, wval, rval)
    real(dp), intent(in)  :: zt(:), col(:,:), zq
    integer,  intent(in)  :: nt, m
    real(dp), intent(out) :: Nval, wval, rval
    integer :: j, c0, jbest
    real(dp) :: dbest, dd
    select case (m)
    case (1); c0 = 1
    case (2); c0 = 4
    case (3); c0 = 7
    case default; c0 = 1
    end select
    Nval = 0.0_dp; wval = 0.0_dp; rval = 0.0_dp
    if (zq < zt(1) .or. zq > zt(nt)) return
    jbest = 1; dbest = abs(zq - zt(1))
    do j = 2, nt
      dd = abs(zq - zt(j))
      if (dd < dbest) then; dbest = dd; jbest = j; end if
    end do
    if (col(jbest,c0) <= 0.0_dp) return        ! mode absent at nearest bin
    Nval = col(jbest, c0)
    wval = col(jbest, c0+1)
    rval = col(jbest, c0+2)
  end subroutine kh80_interp

  !====================================================================
  ! Binned spectrum reconstruction
  !====================================================================
  ! Lognormal: median r0, geometric width sigma_g.  Exact CDF (erf), renormalised.
  subroutine lognormal_bins(Ntot, r0, sigma_g, out)
    real(dp), intent(in)  :: Ntot, r0, sigma_g
    real(dp), intent(out) :: out(:)
    real(dp) :: s, sq2, denom, lo, hi
    integer  :: i
    out = 0.0_dp
    if (Ntot <= 0.0_dp .or. r0 <= 0.0_dp .or. sigma_g <= 1.0_dp) return
    s   = log(sigma_g)
    sq2 = s * sqrt(2.0_dp)
    denom = 0.5_dp * ( erf(log(R_GRID_MAX/r0)/sq2) - erf(log(R_GRID_MIN/r0)/sq2) )
    do i = 1, NBIN
      lo = 0.5_dp * erf(log(r_edge(i)  /r0)/sq2)
      hi = 0.5_dp * erf(log(r_edge(i+1)/r0)/sq2)
      out(i) = Ntot * (hi - lo) / denom
    end do
  end subroutine lognormal_bins

  ! Gaussian in radius: mean r_mean, std sigma_r.  Exact CDF (erf), renormalised
  ! over the (positive) grid.
  subroutine gaussian_bins(Ntot, r_mean, sigma_r, out)
    real(dp), intent(in)  :: Ntot, r_mean, sigma_r
    real(dp), intent(out) :: out(:)
    real(dp) :: sq2, denom, lo, hi
    integer  :: i
    out = 0.0_dp
    if (Ntot <= 0.0_dp .or. sigma_r <= 0.0_dp) return
    sq2 = sigma_r * sqrt(2.0_dp)
    denom = 0.5_dp * ( erf((R_GRID_MAX-r_mean)/sq2) - erf((R_GRID_MIN-r_mean)/sq2) )
    if (denom <= 0.0_dp) return
    do i = 1, NBIN
      lo = 0.5_dp * erf((r_edge(i)  -r_mean)/sq2)
      hi = 0.5_dp * erf((r_edge(i+1)-r_mean)/sq2)
      out(i) = Ntot * (hi - lo) / denom
    end do
  end subroutine gaussian_bins

  subroutine cloud_build_spectrum()
    integer  :: i, k, m
    real(dp) :: ratio, tmp(NBIN)

    if (nz <= 0) then
      write(*,'(a)') 'ERROR: cloud_build_spectrum called before cloud_read'; stop 1
    end if
    if (.not. allocated(Nm)) then
      write(*,'(a)') 'ERROR: cloud_build_spectrum called before cloud_assemble'; stop 1
    end if

    if (allocated(r_edge)) deallocate(r_edge)
    if (allocated(r_bin))  deallocate(r_bin)
    if (allocated(specm))  deallocate(specm)
    if (allocated(spec))   deallocate(spec)
    allocate(r_edge(NBIN+1), r_bin(NBIN), specm(NBIN, nz, NMODES), spec(NBIN, nz))

    ratio = (R_GRID_MAX / R_GRID_MIN) ** (1.0_dp / real(NBIN, dp))
    do i = 1, NBIN+1
      r_edge(i) = R_GRID_MIN * ratio ** real(i-1, dp)
    end do
    do i = 1, NBIN
      r_bin(i) = sqrt(r_edge(i) * r_edge(i+1))
    end do

    specm = 0.0_dp
    do m = 1, NMODES
      if (srcm(m) == SRC_OFF) cycle
      do k = 1, nz
        if (.not. validm(k,m)) cycle
        if (distm(m) == DIST_LOGNORMAL) then
          call lognormal_bins(Nm(k,m), rm(k,m), sm(k,m), tmp)
        else
          call gaussian_bins(Nm(k,m), rm(k,m), sm(k,m), tmp)
        end if
        specm(:,k,m) = tmp
      end do
    end do
    do k = 1, nz
      spec(:,k) = 0.0_dp
      do m = 1, NMODES
        spec(:,k) = spec(:,k) + specm(:,k,m)
      end do
    end do
  end subroutine cloud_build_spectrum

  !====================================================================
  ! Diagnostics
  !====================================================================
  subroutine cloud_diagnostics(unit)
    integer, intent(in) :: unit
    integer :: k, m
    write(unit,'(a)') '# Venus cloud column — assembled per-mode profiles'
    write(unit,'(a,i0)') '# levels = ', nz
    do m = 1, NMODES
      write(unit,'(a,i0,a,a,a,i0)') '#   mode ', m, ' source=', src_name(srcm(m)), &
        '  dist=', distm(m)
    end do
    write(unit,'(a)') '#   (per mode: N[m-3]  r[um]  v=active flag)'
    write(unit,'(a)') '#   k   z[km]    T[K]    '// &
      'N1[m-3]    r1[um] v1   N2[m-3]    r2[um] v2   N3[m-3]    r3[um] v3'
    do k = 1, nz
      write(unit,'(i5,1x,f8.3,1x,f8.2,3(2x,es10.3,1x,f7.3,1x,i1))') &
        k, z(k)/1000.0_dp, T(k), &
        (Nm(k,m), rm(k,m)*1.0e6_dp, merge(1,0,validm(k,m)), m=1,NMODES)
    end do
  end subroutine cloud_diagnostics

  subroutine cloud_summary(unit)
    integer, intent(in) :: unit
    integer  :: m, nv, kb, ktop, ipk
    logical  :: anyhost(nz)

    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a)') ' CLOUD COLUMN SUMMARY (assembled)'
    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a,i0)') '  levels read           : ', nz
    write(unit,'(a,f7.2,a,f7.2,a)') '  altitude range        : ', &
         z(1)/1000.0_dp, ' to ', z(nz)/1000.0_dp, ' km'
    write(unit,'(a,f7.1,a,f7.1,a)') '  temperature range     : ', &
         minval(T), ' to ', maxval(T), ' K'

    anyhost = .false.
    do m = 1, NMODES
      write(unit,'(a)') ' '
      write(unit,'(a,i0,a,a,a,i0)') '  mode ', m, '  source=', src_name(srcm(m)), &
           '  dist=', distm(m)
      if (srcm(m) == SRC_OFF) then
        write(unit,'(a)') '     (off)'
        cycle
      end if
      nv = count(validm(:,m))
      write(unit,'(a,i0)') '     active levels      : ', nv
      if (nv > 0) then
        kb   = minloc(z, dim=1, mask=validm(:,m))
        ktop = maxloc(z, dim=1, mask=validm(:,m))
        ipk  = maxloc(Nm(:,m), dim=1, mask=validm(:,m))
        write(unit,'(a,f7.2,a,f7.2,a)') '     altitude window     : ', &
             z(kb)/1000.0_dp, ' to ', z(ktop)/1000.0_dp, ' km'
        write(unit,'(a,es10.3,a,f6.2,a)') '     peak N              : ', &
             Nm(ipk,m), ' m-3 at ', z(ipk)/1000.0_dp, ' km'
        write(unit,'(a,f8.4,a,f8.4,a)') '     radius range        : ', &
             minval(rm(:,m), mask=validm(:,m))*1.0e6_dp, ' to ', &
             maxval(rm(:,m), mask=validm(:,m))*1.0e6_dp, ' um'
        anyhost = anyhost .or. validm(:,m)
      end if
    end do

    write(unit,'(a)') ' '
    if (count(anyhost) > 0) then
      kb   = minloc(z, dim=1, mask=anyhost)
      ktop = maxloc(z, dim=1, mask=anyhost)
      write(unit,'(a,f7.2,a,f7.2,a)') '  combined host window  : ', &
           z(kb)/1000.0_dp, ' to ', z(ktop)/1000.0_dp, ' km'
      write(unit,'(a,f5.3,a,f5.3,a)') '  WSA range (in cloud)  : ', &
           minval(WSA, mask=anyhost), ' to ', maxval(WSA, mask=anyhost), &
           '  (diagnostic only)'
    else
      write(unit,'(a)') '  combined host window  : NONE'
    end if
    write(unit,'(a)') '----------------------------------------------------------'
  end subroutine cloud_summary

  subroutine cloud_spectrum_check(unit)
    integer, intent(in) :: unit
    integer  :: k, i, m, kpk
    real(dp) :: nrec, rg, errN, errR, maxN, maxR, Npk

    write(unit,'(a)') ' '
    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a)') ' HOST SPECTRUM RECONSTRUCTION CHECK'
    write(unit,'(a,i0,a,f6.3,a,f6.3,a)') '  grid: ', NBIN, ' log bins, ', &
         r_edge(1)*1.0e6_dp, ' to ', r_edge(NBIN+1)*1.0e6_dp, ' um'
    write(unit,'(a)') '  (errN = recovered/input number - 1; errR = recovered/input mean radius - 1)'

    maxN = 0.0_dp; maxR = 0.0_dp
    do m = 1, NMODES
      if (srcm(m) == SRC_OFF) cycle
      do k = 1, nz
        if (.not. validm(k,m)) cycle
        nrec = sum(specm(:,k,m))
        if (nrec <= 0.0_dp) cycle
        ! Recover the radius using each mode's own central statistic, so the
        ! comparison to the input rm is apples-to-apples: geometric mean for
        ! lognormal modes (= input median r0), arithmetic mean for Gaussian
        ! modes (= input mean radius).
        if (distm(m) == DIST_GAUSSIAN) then
          rg = sum(specm(:,k,m) * r_bin) / nrec          ! recovered arithmetic mean radius
        else
          rg = exp(sum(specm(:,k,m) * log(r_bin)) / nrec)! recovered geometric mean radius
        end if
        errN = nrec / Nm(k,m) - 1.0_dp
        errR = rg   / rm(k,m) - 1.0_dp
        maxN = max(maxN, abs(errN))
        maxR = max(maxR, abs(errR))
      end do
    end do
    write(unit,'(a,es9.2,a,es9.2)') '  max |errN| = ', maxN, '   max |errR| = ', maxR
    if (maxN < 1.0e-3_dp .and. maxR < 5.0e-2_dp) then
      write(unit,'(a)') '  RECONSTRUCTION OK'
    else
      write(unit,'(a)') '  NOTE: residual errR can arise from grid truncation of a wide/edge-hugging mode'
    end if

    ! Spectrum dump at the level of peak total spectrum.
    Npk = 0.0_dp; kpk = 1
    do k = 1, nz
      if (sum(spec(:,k)) > Npk) then; Npk = sum(spec(:,k)); kpk = k; end if
    end do
    if (Npk > 0.0_dp) then
      write(unit,'(a)') ' '
      write(unit,'(a,f6.2,a)') '  combined host spectrum at peak level (', &
           z(kpk)/1000.0_dp, ' km), non-empty bins:'
      write(unit,'(a)') '#    r[um]      n_total[m-3]   (n1, n2, n3)'
      do i = 1, NBIN
        if (spec(i,kpk) > 1.0e-2_dp) &
          write(unit,'(2x,f10.4,2x,es12.4,3(1x,es10.2))') r_bin(i)*1.0e6_dp, spec(i,kpk), &
            specm(i,kpk,1), specm(i,kpk,2), specm(i,kpk,3)
      end do
    end if
    write(unit,'(a)') '----------------------------------------------------------'
  end subroutine cloud_spectrum_check

  !====================================================================
  ! Host-capacity engine
  !   Given a cell radius r_cell, query the reconstructed host spectrum
  !   spec(:,k) for droplets large enough to host it (r_host >= r_cell),
  !   draw one weighted by local abundance, and report packing capacity
  !   C = floor(phi*(r_host/r_cell)^3).  Pure queries on the static cloud;
  !   no time integration here.
  !====================================================================

  subroutine set_phi(val)
    real(dp), intent(in) :: val
    if (val > 0.0_dp) phi_pack = val
  end subroutine set_phi

  ! Volume-packing capacity of a host droplet for a given cell radius.
  ! floor() via int() (operand is positive); 0 if the cell does not fit.
  integer function host_capacity(r_host, r_cell) result(C)
    real(dp), intent(in) :: r_host, r_cell
    C = 0
    if (r_cell <= 0.0_dp .or. r_host < r_cell) return
    C = int(phi_pack * (r_host / r_cell)**3)
  end function host_capacity

  ! Total number density [m^-3] of eligible hosts (r_host >= r_cell) at level k.
  real(dp) function host_supply(k, r_cell) result(s)
    integer,  intent(in) :: k
    real(dp), intent(in) :: r_cell
    integer :: i
    s = 0.0_dp
    do i = 1, NBIN
      if (r_bin(i) >= r_cell) s = s + spec(i,k)
    end do
  end function host_supply

  ! Draw one host bin for a germinating cell at level k, weighted by abundance
  ! over the eligible (r_host >= r_cell) part of the spectrum.  Returns the bin
  ! centre radius r_host, its packing capacity, and a found flag.
  subroutine draw_host(k, r_cell, r_host, capac, found)
    integer,  intent(in)  :: k
    real(dp), intent(in)  :: r_cell
    real(dp), intent(out) :: r_host
    integer,  intent(out) :: capac
    logical,  intent(out) :: found
    real(dp) :: tot, u, cum
    integer  :: i, ilast

    found = .false.; r_host = 0.0_dp; capac = 0
    tot = 0.0_dp; ilast = 0
    do i = 1, NBIN
      if (r_bin(i) >= r_cell) then
        tot = tot + spec(i,k)
        if (spec(i,k) > 0.0_dp) ilast = i
      end if
    end do
    if (tot <= 0.0_dp) return

    call random_number(u)
    u = u * tot
    cum = 0.0_dp
    do i = 1, NBIN
      if (r_bin(i) >= r_cell) then
        cum = cum + spec(i,k)
        if (u <= cum .and. spec(i,k) > 0.0_dp) then
          r_host = r_bin(i)
          capac  = host_capacity(r_host, r_cell)
          found  = .true.
          return
        end if
      end if
    end do
    ! Floating-point fallback: u landed just past the running sum -> last bin.
    r_host = r_bin(ilast)
    capac  = host_capacity(r_host, r_cell)
    found  = .true.
  end subroutine draw_host

  ! Verification driver: sweep r_cell at the peak-spectrum level and report
  ! eligible host supply, abundance-weighted <r_host>, mean and max capacity.
  ! A Monte-Carlo check confirms draw_host reproduces the analytic weighting.
  subroutine host_capacity_sweep(unit)
    integer, intent(in) :: unit
    integer,  parameter :: NRC = 7, NDRAW = 200000
    real(dp), parameter :: rc_um(NRC) = &
         [0.2_dp, 0.3_dp, 0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    ! maxC ignores the negligible renormalised-lognormal tail bins (same
    ! abundance floor as the spectrum dump): only count bins with real droplets.
    real(dp), parameter :: SPEC_FLOOR = 1.0e-2_dp   ! [m^-3]
    real(dp) :: Npk, supply, wmean_r, meanC, rc, contrib, mc_r, dmc
    integer  :: i, k, kpk, j, maxC, capj, ndr
    real(dp) :: rh
    logical  :: ok

    ! Peak total-spectrum level (same convention as cloud_spectrum_check).
    Npk = 0.0_dp; kpk = 1
    do k = 1, nz
      if (sum(spec(:,k)) > Npk) then; Npk = sum(spec(:,k)); kpk = k; end if
    end do

    write(unit,'(a)') ' '
    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a)') ' HOST-CAPACITY ENGINE  C(r_cell) = floor(phi*(r_host/r_cell)^3)'
    write(unit,'(a,f4.2,a,f6.2,a)') '  phi = ', phi_pack, &
         ';  peak-spectrum level z = ', z(kpk)/1000.0_dp, ' km'
    write(unit,'(a)') '  r_cell[um]  supply[m-3]   <r_host>[um]  meanC   maxC    MC<r_host>  dMC'
    do i = 1, NRC
      rc = rc_um(i) * 1.0e-6_dp
      supply  = host_supply(kpk, rc)
      if (supply <= 0.0_dp) then
        write(unit,'(2x,f8.2,4x,a)') rc_um(i), 'no eligible host'
        cycle
      end if
      ! Analytic abundance-weighted mean radius and mean capacity over bins.
      wmean_r = 0.0_dp; meanC = 0.0_dp; maxC = 0
      do j = 1, NBIN
        if (r_bin(j) >= rc .and. spec(j,kpk) > 0.0_dp) then
          contrib = spec(j,kpk) / supply
          wmean_r = wmean_r + contrib * r_bin(j)
          capj    = host_capacity(r_bin(j), rc)
          meanC   = meanC + contrib * real(capj, dp)
          if (spec(j,kpk) > SPEC_FLOOR) maxC = max(maxC, capj)
        end if
      end do
      ! Monte-Carlo consistency check on the sampler.
      mc_r = 0.0_dp; ndr = 0
      do j = 1, NDRAW
        call draw_host(kpk, rc, rh, capj, ok)
        if (ok) then; mc_r = mc_r + rh; ndr = ndr + 1; end if
      end do
      if (ndr > 0) mc_r = mc_r / real(ndr, dp)
      dmc = mc_r / wmean_r - 1.0_dp
      write(unit,'(2x,f8.2,2x,es12.4,2x,f9.4,2x,f8.2,2x,i6,4x,f9.4,2x,es9.2)') &
           rc_um(i), supply, wmean_r*1.0e6_dp, meanC, maxC, mc_r*1.0e6_dp, dmc
    end do
    write(unit,'(a)') '  (meanC,maxC near 0-1 => capacity-starved; >>1 => colonies form)'
    write(unit,'(a)') '----------------------------------------------------------'
  end subroutine host_capacity_sweep

  !====================================================================
  ! Settling physics — Stokes drag + Cunningham slip
  !   Terminal fall speed of a single sphere, matched to the per-particle
  !   physics of VPCM cloudvenus/new_cloud_sedim.F (so the 1-D model and a
  !   future 3-D coupling use the same drag law).  Pure functions of the
  !   static column.
  !====================================================================

  ! Dynamic viscosity of CO2 [Pa.s], Johnston & Grilly (1942) via Jones/
  ! Lennard-Jones (VPCM VISCOSITY_CO2; valid ~80-300 K).
  real(dp) function viscosity_co2(temp) result(mu)
    real(dp), intent(in) :: temp
    real(dp) :: numer, denom
    numer = 200.0_dp**(2.27_dp/4.27_dp) - 0.435_dp
    denom = temp **(2.27_dp/4.27_dp) - 0.435_dp
    mu = (numer/denom) * 1015.0_dp * (temp/200.0_dp)**1.5_dp
    mu = mu * 1.0e-8_dp                          ! Poise*1e7 -> Pa.s
  end function viscosity_co2

  ! Gas mean free path [m] (VPCM new_cloud_sedim form):
  !   lambda = (T/P) * 0.707 * R / (4*pi*molrad^2*N_A)   [= k_B T / (sqrt2 pi d^2 P)]
  real(dp) function mean_free_path(temp, pres) result(lam)
    real(dp), intent(in) :: temp, pres
    lam = (temp/pres) * (0.707_dp * R_UNIV / (4.0_dp*PI*CO2_MOLRAD**2 * N_AVO))
  end function mean_free_path

  ! Cunningham slip correction for a sphere radius r in gas of mean free path
  ! lambda.  Cc = 1 + Kn*(A + B*exp(-C/Kn)),  Kn = lambda/r.  >= 1.
  real(dp) function cunningham(r, lambda) result(Cc)
    real(dp), intent(in) :: r, lambda
    real(dp) :: kn
    Cc = 1.0_dp
    if (r <= 0.0_dp) return
    kn = lambda / r
    Cc = 1.0_dp + kn*(SLIP_A + SLIP_B*exp(-SLIP_C/kn))
  end function cunningham

  ! Core terminal settling velocity [m/s], POSITIVE DOWNWARD, from local
  ! atmosphere state (Tq,Pq,rho_a) for a sphere of radius r and density rho_p.
  real(dp) function settling_core(Tq, Pq, rho_a, r, rho_p) result(vset)
    real(dp), intent(in) :: Tq, Pq, rho_a, r, rho_p
    real(dp) :: mu, lam, Cc
    vset = 0.0_dp
    if (r <= 0.0_dp) return
    mu  = viscosity_co2(Tq)
    lam = mean_free_path(Tq, Pq)
    Cc  = cunningham(r, lam)
    vset = (2.0_dp/9.0_dp) * (rho_p - rho_a) * GRAV_VENUS * r*r / mu * Cc
  end function settling_core

  ! Settling at grid level k (uses the ingested VPCM gas density).
  ! Net vertical motion in the lifecycle is w(k) - settling_velocity(...).
  real(dp) function settling_velocity(k, r, rho_p) result(vset)
    integer,  intent(in) :: k
    real(dp), intent(in) :: r, rho_p
    vset = settling_core(T(k), P(k), rho_air(k), r, rho_p)
  end function settling_velocity

  ! Settling at an arbitrary altitude zq [m] (atmosphere fields interpolated).
  real(dp) function settling_velocity_z(zq, r, rho_p) result(vset)
    real(dp), intent(in) :: zq, r, rho_p
    vset = settling_core(col_interp(T, zq), col_interp(P, zq), &
                         col_interp(rho_air, zq), r, rho_p)
  end function settling_velocity_z

  !====================================================================
  ! Column interpolation helpers (organisms live at continuous z)
  !====================================================================

  ! Lower-bracket index klo with z(klo) <= zq < z(klo+1); clamped to [1,nz-1].
  ! Assumes z(:) ascending (VPCM master grid).
  integer function col_locate(zq) result(klo)
    real(dp), intent(in) :: zq
    integer :: ihi, imid
    if (zq <= z(1))  then; klo = 1;    return; end if
    if (zq >= z(nz)) then; klo = nz-1; return; end if
    klo = 1; ihi = nz
    do while (ihi - klo > 1)
      imid = (klo + ihi) / 2
      if (zq >= z(imid)) then; klo = imid; else; ihi = imid; end if
    end do
  end function col_locate

  ! Linear interpolation of a column field f(:) at altitude zq [m] (clamped).
  real(dp) function col_interp(f, zq) result(val)
    real(dp), intent(in) :: f(:), zq
    integer  :: k
    real(dp) :: t
    k   = col_locate(zq)
    t   = (zq - z(k)) / (z(k+1) - z(k))
    t   = max(0.0_dp, min(1.0_dp, t))
    val = f(k) + t * (f(k+1) - f(k))
  end function col_interp

  ! Local vertical gradient of Kzz [1/s] at zq (for the random-walk drift term).
  real(dp) function dKzz_dz(zq) result(dk)
    real(dp), intent(in) :: zq
    integer :: k
    k  = col_locate(zq)
    dk = (Kzz(k+1) - Kzz(k)) / (z(k+1) - z(k))
  end function dKzz_dz

  !====================================================================
  subroutine cloud_cleanup()
    if (allocated(z))        deallocate(z)
    if (allocated(T))        deallocate(T)
    if (allocated(P))        deallocate(P)
    if (allocated(rho_air))  deallocate(rho_air)
    if (allocated(w))        deallocate(w)
    if (allocated(Kzz))      deallocate(Kzz)
    if (allocated(WSA))      deallocate(WSA)
    if (allocated(rho_host)) deallocate(rho_host)
    if (allocated(vN1))      deallocate(vN1)
    if (allocated(vr1))      deallocate(vr1)
    if (allocated(vN2))      deallocate(vN2)
    if (allocated(vr2))      deallocate(vr2)
    if (allocated(vvalid1))  deallocate(vvalid1)
    if (allocated(vvalid2))  deallocate(vvalid2)
    if (allocated(Nm))       deallocate(Nm)
    if (allocated(rm))       deallocate(rm)
    if (allocated(sm))       deallocate(sm)
    if (allocated(validm))   deallocate(validm)
    if (allocated(r_edge))   deallocate(r_edge)
    if (allocated(r_bin))    deallocate(r_bin)
    if (allocated(specm))    deallocate(specm)
    if (allocated(spec))     deallocate(spec)
  end subroutine cloud_cleanup

end module bio_cloud_model
