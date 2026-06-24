!======================================================================
! bio_cloud_model.F90
! Venus cloud biosphere model (Seager-cycle carrying-capacity diagnostic)
!
! STEP 1 (this file, for now): the input layer only.
!   - cloud_read       : ingest the VPCM column profiles from inputs/
!   - cloud_diagnostics: write a per-level table to a log
!   - cloud_summary    : write summary statistics (cloud window, peaks)
! A thin test program (bio_cloud_test, bottom of file) drives these so we
! can verify the columns are read correctly before any biology is added.
!
! Input files (aerial_biosphere/inputs/vpcm_*.txt): single VPCM column,
! 2-column "altitude[km]  value" with a commented (#) header, common grid.
!
! Units: SI throughout (m, kg, s, K).  Altitude is converted km -> m on read.
!
! Host data hygiene (see memory project-venus-inputs):
!   * r1/r2 carry sentinel/fill garbage outside the cloud: 0 below cloud,
!     a constant 1.143147 m fill above it, and spurious huge values at the
!     bottom edge where droplet number ~ 0.  A level's host is only trusted
!     where droplet number exceeds N_DROP_MIN and the radius is physical.
!   * The valid cloud window then falls out at ~49-82 km.
!   * Acidity (WSA) is read but is NOT a habitability criterion in this
!     project; temperature is the sole AHZ axis.  WSA/rho_host are kept for
!     host density / diagnostics only.
!======================================================================
module bio_cloud_model
  implicit none
  private

  integer, parameter :: dp = 8

  ! ---- Cloud column grid (filled by cloud_read) ----
  integer,  save              :: nz = 0
  real(dp), allocatable, save :: z(:)        ! altitude        [m]
  real(dp), allocatable, save :: T(:)        ! temperature     [K]
  real(dp), allocatable, save :: P(:)        ! pressure        [Pa]
  real(dp), allocatable, save :: rho_air(:)  ! gas density     [kg/m^3]
  real(dp), allocatable, save :: w(:)        ! vertical wind   [m/s]
  real(dp), allocatable, save :: Kzz(:)      ! eddy diffusivity[m^2/s]
  real(dp), allocatable, save :: N1(:)       ! mode-1 droplet number [m^-3]
  real(dp), allocatable, save :: r1(:)       ! mode-1 mean radius    [m]
  real(dp), allocatable, save :: N2(:)       ! mode-2 droplet number [m^-3]
  real(dp), allocatable, save :: r2(:)       ! mode-2 mean radius    [m]
  real(dp), allocatable, save :: WSA(:)      ! H2SO4 mass fraction   [-]
  real(dp), allocatable, save :: rho_host(:) ! droplet density       [kg/m^3]
  logical,  allocatable, save :: valid1(:)   ! mode-1 host trustworthy at level
  logical,  allocatable, save :: valid2(:)   ! mode-2 host trustworthy at level

  ! ---- Reconstructed binned host spectrum (filled by cloud_build_spectrum) ----
  ! Each mode is a log-normal of fixed width (SIGMA1/SIGMA2) and median radius
  ! r1/r2 (the VPCM rmean == geometric median r0, verified via alpha_k).  We bin
  ! it onto a shared log-spaced radius grid, conserving droplet number exactly.
  integer,  parameter, public :: NBIN = 80
  real(dp), parameter         :: R_GRID_MIN = 1.0e-8_dp   ! grid lower edge [m] (0.01 um)
  real(dp), parameter         :: R_GRID_MAX = 1.0e-5_dp   ! grid upper edge [m] (10 um)
  real(dp), allocatable, save :: r_edge(:)   ! bin edges   [m]  (NBIN+1)
  real(dp), allocatable, save :: r_bin(:)    ! bin centres [m]  (NBIN)
  real(dp), allocatable, save :: spec1(:,:)  ! mode-1 number/bin [m^-3]  (NBIN, nz)
  real(dp), allocatable, save :: spec2(:,:)  ! mode-2 number/bin [m^-3]  (NBIN, nz)
  real(dp), allocatable, save :: spec(:,:)   ! total host spectrum [m^-3] (NBIN, nz)

  ! ---- Fixed log-normal widths of the VPCM modes (conf_phys.F90 defaults) ----
  real(dp), parameter, public :: SIGMA1 = 1.56_dp   ! mode 1
  real(dp), parameter, public :: SIGMA2 = 1.29_dp   ! mode 2

  ! ---- Host data-hygiene masking thresholds ----
  real(dp), parameter, public :: N_DROP_MIN = 1.0e3_dp   ! min trustworthy droplet number [m^-3]
  ! Per-mode upper radius guards.  A mode's mean radius should sit near its
  ! geometric radius (ri1=0.33 um, ri2=1.0 um) and grow only modestly in the
  ! cloud deck (mode1 stays sub-um, mode2 ~1 um).  Radii above these caps are
  ! outliers: the constant 1.143147 m fill above the cloud, and transitional
  ! cloud-base levels where r=(M3/M0)^(1/3) spikes (e.g. mode1 ~6 um at 52 km
  ! while the mode1 norm is 0.05-0.43 um).  These caps subsume the old fill
  ! rejection (both << 1 m), so a separate physical cap is no longer needed.
  real(dp), parameter, public :: R1_MAX = 1.0e-6_dp   ! mode-1 max trusted mean radius [m] (~3 x ri1)
  real(dp), parameter, public :: R2_MAX = 3.0e-6_dp   ! mode-2 max trusted mean radius [m] (~3 x ri2)

  public :: cloud_read, cloud_diagnostics, cloud_summary, cloud_cleanup
  public :: cloud_build_spectrum, cloud_spectrum_check
  public :: nz, z, T, P, rho_air, w, Kzz, N1, r1, N2, r2, WSA, rho_host
  public :: valid1, valid2
  public :: r_bin, r_edge, spec1, spec2, spec

contains

  !--------------------------------------------------------------------
  ! Count the data (non-comment, non-blank) rows in a column file.
  !--------------------------------------------------------------------
  integer function count_rows(path) result(n)
    character(*), intent(in) :: path
    integer :: u, ios
    character(len=512) :: line
    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(a)') 'ERROR: cannot open '//trim(path)
      stop 1
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

  !--------------------------------------------------------------------
  ! Read a 2-column "altitude[km]  value" file into zc (km) and vc.
  ! Arrays must be pre-sized to the expected number of rows.
  !--------------------------------------------------------------------
  subroutine read_column(path, zc, vc)
    character(*), intent(in)  :: path
    real(dp),     intent(out) :: zc(:), vc(:)
    integer :: u, ios, k
    character(len=512) :: line
    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(a)') 'ERROR: cannot open '//trim(path)
      stop 1
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
        write(*,'(a,i0)') 'ERROR: parse failure in '//trim(path)//' at data row ', k
        stop 1
      end if
    end do
    close(u)
    if (k /= size(zc)) then
      write(*,'(a,i0,a,i0)') 'ERROR: '//trim(path)//' row count ', k, &
                             ' /= expected ', size(zc)
      stop 1
    end if
  end subroutine read_column

  !--------------------------------------------------------------------
  ! Read a field, asserting its altitude grid matches the master grid.
  !--------------------------------------------------------------------
  subroutine read_field(indir, fname, dest, zkm_ref)
    character(*), intent(in)  :: indir, fname
    real(dp),     intent(out) :: dest(:)
    real(dp),     intent(in)  :: zkm_ref(:)
    real(dp) :: zkm(size(dest))
    real(dp) :: dzmax
    call read_column(trim(indir)//'/'//trim(fname), zkm, dest)
    dzmax = maxval(abs(zkm - zkm_ref))
    if (dzmax > 1.0e-3_dp) then
      write(*,'(a,es10.3)') 'ERROR: '//trim(fname)// &
        ' altitude grid mismatch, max |dz| [km] = ', dzmax
      stop 1
    end if
  end subroutine read_field

  !--------------------------------------------------------------------
  ! Ingest all VPCM column profiles from directory `indir`.
  !--------------------------------------------------------------------
  subroutine cloud_read(indir)
    character(*), intent(in) :: indir
    real(dp), allocatable :: zkm(:)
    integer :: k

    ! Number of levels is set by the temperature file; all share the grid.
    nz = count_rows(trim(indir)//'/vpcm_temp.txt')
    if (nz < 2) then
      write(*,'(a)') 'ERROR: fewer than 2 data rows found in vpcm_temp.txt'
      stop 1
    end if

    call cloud_cleanup()
    allocate(z(nz), T(nz), P(nz), rho_air(nz), w(nz), Kzz(nz))
    allocate(N1(nz), r1(nz), N2(nz), r2(nz), WSA(nz), rho_host(nz))
    allocate(valid1(nz), valid2(nz))
    allocate(zkm(nz))

    ! Master grid from temperature file (km in file -> m internally).
    call read_column(trim(indir)//'/vpcm_temp.txt', zkm, T)
    z = zkm * 1000.0_dp

    ! Remaining fields, each checked against the master grid.
    call read_field(indir, 'vpcm_pres.txt',         P,        zkm)
    call read_field(indir, 'vpcm_rho.txt',          rho_air,  zkm)
    call read_field(indir, 'vpcm_vitwz.txt',        w,        zkm)
    call read_field(indir, 'vpcm_kz.txt',           Kzz,      zkm)
    call read_field(indir, 'vpcm_m0_mode1drop.txt', N1,       zkm)
    call read_field(indir, 'vpcm_r1.txt',           r1,       zkm)
    call read_field(indir, 'vpcm_m0_mode2drop.txt', N2,       zkm)
    call read_field(indir, 'vpcm_r2.txt',           r2,       zkm)
    call read_field(indir, 'vpcm_WSAVtab.txt',      WSA,      zkm)
    call read_field(indir, 'vpcm_rhotab.txt',       rho_host, zkm)

    ! Host data-hygiene mask: trust a mode's host only where droplet number
    ! is substantial AND the mean radius is positive and below the per-mode
    ! outlier guard.  This rejects (i) the 0 fill below cloud, (ii) the
    ! 1.143147 m fill above cloud, (iii) the near-zero-number bottom-edge
    ! radius blow-ups, and (iv) substantial-number but anomalously large
    ! cloud-base radii (the ~6 um mode-1 spike at 52 km).
    do k = 1, nz
      valid1(k) = (N1(k) > N_DROP_MIN) .and. (r1(k) > 0.0_dp) .and. (r1(k) < R1_MAX)
      valid2(k) = (N2(k) > N_DROP_MIN) .and. (r2(k) > 0.0_dp) .and. (r2(k) < R2_MAX)
    end do

    deallocate(zkm)
  end subroutine cloud_read

  !--------------------------------------------------------------------
  ! Write the full per-level diagnostic table to `unit`.
  !--------------------------------------------------------------------
  subroutine cloud_diagnostics(unit)
    integer, intent(in) :: unit
    integer :: k
    write(unit,'(a)') '# Venus cloud column — ingested VPCM profiles'
    write(unit,'(a,i0,a)') '# levels = ', nz, &
      '   (radii reported in microns; v1/v2 = host-trusted flag)'
    write(unit,'(a)') '#'
    write(unit,'(a)') &
      '#   k   z[km]    T[K]       P[Pa]    rho_air[kg/m3]   w[m/s]      Kzz[m2/s]'// &
      '   N1[m-3]     r1[um]  v1    N2[m-3]     r2[um]  v2    WSA   rho_h[kg/m3]'
    do k = 1, nz
      write(unit,'(i5,1x,f8.3,1x,f9.2,3x,es10.3,3x,es11.3,3x,es11.3,3x,es10.3,&
                   &3x,es10.3,1x,es10.3,1x,i2,3x,es10.3,1x,es10.3,1x,i2,3x,f5.3,3x,f7.1)') &
        k, z(k)/1000.0_dp, T(k), P(k), rho_air(k), w(k), Kzz(k), &
        N1(k), r1(k)*1.0e6_dp, merge(1,0,valid1(k)), &
        N2(k), r2(k)*1.0e6_dp, merge(1,0,valid2(k)), WSA(k), rho_host(k)
    end do
  end subroutine cloud_diagnostics

  !--------------------------------------------------------------------
  ! Write summary statistics (valid window, peaks, ranges) to `unit`.
  !--------------------------------------------------------------------
  subroutine cloud_summary(unit)
    integer, intent(in) :: unit
    logical :: anyhost(nz)
    integer :: kb, ktop, i1, i2
    integer :: nv1, nv2, ng1, ng2

    anyhost = valid1 .or. valid2
    nv1 = count(valid1)
    nv2 = count(valid2)
    ! Levels with substantial droplet number but rejected by the radius guard.
    ng1 = count( (N1 > N_DROP_MIN) .and. (r1 > 0.0_dp) .and. (r1 >= R1_MAX) )
    ng2 = count( (N2 > N_DROP_MIN) .and. (r2 > 0.0_dp) .and. (r2 >= R2_MAX) )

    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a)') ' CLOUD COLUMN SUMMARY'
    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a,i0)')        '  levels read           : ', nz
    write(unit,'(a,f7.2,a,f7.2,a)') '  altitude range        : ', &
         z(1)/1000.0_dp, ' to ', z(nz)/1000.0_dp, ' km'
    write(unit,'(a,f7.1,a,f7.1,a)') '  temperature range     : ', &
         minval(T), ' to ', maxval(T), ' K'
    write(unit,'(a,i0,a,i0)')   '  trusted host levels   : mode1 = ', nv1, &
         '   mode2 = ', nv2
    write(unit,'(a,f5.2,a,f5.2,a)') '  radius outlier guards : mode1 < ', &
         R1_MAX*1.0e6_dp, ' um   mode2 < ', R2_MAX*1.0e6_dp, ' um'
    write(unit,'(a,i0,a,i0)')   '  guard-rejected levels : mode1 = ', ng1, &
         '   mode2 = ', ng2

    if (count(anyhost) > 0) then
      kb   = minloc(z, dim=1, mask=anyhost)
      ktop = maxloc(z, dim=1, mask=anyhost)
      write(unit,'(a,f7.2,a,f7.2,a)') '  cloud (host) window   : ', &
           z(kb)/1000.0_dp, ' to ', z(ktop)/1000.0_dp, ' km'
    else
      write(unit,'(a)') '  cloud (host) window   : NONE (no trusted host levels!)'
    end if

    if (nv1 > 0) then
      i1 = maxloc(N1, dim=1, mask=valid1)
      write(unit,'(a,es10.3,a,f6.2,a)') '  peak N1               : ', &
           N1(i1), ' m-3 at ', z(i1)/1000.0_dp, ' km'
      write(unit,'(a,f7.4,a,f7.4,a)') '  mode-1 radius range   : ', &
           minval(r1, mask=valid1)*1.0e6_dp, ' to ', &
           maxval(r1, mask=valid1)*1.0e6_dp, ' um'
    end if
    if (nv2 > 0) then
      i2 = maxloc(N2, dim=1, mask=valid2)
      write(unit,'(a,es10.3,a,f6.2,a)') '  peak N2               : ', &
           N2(i2), ' m-3 at ', z(i2)/1000.0_dp, ' km'
      write(unit,'(a,f7.4,a,f7.4,a)') '  mode-2 radius range   : ', &
           minval(r2, mask=valid2)*1.0e6_dp, ' to ', &
           maxval(r2, mask=valid2)*1.0e6_dp, ' um'
    end if
    if (count(anyhost) > 0) then
      write(unit,'(a,f5.3,a,f5.3,a)') '  WSA range (in cloud)  : ', &
           minval(WSA, mask=anyhost), ' to ', maxval(WSA, mask=anyhost), &
           '  (diagnostic only; not a habitability axis)'
      write(unit,'(a,f7.1,a,f7.1,a)') '  host density (in cld) : ', &
           minval(rho_host, mask=anyhost), ' to ', &
           maxval(rho_host, mask=anyhost), ' kg/m3'
    end if
    write(unit,'(a)') '----------------------------------------------------------'
  end subroutine cloud_summary

  !--------------------------------------------------------------------
  ! Distribute Ntot droplets of a log-normal mode (median radius r0,
  ! geometric width sigma_g) across the radius grid, using the exact
  ! log-normal CDF (erf) over each bin.  Renormalised so the binned
  ! counts sum to Ntot exactly over the finite grid.
  !--------------------------------------------------------------------
  subroutine lognormal_bins(Ntot, r0, sigma_g, out)
    real(dp), intent(in)  :: Ntot, r0, sigma_g
    real(dp), intent(out) :: out(:)
    real(dp) :: s, sq2, denom, lo, hi
    integer  :: i
    out = 0.0_dp
    if (Ntot <= 0.0_dp .or. r0 <= 0.0_dp) return
    s   = log(sigma_g)
    sq2 = s * sqrt(2.0_dp)
    ! Total probability captured by the finite grid (for exact renormalisation).
    denom = 0.5_dp * ( erf(log(R_GRID_MAX/r0)/sq2) - erf(log(R_GRID_MIN/r0)/sq2) )
    do i = 1, NBIN
      lo = 0.5_dp * erf(log(r_edge(i)  /r0)/sq2)
      hi = 0.5_dp * erf(log(r_edge(i+1)/r0)/sq2)
      out(i) = Ntot * (hi - lo) / denom
    end do
  end subroutine lognormal_bins

  !--------------------------------------------------------------------
  ! Build the shared radius grid and reconstruct the per-level binned
  ! host spectrum for each mode (only at trusted levels).
  !--------------------------------------------------------------------
  subroutine cloud_build_spectrum()
    integer  :: i, k
    real(dp) :: ratio

    if (nz <= 0) then
      write(*,'(a)') 'ERROR: cloud_build_spectrum called before cloud_read'
      stop 1
    end if

    if (allocated(r_edge)) deallocate(r_edge)
    if (allocated(r_bin))  deallocate(r_bin)
    if (allocated(spec1))  deallocate(spec1)
    if (allocated(spec2))  deallocate(spec2)
    if (allocated(spec))   deallocate(spec)
    allocate(r_edge(NBIN+1), r_bin(NBIN))
    allocate(spec1(NBIN, nz), spec2(NBIN, nz), spec(NBIN, nz))

    ! Log-spaced grid: edges geometrically spaced, centres as geometric means.
    ratio = (R_GRID_MAX / R_GRID_MIN) ** (1.0_dp / real(NBIN, dp))
    do i = 1, NBIN+1
      r_edge(i) = R_GRID_MIN * ratio ** real(i-1, dp)
    end do
    do i = 1, NBIN
      r_bin(i) = sqrt(r_edge(i) * r_edge(i+1))
    end do

    spec1 = 0.0_dp
    spec2 = 0.0_dp
    do k = 1, nz
      if (valid1(k)) call lognormal_bins(N1(k), r1(k), SIGMA1, spec1(:,k))
      if (valid2(k)) call lognormal_bins(N2(k), r2(k), SIGMA2, spec2(:,k))
      spec(:,k) = spec1(:,k) + spec2(:,k)
    end do
  end subroutine cloud_build_spectrum

  !--------------------------------------------------------------------
  ! Verify the reconstruction: per trusted level/mode, the binned counts
  ! should recover the input droplet number and geometric-mean radius.
  ! Writes per-level errors plus a spectrum dump at the peak-N1 level.
  !--------------------------------------------------------------------
  subroutine cloud_spectrum_check(unit)
    integer, intent(in) :: unit
    integer  :: k, i, kpk
    real(dp) :: nrec, rg, errN, errR, maxN, maxR
    real(dp) :: sumn, sumln

    write(unit,'(a)') ' '
    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a)') ' HOST SPECTRUM RECONSTRUCTION CHECK'
    write(unit,'(a,i0,a,f6.3,a,f6.3,a)') '  grid: ', NBIN, ' log bins, ', &
         r_edge(1)*1.0e6_dp, ' to ', r_edge(NBIN+1)*1.0e6_dp, ' um'
    write(unit,'(a)') '  (errN = recovered/input droplet number - 1;'// &
                      '  errR = recovered/input geom-mean radius - 1)'
    write(unit,'(a)') '#   z[km]  mode   N_in[m-3]   N_rec[m-3]   errN'// &
                      '      r0_in[um]  rg_rec[um]   errR'

    maxN = 0.0_dp
    maxR = 0.0_dp
    do k = 1, nz
      call check_one(k, 1, valid1(k), N1(k), r1(k), spec1(:,k), maxN, maxR, unit)
      call check_one(k, 2, valid2(k), N2(k), r2(k), spec2(:,k), maxN, maxR, unit)
    end do
    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a,es9.2,a,es9.2)') '  max |errN| = ', maxN, '   max |errR| = ', maxR
    if (maxN < 1.0e-3_dp .and. maxR < 1.0e-2_dp) then
      write(unit,'(a)') '  RECONSTRUCTION OK (number conserved; radius recovered)'
    else
      write(unit,'(a)') '  WARNING: reconstruction error larger than expected'
    end if

    ! Spectrum dump at the peak-N1 level.
    if (count(valid1) > 0) then
      kpk = maxloc(N1, dim=1, mask=valid1)
      write(unit,'(a)') ' '
      write(unit,'(a,f6.2,a)') '  combined host spectrum at peak-N1 level (', &
           z(kpk)/1000.0_dp, ' km), non-empty bins:'
      write(unit,'(a)') '#    r[um]        n[m-3]'
      do i = 1, NBIN
        if (spec(i,kpk) > 1.0e-3_dp) &
          write(unit,'(2x,f10.4,2x,es12.4)') r_bin(i)*1.0e6_dp, spec(i,kpk)
      end do
    end if
    write(unit,'(a)') '----------------------------------------------------------'
    return

  contains
    subroutine check_one(kk, mode, ok, Nin, r0in, sp, mN, mR, u)
      integer,  intent(in)    :: kk, mode, u
      logical,  intent(in)    :: ok
      real(dp), intent(in)    :: Nin, r0in, sp(:)
      real(dp), intent(inout) :: mN, mR
      if (.not. ok) return
      nrec  = sum(sp)
      sumn  = nrec
      sumln = sum(sp * log(r_bin))
      rg    = exp(sumln / sumn)          ! recovered geometric-mean radius
      errN  = nrec / Nin - 1.0_dp
      errR  = rg   / r0in - 1.0_dp
      mN = max(mN, abs(errN))
      mR = max(mR, abs(errR))
      write(u,'(2x,f7.2,3x,i1,3x,es11.3,1x,es11.3,1x,f8.4,3x,f8.4,3x,f8.4,3x,f8.4)') &
           z(kk)/1000.0_dp, mode, Nin, nrec, errN, r0in*1.0e6_dp, rg*1.0e6_dp, errR
    end subroutine check_one
  end subroutine cloud_spectrum_check

  !--------------------------------------------------------------------
  subroutine cloud_cleanup()
    if (allocated(z))        deallocate(z)
    if (allocated(T))        deallocate(T)
    if (allocated(P))        deallocate(P)
    if (allocated(rho_air))  deallocate(rho_air)
    if (allocated(w))        deallocate(w)
    if (allocated(Kzz))      deallocate(Kzz)
    if (allocated(N1))       deallocate(N1)
    if (allocated(r1))       deallocate(r1)
    if (allocated(N2))       deallocate(N2)
    if (allocated(r2))       deallocate(r2)
    if (allocated(WSA))      deallocate(WSA)
    if (allocated(rho_host)) deallocate(rho_host)
    if (allocated(valid1))   deallocate(valid1)
    if (allocated(valid2))   deallocate(valid2)
    if (allocated(r_edge))   deallocate(r_edge)
    if (allocated(r_bin))    deallocate(r_bin)
    if (allocated(spec1))    deallocate(spec1)
    if (allocated(spec2))    deallocate(spec2)
    if (allocated(spec))     deallocate(spec)
  end subroutine cloud_cleanup

end module bio_cloud_model


!======================================================================
! Temporary test driver: read the columns, dump diagnostics to a log,
! and echo the summary to stdout.
!   ./bio_cloud              (reads ../inputs)
!   ./bio_cloud <input_dir>
!======================================================================
program bio_cloud_test
  use bio_cloud_model
  implicit none
  character(len=256) :: indir
  character(len=256) :: logfile
  integer :: ulog

  indir   = '../inputs'
  logfile = 'cloud_read.log'
  if (command_argument_count() >= 1) call get_command_argument(1, indir)

  write(*,'(a)') '============================================'
  write(*,'(a)') ' Venus cloud model — input reader (step 1)'
  write(*,'(a)') '============================================'
  write(*,'(a)') ' input dir : '//trim(indir)
  write(*,'(a)') ' log file  : '//trim(logfile)

  call cloud_read(trim(indir))
  call cloud_build_spectrum()

  open(newunit=ulog, file=trim(logfile), status='replace', action='write')
  call cloud_diagnostics(ulog)
  call cloud_summary(ulog)
  call cloud_spectrum_check(ulog)
  close(ulog)

  ! Echo the summary + reconstruction check to the terminal too.
  call cloud_summary(6)
  call cloud_spectrum_check(6)

  write(*,'(a)') ' Done.  Full per-level table in '//trim(logfile)//'.'

  call cloud_cleanup()

end program bio_cloud_test
