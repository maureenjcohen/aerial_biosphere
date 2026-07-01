!======================================================================
! bio_venus.F90
! Venus cloud biosphere — BIOLOGY MODULE.
!
! The organisms and their dynamics, layered on the static environment in
! module bio_cloud_model (cloud spectrum, settling, column interpolation):
!   - structure-of-arrays organism store + 3-state machine (per plan s.2.2)
!   - Lagrangian transport of free spores (advection + eddy-diffusion walk)
!   - verification drivers: settling_check, transport_test
!
! Built incrementally: STEP 1 (here) = organism store + transport only
! (states DEAD/DORMANT).  Germination/ACTIVE/DEPOT, fission and mutation
! arrive in later steps.
!
! Units: SI throughout (m, kg, s, K).
!======================================================================
module bio_venus
  use bio_cloud_model
  implicit none
  private

  integer, parameter :: dp = 8

  ! ---- Fixed desiccated-spore density [kg/m^3] (a knob; plan section 6) ----
  real(dp), parameter, public :: RHO_SPORE_DEFAULT = 1200.0_dp
  real(dp), save,      public :: rho_spore = RHO_SPORE_DEFAULT
  ! ---- Hydrated (in-droplet, actively-dividing) cell density [kg/m^3] ----
  ! Distinct from the DESICCATED spore density above: a live cell bathed in liquid
  ! is near water (~1000-1100).  Used only for the combined-droplet settling
  ! verdict in the depot return-leg test (rho_spore governs free-spore settling).
  real(dp), save, public :: rho_cell_wet = 1050.0_dp

  ! ---- Organism store (structure-of-arrays) + 3-state machine ----
  ! State machine per plan section 2.2.  Step-1 transport uses DEAD/DORMANT
  ! only; ACTIVE/DEPOT transitions arrive with germination (step 2).
  integer, parameter, public :: ST_DEAD = 0, ST_ACTIVE = 1, ST_DORMANT = 2, ST_DEPOT = 3
  integer,  save              :: max_orgs = 0      ! allocated slots
  integer,  save              :: n_orgs   = 0      ! slots in use (high-water mark)
  integer,  allocatable, save :: o_state(:)        ! ST_*
  real(dp), allocatable, save :: o_z(:)            ! altitude [m]
  real(dp), allocatable, save :: o_rcell(:)        ! cell radius [m] (heritable; fixed pre-step-3)
  real(dp), allocatable, save :: o_age(:)          ! time as DORMANT spore (the Y-day clock) [s]
  real(dp), allocatable, save :: o_rhost(:)        ! host droplet radius [m] (0 unless ACTIVE)
  integer,  allocatable, save :: o_capac(:)        ! host packing capacity
  real(dp), allocatable, save :: o_X(:)            ! reproduction half-life [s] (heritable)
  integer,  allocatable, save :: o_host(:)         ! host registry index (0 unless ACTIVE)
  integer,  allocatable, save :: o_birthstep(:)    ! step created (same-step division guard)
  real(dp), allocatable, save :: o_uconv(:)        ! per-parcel OU convective state [-] (free spores)
  integer,  save              :: n_lost_top = 0, n_lost_bot = 0   ! boundary-loss tallies
  integer,  save              :: n_germ = 0, n_dead_dormant = 0   ! germinations / Y-clock deaths
  integer,  save              :: n_birth = 0, n_eject = 0         ! surviving divisions / ejected daughters (died)
  integer,  save              :: n_dead_cell = 0, n_host_die = 0  ! ACTIVE-cell deaths / host-droplet deaths
  ! Depot return-leg test: colony deaths whose host droplet was heavy enough to
  ! sediment to the depot (rain-out, cells -> DEPOT) vs desiccated in place (cells
  ! -> DORMANT in cloud), and cells delivered to the depot bank.
  integer,  save              :: n_rainout = 0, n_desicc_place = 0, n_cells_depot = 0
  real(dp), save              :: rmin_settle = 0.0_dp   ! smallest host r that reached depot [m]
  real(dp), save              :: rmax_fail   = 0.0_dp   ! largest host r that failed [m] (threshold sits between)
  real(dp), save              :: sum_rho_rain = 0.0_dp  ! Sum of combined density rho_eff over rain-out
                                                        ! colonies -> measured mean = sum/n_rainout

  ! ---- Host registry (the colony droplets; only cell-bearing droplets exist) ----
  ! Hosts are the transported ACTIVE entity; member cells sync to h_z and read
  ! h_fate, so destruction needs no per-host membership lists.
  integer, parameter :: HF_ALIVE = 0, HF_DESICCATE = 1, HF_DEAD_BOT = 2, HF_DEAD_TOP = 3
  integer, parameter :: HF_RAINOUT = 4   ! colony heavy enough to sediment to depot -> cells DEPOT
  integer,  save              :: max_hosts = 0, n_hosts = 0
  real(dp), allocatable, save :: h_rhost(:)        ! droplet radius [m]
  integer,  allocatable, save :: h_capac(:)        ! packing capacity C
  integer,  allocatable, save :: h_occ(:)          ! occupancy (0 = free slot)
  real(dp), allocatable, save :: h_z(:)            ! altitude [m]
  integer,  allocatable, save :: h_lev(:)          ! current grid level
  integer,  allocatable, save :: h_fate(:)         ! HF_* this step
  real(dp), allocatable, save :: h_uconv(:)        ! per-colony OU convective state [-]

  ! Per-level total droplet number density [m^-3] = sum over the spectrum.
  ! "Liquid present" (a host can exist) where this exceeds N_LIQ_MIN.
  real(dp), allocatable, save :: n_total(:)
  real(dp), parameter        :: N_LIQ_MIN = 1.0e3_dp   ! [m^-3]

  ! ---- Carrying-capacity ceiling (reference-volume, per level) ----
  ! avail_lev(k) = n_total(k) * A_ref * dz(k) = droplets available in the
  ! simulated air parcel at level k.  occ_lev(k) = droplets currently occupied.
  ! Population scales with A_ref (a knob; absolute count is NOT a prediction).
  real(dp), allocatable, save :: dz_lev(:)         ! layer thickness [m]
  real(dp), allocatable, save :: avail_lev(:)      ! available droplets per level
  integer,  allocatable, save :: occ_lev(:)        ! occupied droplets per level
  real(dp), save, public :: A_ref = 1.0e-8_dp      ! reference area [m^2] (population scale)

  ! ---- Evolution / reproduction parameters (knobs) ----
  real(dp), save, public :: X_init_s   = 86400.0_dp   ! seed reproduction half-life [s] (1 day)
  real(dp), save, public :: mut_sigma_r = 0.10_dp     ! lognormal mutation width, r_cell
  real(dp), save, public :: mut_sigma_X = 0.10_dp     ! lognormal mutation width, X
  real(dp), save, public :: r_min      = 0.2e-6_dp    ! min viable cell radius [m] (Seager floor)
  real(dp), save, public :: r_max      = 2.0e-5_dp    ! cap to keep r_cell on the grid [m]
  real(dp), save, public :: X_min_s    = 3600.0_dp    ! min reproduction half-life [s]
  real(dp), save, public :: X_max_s    = 8.64e6_dp    ! max reproduction half-life [s] (100 days)
  ! ---- Mortality / turnover (closes the Seager cycle) ----
  ! cell_half_s: baseline ACTIVE-cell half-life (Yates-style); death frees a colony
  ! slot -> turnover keeps fission (and so evolution) running at capacity.
  ! host_half_s: host-droplet lifetime; on host death its cells DESICCATE back to
  ! DORMANT spores -> replenishes the seed bank (active -> spore -> germinate loop).
  real(dp), save, public :: cell_half_s = 8.64e5_dp   ! ACTIVE-cell half-life [s] (10 days)
  real(dp), save, public :: host_half_s = 4.32e5_dp   ! host-droplet lifetime [s] (5 days)
  integer,  save              :: this_step = 0

  ! ---- Trait-distribution diagnostic (optional) ----
  ! When on, evolve_test writes per-snapshot histograms of the two heritable
  ! traits (cell radius r, reproduction half-life X) over the ACTIVE population,
  ! plus a summary CSV, so the distribution (not just the mean) can be tracked.
  logical, save, public :: trait_dump = .false.
  integer, parameter    :: NTBIN = 24                 ! histogram bins per trait
  ! Amortized free-lists for O(1) slot reuse (rebuilt by one scan when drained).
  integer,  allocatable, save :: org_free(:),  host_free(:)
  integer,  save              :: org_free_n  = 0, org_free_ptr  = 1
  integer,  save              :: host_free_n = 0, host_free_ptr = 1

  ! Altitude bands [m] (knobs; defaults from plan section 6 / venus_phase2_plan).
  real(dp), save, public :: z_lethal_floor = 33.0e3_dp   ! below haze -> DEAD
  real(dp), save, public :: z_depot_lo     = 33.0e3_dp   ! depot/haze band bottom
  real(dp), save, public :: z_depot_hi     = 48.0e3_dp   ! depot/haze band top = cloud base
  real(dp), save, public :: z_domain_top   = 85.0e3_dp   ! above cloud top -> lost (DEAD)

  public :: set_rho_spore, set_rng_seed
  public :: bio_init, seed_spores, transport_step, transport_test, bio_cleanup
  public :: lifecycle_step, lifecycle_test, settling_check
  public :: step3_init, evolve_step, evolve_test

contains

  subroutine set_rho_spore(val)
    real(dp), intent(in) :: val
    if (val > 0.0_dp) rho_spore = val
  end subroutine set_rho_spore

  ! Seed the intrinsic RNG deterministically from a single integer, so ensemble
  ! members (different iseed) are independent yet each run is reproducible.
  subroutine set_rng_seed(iseed)
    integer, intent(in) :: iseed
    integer :: n, i
    integer, allocatable :: s(:)
    call random_seed(size=n)
    allocate(s(n))
    do i = 1, n
      s(i) = iseed * 1009 + i * 2718 + 12345
    end do
    call random_seed(put=s)
  end subroutine set_rng_seed

  ! Box-Muller standard normal variate.
  subroutine normal_rand(g)
    real(dp), intent(out) :: g
    real(dp) :: u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, 1.0e-15_dp)
    g  = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * PI * u2)
  end subroutine normal_rand

  ! Build the per-level total droplet density from the (already reconstructed)
  ! host spectrum.  Call after cloud_build_spectrum, before running a lifecycle.
  subroutine bio_init()
    integer :: k
    if (allocated(n_total)) deallocate(n_total)
    allocate(n_total(nz))
    do k = 1, nz
      n_total(k) = sum(spec(:,k))
    end do
  end subroutine bio_init

  ! Is there liquid (any host droplets) at altitude zq?
  logical function liquid_here(zq) result(yes)
    real(dp), intent(in) :: zq
    yes = .false.
    if (.not. allocated(n_total)) return
    yes = col_interp(n_total, zq) > N_LIQ_MIN
  end function liquid_here

  ! (Re)allocate the organism store to n slots, all DEAD.
  ! Allocate the organism store to nmax slots (the CAPACITY), all DEAD.
  subroutine org_alloc(nmax)
    integer, intent(in) :: nmax
    call free_orgs()
    max_orgs = nmax
    allocate(o_state(nmax), o_z(nmax), o_rcell(nmax), o_age(nmax), o_rhost(nmax), &
             o_capac(nmax), o_X(nmax), o_host(nmax), o_birthstep(nmax), o_uconv(nmax), org_free(nmax))
    o_state = ST_DEAD; o_z = 0.0_dp; o_rcell = 0.0_dp; o_age = 0.0_dp
    o_rhost = 0.0_dp; o_capac = 0; o_X = X_init_s; o_host = 0; o_birthstep = 0; o_uconv = 0.0_dp
    n_orgs = 0; n_lost_top = 0; n_lost_bot = 0; n_germ = 0; n_dead_dormant = 0
    n_birth = 0; n_eject = 0; n_dead_cell = 0; n_host_die = 0
    n_rainout = 0; n_desicc_place = 0; n_cells_depot = 0
    rmin_settle = 0.0_dp; rmax_fail = 0.0_dp; sum_rho_rain = 0.0_dp
    org_free_n = 0; org_free_ptr = 1
  end subroutine org_alloc

  subroutine free_orgs()
    if (allocated(o_state))     deallocate(o_state)
    if (allocated(o_z))         deallocate(o_z)
    if (allocated(o_rcell))     deallocate(o_rcell)
    if (allocated(o_age))       deallocate(o_age)
    if (allocated(o_rhost))     deallocate(o_rhost)
    if (allocated(o_capac))     deallocate(o_capac)
    if (allocated(o_X))         deallocate(o_X)
    if (allocated(o_host))      deallocate(o_host)
    if (allocated(o_birthstep)) deallocate(o_birthstep)
    if (allocated(o_uconv))     deallocate(o_uconv)
    if (allocated(org_free))    deallocate(org_free)
    max_orgs = 0; n_orgs = 0
  end subroutine free_orgs

  subroutine bio_cleanup()
    call free_orgs()
    if (allocated(h_rhost))   deallocate(h_rhost)
    if (allocated(h_capac))   deallocate(h_capac)
    if (allocated(h_occ))     deallocate(h_occ)
    if (allocated(h_z))       deallocate(h_z)
    if (allocated(h_lev))     deallocate(h_lev)
    if (allocated(h_fate))    deallocate(h_fate)
    if (allocated(h_uconv))   deallocate(h_uconv)
    if (allocated(host_free)) deallocate(host_free)
    if (allocated(n_total))   deallocate(n_total)
    if (allocated(dz_lev))    deallocate(dz_lev)
    if (allocated(avail_lev)) deallocate(avail_lev)
    if (allocated(occ_lev))   deallocate(occ_lev)
    max_hosts = 0; n_hosts = 0
  end subroutine bio_cleanup

  ! Seed n DORMANT spores of radius rcell, uniformly in altitude band [zlo,zhi].
  ! Allocates the store with headroom (nmax >= n) so daughters have free slots.
  subroutine seed_spores(n, rcell, zlo, zhi, nmax)
    integer,  intent(in)           :: n
    real(dp), intent(in)           :: rcell, zlo, zhi
    integer,  intent(in), optional :: nmax
    integer  :: i, cap
    real(dp) :: u
    cap = n
    if (present(nmax)) cap = max(nmax, n)
    call org_alloc(cap)
    n_orgs = n
    do i = 1, n
      o_state(i) = ST_DORMANT
      call random_number(u)
      o_z(i)     = zlo + u * (zhi - zlo)
      o_rcell(i) = rcell
      o_age(i)   = 0.0_dp
      o_X(i)     = X_init_s
      if (do_convection) call normal_rand(o_uconv(i))   ! stationary N(0,1) OU state (no RNG draw if off)
    end do
  end subroutine seed_spores

  ! Advance every live organism one step dt by free-spore transport:
  !   dz = (w - v_settle + dKzz/dz)*dt + sqrt(2*Kzz*dt)*N(0,1)
  ! The dKzz/dz drift keeps the random walk well-mixed in non-uniform Kzz.
  ! Absorbing boundaries: below z_lethal_floor or above z_domain_top -> DEAD.
  subroutine transport_step(dt)
    real(dp), intent(in) :: dt
    integer  :: i
    real(dp) :: zq, wq, kq, vset, drift, noise, g
    do i = 1, n_orgs
      if (o_state(i) == ST_DEAD) cycle
      zq    = o_z(i)
      wq    = col_interp(w, zq)
      kq    = max(col_interp(Kzz, zq), 0.0_dp)
      vset  = settling_velocity_z(zq, o_rcell(i), rho_spore)
      drift = (wq - vset + dKzz_dz(zq)) * dt
      call normal_rand(g)
      noise = sqrt(2.0_dp * kq * dt) * g
      zq    = zq + drift + noise
      if (zq <= z_lethal_floor) then
        o_state(i) = ST_DEAD; n_lost_bot = n_lost_bot + 1
      else if (zq >= z_domain_top) then
        o_state(i) = ST_DEAD; n_lost_top = n_lost_top + 1
      else
        o_z(i)   = zq
        o_age(i) = o_age(i) + dt
      end if
    end do
  end subroutine transport_step

  ! Tally live organisms by altitude band + altitude stats.
  subroutine org_stats(nalive, ndepot, ncloud, zmean, zmin, zmax)
    integer,  intent(out) :: nalive, ndepot, ncloud
    real(dp), intent(out) :: zmean, zmin, zmax
    integer  :: i
    real(dp) :: zsum
    nalive = 0; ndepot = 0; ncloud = 0; zsum = 0.0_dp
    zmin = huge(1.0_dp); zmax = -huge(1.0_dp)
    do i = 1, n_orgs
      if (o_state(i) == ST_DEAD) cycle
      nalive = nalive + 1
      zsum   = zsum + o_z(i)
      zmin   = min(zmin, o_z(i)); zmax = max(zmax, o_z(i))
      if (o_z(i) < z_depot_hi) then; ndepot = ndepot + 1; else; ncloud = ncloud + 1; end if
    end do
    if (nalive > 0) then; zmean = zsum / real(nalive, dp); else; zmean = 0.0_dp; end if
  end subroutine org_stats

  ! Verification driver: settling vs radius at the peak-spectrum level, for a
  ! host droplet (rho_host) and a desiccated spore (rho_spore).  Reports Kn and
  ! the Stokes-only vs slip-corrected speeds so the slip regime is explicit.
  subroutine settling_check(unit)
    integer, intent(in) :: unit
    integer,  parameter :: NR = 7
    real(dp), parameter :: r_um(NR) = &
         [0.1_dp, 0.2_dp, 0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp, 10.0_dp]
    real(dp) :: Npk, lam, mu, kn, Cc, vst, vsl, rr, rho_h
    integer  :: i, k, kpk

    Npk = 0.0_dp; kpk = 1
    do k = 1, nz
      if (sum(spec(:,k)) > Npk) then; Npk = sum(spec(:,k)); kpk = k; end if
    end do
    lam   = mean_free_path(T(kpk), P(kpk))
    mu    = viscosity_co2(T(kpk))
    rho_h = rho_host(kpk)
    if (rho_h <= 0.0_dp) rho_h = 1900.0_dp        ! H2SO4 droplet fallback

    write(unit,'(a)') ' '
    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a)') ' SETTLING  v = (2/9)(rho_p-rho_air) g r^2/mu * Cc(r)  [Stokes+slip]'
    write(unit,'(a,f6.2,a,f7.2,a,es10.3,a)') '  z = ', z(kpk)/1000.0_dp, &
         ' km,  T = ', T(kpk), ' K,  P = ', P(kpk), ' Pa'
    write(unit,'(a,es10.3,a,es10.3,a,f7.4,a)') '  rho_air = ', rho_air(kpk), &
         ' kg/m3,  mu = ', mu, ' Pa.s,  lambda = ', lam*1.0e6_dp, ' um'
    write(unit,'(a,f6.0,a,f6.0,a)') '  rho_host = ', rho_h, &
         ' kg/m3,  rho_spore = ', rho_spore, ' kg/m3'
    write(unit,'(a)') '   r[um]    Kn      Cc     host: vStokes  vSlip [mm/s]   spore: vSlip [mm/s]'
    do i = 1, NR
      rr  = r_um(i) * 1.0e-6_dp
      kn  = lam / rr
      Cc  = cunningham(rr, lam)
      vst = (2.0_dp/9.0_dp) * (rho_h - rho_air(kpk)) * GRAV_VENUS * rr*rr / mu
      vsl = settling_velocity(kpk, rr, rho_h)
      write(unit,'(2x,f7.2,2x,f7.3,2x,f6.2,4x,es11.3,2x,es11.3,4x,es11.3)') &
           r_um(i), kn, Cc, vst*1.0e3_dp, vsl*1.0e3_dp, &
           settling_velocity(kpk, rr, rho_spore)*1.0e3_dp
    end do
    write(unit,'(a)') '  (Cc->1 as r>>lambda; slip raises fall speed most for small r)'
    write(unit,'(a)') '----------------------------------------------------------'
  end subroutine settling_check

  ! Step-1 verification: seed depot spores, advect, watch the altitude
  ! distribution evolve (loft vs settle vs boundary loss).  No biology.
  subroutine transport_test(unit, n, rcell, dt, nsteps, nout)
    integer,  intent(in) :: unit, n, nsteps, nout
    real(dp), intent(in) :: rcell, dt
    integer  :: step, nalive, ndepot, ncloud
    real(dp) :: zmean, zmin, zmax, zseed, wseed, vseed

    call seed_spores(n, rcell, z_depot_lo, z_depot_hi)
    zseed = 0.5_dp * (z_depot_lo + z_depot_hi)
    wseed = col_interp(w, zseed)
    vseed = settling_velocity_z(zseed, rcell, rho_spore)

    write(unit,'(a)') ' '
    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a)') ' TRANSPORT TEST (step 1: free-spore kinematics, no biology)'
    write(unit,'(a,i0,a,f6.3,a)') '  ', n, ' DORMANT spores, r_cell = ', &
         rcell*1.0e6_dp, ' um'
    write(unit,'(a,f5.1,a,f5.1,a)') '  seeded uniformly in depot band ', &
         z_depot_lo/1000.0_dp, '-', z_depot_hi/1000.0_dp, ' km'
    write(unit,'(a,es10.3,a,es10.3,a)') '  at band mid: w = ', wseed, &
         ' m/s,  v_settle = ', vseed, ' m/s (net up if w>v_settle)'
    write(unit,'(a,es9.2,a,i0,a,es9.2,a)') '  dt = ', dt, ' s x ', nsteps, &
         ' steps = ', dt*real(nsteps,dp)/86400.0_dp, ' Earth-days'
    write(unit,'(a)') '   t[day]  nalive  ndepot  ncloud   <z>[km]  zmin   zmax   lostBot lostTop'
    do step = 0, nsteps
      if (mod(step, nout) == 0 .or. step == nsteps) then
        call org_stats(nalive, ndepot, ncloud, zmean, zmin, zmax)
        write(unit,'(2x,f7.2,2x,i7,2x,i6,2x,i6,2x,f7.2,2x,f6.2,2x,f6.2,2x,i7,1x,i7)') &
             dt*real(step,dp)/86400.0_dp, nalive, ndepot, ncloud, &
             zmean/1000.0_dp, zmin/1000.0_dp, zmax/1000.0_dp, n_lost_bot, n_lost_top
      end if
      if (step < nsteps) call transport_step(dt)
    end do
    write(unit,'(a)') '----------------------------------------------------------'
  end subroutine transport_test

  !====================================================================
  ! Step-2 lifecycle: 3-state machine with germination + transitions
  !   (no division yet — that is step 3, so occupancy stays 1 and the host
  !    reduces to the per-cell radius o_rhost).  Per-step algorithm follows
  !    planning/venus_phase2_plan.md section 5.
  !====================================================================

  ! Advance every live organism one step dt; Y_s = dormancy survival time [s].
  subroutine lifecycle_step(dt, Y_s)
    real(dp), intent(in) :: dt, Y_s
    integer  :: i, k, capac
    real(dp) :: zq, wq, kq, vset, drift, noise, g, rho_h, r_host
    logical  :: found, has_liq, in_haze

    do i = 1, n_orgs
      if (o_state(i) == ST_DEAD) cycle

      ! --- 1. transport by state ---
      zq = o_z(i)
      wq = col_interp(w, zq)
      if (o_state(i) == ST_ACTIVE) then
        ! ride the host droplet (H2SO4 density); advection + host settling
        rho_h = col_interp(rho_host, zq)
        if (rho_h <= 0.0_dp) rho_h = 1900.0_dp
        vset = settling_velocity_z(zq, o_rhost(i), rho_h)
        zq   = zq + (wq - vset) * dt
      else                                   ! DORMANT or DEPOT: own settling + eddy walk
        kq    = max(col_interp(Kzz, zq), 0.0_dp)
        vset  = settling_velocity_z(zq, o_rcell(i), rho_spore)
        drift = (wq - vset + dKzz_dz(zq)) * dt
        call normal_rand(g)
        noise = sqrt(2.0_dp * kq * dt) * g
        zq    = zq + drift + noise
      end if

      ! --- 2. domain boundaries (lethal) ---
      if (zq <= z_lethal_floor) then
        o_state(i) = ST_DEAD; n_lost_bot = n_lost_bot + 1; cycle
      else if (zq >= z_domain_top) then
        o_state(i) = ST_DEAD; n_lost_top = n_lost_top + 1; cycle
      end if
      o_z(i) = zq

      in_haze = zq < z_depot_hi
      has_liq = liquid_here(zq)

      ! --- 3. state transitions / germination ---
      select case (o_state(i))
      case (ST_ACTIVE)
        ! host destroyed when carried to a level with no liquid -> desiccate
        if (.not. has_liq) then
          o_state(i) = ST_DORMANT; o_rhost(i) = 0.0_dp; o_capac(i) = 0; o_age(i) = 0.0_dp
        end if

      case (ST_DORMANT)
        if (in_haze) then
          o_state(i) = ST_DEPOT; o_age(i) = 0.0_dp          ! settle into depot; clock off
        else
          o_age(i) = o_age(i) + dt                          ! Y-day clock runs in transit zone
          if (o_age(i) > Y_s) then
            o_state(i) = ST_DEAD; n_dead_dormant = n_dead_dormant + 1
          else if (has_liq) then                            ! 4. germination attempt
            k = col_level(zq)
            call draw_host(k, o_rcell(i), r_host, capac, found)
            if (found) then
              o_state(i) = ST_ACTIVE; o_rhost(i) = r_host
              o_capac(i) = capac;     o_age(i)   = 0.0_dp
              n_germ = n_germ + 1
            end if
          end if
        end if

      case (ST_DEPOT)
        if (.not. in_haze) then
          o_state(i) = ST_DORMANT; o_age(i) = 0.0_dp        ! up-leak into transit; clock on
        end if                                              ! else stay DEPOT (clock off)
      end select
    end do
  end subroutine lifecycle_step

  ! Per-state population counts + mean altitudes.
  subroutine lifecycle_stats(nact, ndor, ndep, zact, zall)
    integer,  intent(out) :: nact, ndor, ndep
    real(dp), intent(out) :: zact, zall
    integer  :: i, nall
    real(dp) :: sact, sall
    nact = 0; ndor = 0; ndep = 0; nall = 0; sact = 0.0_dp; sall = 0.0_dp
    do i = 1, n_orgs
      select case (o_state(i))
      case (ST_ACTIVE);  nact = nact + 1; sact = sact + o_z(i)
      case (ST_DORMANT); ndor = ndor + 1
      case (ST_DEPOT);   ndep = ndep + 1
      case default;      cycle                              ! DEAD
      end select
      nall = nall + 1; sall = sall + o_z(i)
    end do
    if (nact > 0) then; zact = sact / real(nact, dp); else; zact = 0.0_dp; end if
    if (nall > 0) then; zall = sall / real(nall, dp); else; zall = 0.0_dp; end if
  end subroutine lifecycle_stats

  ! Step-2 verification: seed spores, run the full state machine, report the
  ! per-state populations and cumulative germinations / deaths over time.
  ! (No reproduction yet, so total alive only decays — step 3 adds fission.)
  subroutine lifecycle_test(unit, n, rcell, dt, nsteps, nout, Ydays, zlo, zhi)
    integer,  intent(in) :: unit, n, nsteps, nout
    real(dp), intent(in) :: rcell, dt, Ydays, zlo, zhi
    integer  :: step, nact, ndor, ndep
    real(dp) :: zact, zall, Y_s

    Y_s = Ydays * 86400.0_dp
    call seed_spores(n, rcell, zlo, zhi)

    write(unit,'(a)') ' '
    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a)') ' LIFECYCLE TEST (step 2: 3-state machine + germination, no fission)'
    write(unit,'(a,i0,a,f6.3,a)') '  ', n, ' spores, r_cell = ', rcell*1.0e6_dp, ' um'
    write(unit,'(a,f5.1,a,f5.1,a,f6.2,a)') '  seeded in ', zlo/1000.0_dp, '-', &
         zhi/1000.0_dp, ' km;  Y(dormancy) = ', Ydays, ' Earth-days'
    write(unit,'(a,es9.2,a,i0,a,f7.1,a)') '  dt = ', dt, ' s x ', nsteps, &
         ' = ', dt*real(nsteps,dp)/86400.0_dp, ' Earth-days'
    write(unit,'(a)') '   t[day]  nACT  nDOR  nDEP  alive  <z>act <z>all  germ   dDorm lostB lostT'
    do step = 0, nsteps
      if (mod(step, nout) == 0 .or. step == nsteps) then
        call lifecycle_stats(nact, ndor, ndep, zact, zall)
        write(unit,'(2x,f7.2,3(1x,i5),1x,i6,2(1x,f6.2),1x,i7,3(1x,i5))') &
             dt*real(step,dp)/86400.0_dp, nact, ndor, ndep, nact+ndor+ndep, &
             zact/1000.0_dp, zall/1000.0_dp, n_germ, n_dead_dormant, n_lost_bot, n_lost_top
      end if
      if (step < nsteps) call lifecycle_step(dt, Y_s)
    end do
    write(unit,'(a)') '----------------------------------------------------------'
  end subroutine lifecycle_test

  !====================================================================
  ! Step-3 evolution: fission + mutation + per-host capacity/ejection +
  !   reference-volume carrying-capacity ceiling + slot reuse.
  !   Per-step algorithm follows planning/venus_phase2_plan.md sections 2.3, 5.
  !====================================================================

  ! Build per-level layer thickness, available-droplet ceiling, and the host
  ! registry.  Call AFTER seed_spores (needs max_orgs) and bio_init (n_total).
  subroutine step3_init()
    integer :: k
    if (.not. allocated(n_total)) call bio_init()
    if (allocated(dz_lev))    deallocate(dz_lev)
    if (allocated(avail_lev)) deallocate(avail_lev)
    if (allocated(occ_lev))   deallocate(occ_lev)
    allocate(dz_lev(nz), avail_lev(nz), occ_lev(nz))
    do k = 1, nz
      if (k == 1) then
        dz_lev(k) = z(2) - z(1)
      else if (k == nz) then
        dz_lev(k) = z(nz) - z(nz-1)
      else
        dz_lev(k) = 0.5_dp * (z(k+1) - z(k-1))
      end if
      avail_lev(k) = n_total(k) * A_ref * dz_lev(k)
      occ_lev(k)   = 0
    end do
    if (allocated(h_rhost))   deallocate(h_rhost, h_capac, h_occ, h_z, h_lev, h_fate, h_uconv)
    if (allocated(host_free)) deallocate(host_free)
    max_hosts = max_orgs
    allocate(h_rhost(max_hosts), h_capac(max_hosts), h_occ(max_hosts), &
             h_z(max_hosts), h_lev(max_hosts), h_fate(max_hosts), h_uconv(max_hosts), host_free(max_hosts))
    h_rhost = 0.0_dp; h_capac = 0; h_occ = 0; h_z = 0.0_dp; h_lev = 0; h_fate = HF_ALIVE; h_uconv = 0.0_dp
    n_hosts = 0; this_step = 0
    org_free_n = 0; org_free_ptr = 1; host_free_n = 0; host_free_ptr = 1
  end subroutine step3_init

  ! Return a free organism slot: extend the high-water mark, else serve a DEAD
  ! slot from the free-list (rebuilt by one O(N) scan when drained -> amortized
  ! O(1)).  0 if the store is genuinely full of live organisms.
  integer function new_org_slot() result(j)
    integer :: i
    if (n_orgs < max_orgs) then
      n_orgs = n_orgs + 1; j = n_orgs; return
    end if
    if (org_free_ptr > org_free_n) then              ! drained -> rebuild
      org_free_n = 0
      do i = 1, max_orgs
        if (o_state(i) == ST_DEAD) then; org_free_n = org_free_n + 1; org_free(org_free_n) = i; end if
      end do
      org_free_ptr = 1
    end if
    if (org_free_ptr > org_free_n) then; j = 0; return; end if   ! truly full
    j = org_free(org_free_ptr); org_free_ptr = org_free_ptr + 1
  end function new_org_slot

  ! Return a free host slot (free = occupancy 0); same amortized free-list.
  integer function new_host_slot() result(j)
    integer :: i
    if (n_hosts < max_hosts) then
      n_hosts = n_hosts + 1; j = n_hosts; return
    end if
    if (host_free_ptr > host_free_n) then
      host_free_n = 0
      do i = 1, max_hosts
        if (h_occ(i) == 0) then; host_free_n = host_free_n + 1; host_free(host_free_n) = i; end if
      end do
      host_free_ptr = 1
    end if
    if (host_free_ptr > host_free_n) then; j = 0; return; end if
    j = host_free(host_free_ptr); host_free_ptr = host_free_ptr + 1
  end function new_host_slot

  ! Lognormal multiplicative mutation, clamped to [vmin,vmax].
  real(dp) function mutate_log(val, sigma, vmin, vmax) result(out)
    real(dp), intent(in) :: val, sigma, vmin, vmax
    real(dp) :: g
    call normal_rand(g)
    out = val * exp(sigma * g)
    out = max(vmin, min(vmax, out))
  end function mutate_log

  ! Combined (cell-laden) droplet density [kg/m^3] for host h: volume-weighted mix
  ! of hydrated cells and H2SO4 solution.  A full colony's cell volume fraction is
  ! ~phi, so f_cells = phi * occ/capac (cells are LIGHTER than acid, so a fuller
  ! colony is slightly less dense -> settles marginally slower).
  real(dp) function host_density_eff(h, zq) result(rho_eff)
    integer,  intent(in) :: h
    real(dp), intent(in) :: zq
    real(dp) :: rho_liq, f_cells
    rho_liq = col_interp(rho_host, zq); if (rho_liq <= 0.0_dp) rho_liq = 1900.0_dp
    f_cells = 0.0_dp
    if (h_capac(h) > 0) f_cells = phi_pack * real(h_occ(h), dp) / real(h_capac(h), dp)
    rho_eff = f_cells * rho_cell_wet + (1.0_dp - f_cells) * rho_liq
  end function host_density_eff

  ! Depot return-leg feasibility: would a droplet of radius r, density rho, released
  ! at z_from, sediment all the way to the depot top (z_depot_hi) against the
  ! vertical wind?  True iff its settling speed beats the upward wind at EVERY level
  ! on the path down (net-downward everywhere -> not trapped by the cloud-base
  ! updraft).  This is the physical test of the Seager cycle's return leg.
  logical function settles_to_depot(z_from, r, rho) result(yes)
    real(dp), intent(in) :: z_from, r, rho
    integer :: k, klo, khi
    yes = .true.
    if (z_from <= z_depot_hi) return          ! already at/below depot top
    klo = col_level(z_depot_hi)
    khi = col_level(z_from)
    do k = klo, khi
      if (settling_velocity(k, r, rho) <= w(k)) then; yes = .false.; return; end if
    end do
  end function settles_to_depot

  ! Classify a dying colony (host h) as rain-out (heavy enough to reach the depot)
  ! or in-place desiccation, and record the verdict + host-radius envelope.
  integer function desiccation_fate(h) result(fate)
    integer, intent(in) :: h
    real(dp) :: rho
    rho = host_density_eff(h, h_z(h))
    if (settles_to_depot(h_z(h), h_rhost(h), rho)) then
      fate = HF_RAINOUT;   n_rainout = n_rainout + 1
      sum_rho_rain = sum_rho_rain + rho          ! measured combined density of this rained-out colony
      if (rmin_settle == 0.0_dp) then; rmin_settle = h_rhost(h); else; rmin_settle = min(rmin_settle, h_rhost(h)); end if
    else
      fate = HF_DESICCATE; n_desicc_place = n_desicc_place + 1
      rmax_fail = max(rmax_fail, h_rhost(h))
    end if
  end function desiccation_fate

  ! One evolution step.  Passes: A transport free spores; B transport hosts, set
  ! fate (incl. finite-lifetime droplet death -> desiccate); C apply host fate to
  ! ACTIVE cells; C2 baseline cell mortality; D germination (ceiling-limited);
  ! E fission + mutation (a daughter ejected from a full colony dies).
  subroutine evolve_step(dt, Y_s)
    real(dp), intent(in) :: dt, Y_s
    integer  :: i, k, h, j, nnow, capac, s, nsub
    real(dp) :: zq, wq, kq, vset, drift, noise, g, rho_h, r_host, p, u, dt_sub, wconv
    logical  :: found, in_haze, dead_bot, dead_top, dried

    this_step = this_step + 1
    nsub   = merge(conv_nsub, 1, do_convection)     ! sub-step transport only if convecting
    dt_sub = dt / real(nsub, dp)

    ! ---- A: DORMANT/DEPOT transport (sub-stepped w/ convection) + transitions ----
    do i = 1, n_orgs
      if (o_state(i) /= ST_DORMANT .and. o_state(i) /= ST_DEPOT) cycle
      zq = o_z(i); dead_bot = .false.; dead_top = .false.
      do s = 1, nsub
        wq   = col_interp(w, zq)
        kq   = max(col_interp(Kzz, zq), 0.0_dp)
        vset = settling_velocity_z(zq, o_rcell(i), rho_spore)
        wconv = 0.0_dp
        if (do_convection) then
          call normal_rand(g); o_uconv(i) = conv_ar1(o_uconv(i), dt_sub, g)
          wconv = conv_sigma * conv_envelope(zq) * o_uconv(i)
        end if
        drift = (wq - vset + dKzz_dz(zq) + wconv) * dt_sub
        call normal_rand(g); noise = sqrt(2.0_dp * kq * dt_sub) * g
        zq = zq + drift + noise
        if (zq <= z_lethal_floor) then; dead_bot = .true.; exit
        else if (zq >= z_domain_top) then; dead_top = .true.; exit; end if
      end do
      if (dead_bot) then; o_state(i) = ST_DEAD; n_lost_bot = n_lost_bot + 1; cycle; end if
      if (dead_top) then; o_state(i) = ST_DEAD; n_lost_top = n_lost_top + 1; cycle; end if
      o_z(i)  = zq
      in_haze = zq < z_depot_hi
      if (o_state(i) == ST_DORMANT) then
        if (in_haze) then
          o_state(i) = ST_DEPOT; o_age(i) = 0.0_dp
        else
          o_age(i) = o_age(i) + dt
          if (o_age(i) > Y_s) then; o_state(i) = ST_DEAD; n_dead_dormant = n_dead_dormant + 1; end if
        end if
      else                                    ! DEPOT
        if (.not. in_haze) then; o_state(i) = ST_DORMANT; o_age(i) = 0.0_dp; end if
      end if
    end do

    ! ---- B: host transport (sub-stepped w/ convection), fate, occ_lev bookkeeping ----
    do h = 1, n_hosts
      if (h_occ(h) == 0) cycle
      zq = h_z(h); dead_bot = .false.; dead_top = .false.; dried = .false.
      do s = 1, nsub
        wq = col_interp(w, zq)
        if (do_convection) then
          rho_h = host_density_eff(h, zq)             ! cell-laden droplet density
        else
          rho_h = col_interp(rho_host, zq); if (rho_h <= 0.0_dp) rho_h = 1900.0_dp
        end if
        vset = settling_velocity_z(zq, h_rhost(h), rho_h)
        wconv = 0.0_dp
        if (do_convection) then
          call normal_rand(g); h_uconv(h) = conv_ar1(h_uconv(h), dt_sub, g)
          wconv = conv_sigma * conv_envelope(zq) * h_uconv(h)
        end if
        zq = zq + (wq - vset + wconv) * dt_sub
        if (zq <= z_lethal_floor) then; dead_bot = .true.; exit
        else if (zq >= z_domain_top) then; dead_top = .true.; exit
        else if (.not. liquid_here(zq)) then; dried = .true.; exit; end if
      end do
      if (dead_bot) then
        h_fate(h) = HF_DEAD_BOT; occ_lev(h_lev(h)) = occ_lev(h_lev(h)) - 1
      else if (dead_top) then
        h_fate(h) = HF_DEAD_TOP; occ_lev(h_lev(h)) = occ_lev(h_lev(h)) - 1
      else if (dried) then
        h_z(h) = zq                             ! desiccates where it arrived (convection may have
        ! carried it into/toward the depot); with convection the transport handles the return, so
        ! cells just desiccate in place (-> DORMANT); without it, use the settling-threshold proxy.
        if (do_convection) then; h_fate(h) = HF_DESICCATE
        else;                     h_fate(h) = desiccation_fate(h); end if
        occ_lev(h_lev(h)) = occ_lev(h_lev(h)) - 1
      else
        h_fate(h) = HF_ALIVE; h_z(h) = zq
        k = col_level(zq)
        if (k /= h_lev(h)) then
          occ_lev(h_lev(h)) = occ_lev(h_lev(h)) - 1
          occ_lev(k)        = occ_lev(k) + 1
          h_lev(h)          = k
        end if
        ! droplet finite lifetime (memoryless): evaporate.
        p = 1.0_dp - exp(-dt * log(2.0_dp) / host_half_s)
        call random_number(u)
        if (u < p) then
          if (do_convection) then; h_fate(h) = HF_DESICCATE
          else;                     h_fate(h) = desiccation_fate(h); end if
          occ_lev(h_lev(h)) = occ_lev(h_lev(h)) - 1
          n_host_die = n_host_die + 1
        end if
      end if
    end do

    ! ---- C: ACTIVE cells follow / release from their host ----
    do i = 1, n_orgs
      if (o_state(i) /= ST_ACTIVE) cycle
      h = o_host(i)
      select case (h_fate(h))
      case (HF_ALIVE)
        o_z(i) = h_z(h)
      case (HF_DESICCATE)
        o_state(i) = ST_DORMANT; o_z(i) = h_z(h)
        o_host(i) = 0; o_rhost(i) = 0.0_dp; o_capac(i) = 0; o_age(i) = 0.0_dp
        o_uconv(i) = h_uconv(h)                  ! released spore rides the colony's convective parcel
      case (HF_RAINOUT)                          ! sedimented to the haze depot
        o_state(i) = ST_DEPOT; o_z(i) = 0.5_dp * (z_depot_lo + z_depot_hi)
        o_host(i) = 0; o_rhost(i) = 0.0_dp; o_capac(i) = 0; o_age(i) = 0.0_dp
        n_cells_depot = n_cells_depot + 1
      case (HF_DEAD_BOT)
        o_state(i) = ST_DEAD; n_lost_bot = n_lost_bot + 1; o_host(i) = 0
      case (HF_DEAD_TOP)
        o_state(i) = ST_DEAD; n_lost_top = n_lost_top + 1; o_host(i) = 0
      end select
    end do
    do h = 1, n_hosts                          ! free destroyed hosts; reset fate
      if (h_occ(h) > 0 .and. h_fate(h) /= HF_ALIVE) h_occ(h) = 0
      h_fate(h) = HF_ALIVE
    end do

    ! ---- C2: baseline ACTIVE-cell mortality (memoryless half-life) ----
    ! Frees a colony slot; when a colony loses its last cell the droplet reverts to
    ! an empty cloud particle and its level slot returns to the available pool.
    do i = 1, n_orgs
      if (o_state(i) /= ST_ACTIVE) cycle
      if (o_birthstep(i) >= this_step) cycle   ! germinated/born this step
      p = 1.0_dp - exp(-dt * log(2.0_dp) / cell_half_s)
      call random_number(u)
      if (u >= p) cycle
      h = o_host(i)
      o_state(i) = ST_DEAD; n_dead_cell = n_dead_cell + 1
      if (h > 0) then
        h_occ(h) = h_occ(h) - 1
        if (h_occ(h) == 0) occ_lev(h_lev(h)) = occ_lev(h_lev(h)) - 1
      end if
      o_host(i) = 0; o_rhost(i) = 0.0_dp; o_capac(i) = 0
    end do

    ! ---- D: germination (DORMANT -> ACTIVE), ceiling-limited ----
    do i = 1, n_orgs
      if (o_state(i) /= ST_DORMANT) cycle
      zq = o_z(i)
      if (zq < z_depot_hi) cycle               ! haze: becomes depot next step
      if (.not. liquid_here(zq)) cycle
      k = col_level(zq)
      if (real(occ_lev(k)+1, dp) > avail_lev(k)) cycle   ! no whole host droplet available here
      call draw_host(k, o_rcell(i), r_host, capac, found)
      if (.not. found) cycle
      j = new_host_slot()
      if (j == 0) cycle
      h_rhost(j) = r_host; h_capac(j) = capac; h_occ(j) = 1
      h_z(j) = zq; h_lev(j) = k; h_fate(j) = HF_ALIVE
      if (do_convection) then; call normal_rand(g); h_uconv(j) = g; else; h_uconv(j) = 0.0_dp; end if
      occ_lev(k) = occ_lev(k) + 1
      o_state(i) = ST_ACTIVE; o_host(i) = j; o_rhost(i) = r_host; o_capac(i) = capac
      o_age(i) = 0.0_dp; o_birthstep(i) = this_step
      n_germ = n_germ + 1
    end do

    ! ---- E: fission + mutation; a daughter ejected from a full colony DIES ----
    ! A division's surplus daughter is a bare cell, not a spore: expelled from the
    ! host droplet it has no coat and dies (it is simply not instantiated). So a
    ! full colony stops growing and the population is bounded by host capacity,
    ! with no ejected-spore contribution to the dormant pool.
    nnow = n_orgs
    do i = 1, nnow
      if (o_state(i) /= ST_ACTIVE) cycle
      if (o_birthstep(i) >= this_step) cycle   ! germinated/born this step
      p = 1.0_dp - exp(-dt * log(2.0_dp) / o_X(i))
      call random_number(u)
      if (u >= p) cycle
      h = o_host(i)
      if (h_occ(h) >= h_capac(h)) then         ! colony full -> daughter ejected, dies
        n_eject = n_eject + 1
        cycle
      end if
      j = new_org_slot()
      if (j == 0) cycle                        ! store full -> drop birth
      o_rcell(j)     = mutate_log(o_rcell(i), mut_sigma_r, r_min, r_max)
      o_X(j)         = mutate_log(o_X(i),     mut_sigma_X, X_min_s, X_max_s)
      o_age(j)       = 0.0_dp
      o_birthstep(j) = this_step
      o_state(j) = ST_ACTIVE; o_host(j) = h; o_z(j) = h_z(h)
      o_rhost(j) = h_rhost(h); o_capac(j) = h_capac(h)
      h_occ(h)   = h_occ(h) + 1
      n_birth    = n_birth + 1
    end do
  end subroutine evolve_step

  ! Population by state + colony stats + mean ACTIVE traits.
  subroutine evolve_stats(nact, ndor, ndep, nhost, ncol, maxocc, rmean_um, Xmean_d, zact)
    integer,  intent(out) :: nact, ndor, ndep, nhost, ncol, maxocc
    real(dp), intent(out) :: rmean_um, Xmean_d, zact
    integer  :: i, h
    real(dp) :: sr, sx, sz
    nact = 0; ndor = 0; ndep = 0; sr = 0.0_dp; sx = 0.0_dp; sz = 0.0_dp
    do i = 1, n_orgs
      select case (o_state(i))
      case (ST_ACTIVE)
        nact = nact + 1
        sr = sr + o_rcell(i); sx = sx + o_X(i); sz = sz + o_z(i)
      case (ST_DORMANT); ndor = ndor + 1
      case (ST_DEPOT);   ndep = ndep + 1
      end select
    end do
    nhost = 0; ncol = 0; maxocc = 0
    do h = 1, n_hosts
      if (h_occ(h) > 0) then
        nhost = nhost + 1
        if (h_occ(h) > 1)      ncol   = ncol + 1
        if (h_occ(h) > maxocc) maxocc = h_occ(h)
      end if
    end do
    if (nact > 0) then
      rmean_um = sr / real(nact, dp) * 1.0e6_dp
      Xmean_d  = sx / real(nact, dp) / 86400.0_dp
      zact     = sz / real(nact, dp) / 1000.0_dp
    else
      rmean_um = 0.0_dp; Xmean_d = 0.0_dp; zact = 0.0_dp
    end if
  end subroutine evolve_stats

  ! Log-spaced bin edges (NTBIN bins) spanning [lo,hi]; edges(1..NTBIN+1).
  subroutine log_edges(lo, hi, edges)
    real(dp), intent(in)  :: lo, hi
    real(dp), intent(out) :: edges(NTBIN+1)
    integer  :: i
    real(dp) :: llo, lhi
    llo = log(lo); lhi = log(hi)
    do i = 0, NTBIN
      edges(i+1) = exp(llo + (lhi - llo) * real(i, dp) / real(NTBIN, dp))
    end do
  end subroutine log_edges

  ! Log-bin index for value v given the low/high edges (clamped to [1,NTBIN]).
  integer function log_bin(v, lo, hi) result(b)
    real(dp), intent(in) :: v, lo, hi
    real(dp) :: llo, lhi
    llo = log(lo); lhi = log(hi)
    b = 1 + int(real(NTBIN, dp) * (log(v) - llo) / (lhi - llo))
    if (b < 1)     b = 1
    if (b > NTBIN) b = NTBIN
  end function log_bin

  ! One snapshot: bin the ACTIVE population's traits into r/X histograms (units
  ! and mean/std to the summary), scaled by the mode multiplicity so a bin count
  ! reflects the true number of cells (a full colony carries h_occ identical
  ! members already, since each is its own org).  Writes one row per file.
  subroutine trait_snapshot(t_day, ur, ux, us, nact, ndor, ndep)
    real(dp), intent(in) :: t_day
    integer,  intent(in) :: ur, ux, us, nact, ndor, ndep
    integer  :: i, b, rc(NTBIN), xc(NTBIN)
    real(dp) :: sr, sr2, sx, sx2, rmean, rstd, xmean, xstd, rv, xv
    rc = 0; xc = 0; sr = 0.0_dp; sr2 = 0.0_dp; sx = 0.0_dp; sx2 = 0.0_dp
    do i = 1, n_orgs
      if (o_state(i) /= ST_ACTIVE) cycle
      rv = o_rcell(i) * 1.0e6_dp                    ! um
      xv = o_X(i) / 86400.0_dp                      ! days
      b = log_bin(o_rcell(i), r_min, r_max); rc(b) = rc(b) + 1
      b = log_bin(o_X(i),    X_min_s, X_max_s); xc(b) = xc(b) + 1
      sr = sr + rv; sr2 = sr2 + rv*rv
      sx = sx + xv; sx2 = sx2 + xv*xv
    end do
    if (nact > 0) then
      rmean = sr / real(nact, dp); rstd = sqrt(max(0.0_dp, sr2/real(nact,dp) - rmean*rmean))
      xmean = sx / real(nact, dp); xstd = sqrt(max(0.0_dp, sx2/real(nact,dp) - xmean*xmean))
    else
      rmean = 0.0_dp; rstd = 0.0_dp; xmean = 0.0_dp; xstd = 0.0_dp
    end if
    write(ur,'(f9.3,24(",",i0))') t_day, rc
    write(ux,'(f9.3,24(",",i0))') t_day, xc
    write(us,'(f9.3,3(",",i0),4(",",es12.5))') t_day, nact, ndor, ndep, &
         rmean, rstd, xmean, xstd
  end subroutine trait_snapshot

  ! Step-3 verification: seed spores, run fission+mutation+ceiling, report the
  ! per-state populations, colony stats and evolving mean traits over time.
  subroutine evolve_test(unit, nseed, rcell0, dt, nsteps, nout, Ydays, zlo, zhi, maxorgs)
    integer,  intent(in) :: unit, nseed, nsteps, nout, maxorgs
    real(dp), intent(in) :: rcell0, dt, Ydays, zlo, zhi
    integer  :: step, nact, ndor, ndep, nhost, ncol, maxocc
    real(dp) :: rmean_um, Xmean_d, zact, Y_s
    integer  :: ur, ux, us, b
    real(dp) :: redg(NTBIN+1), xedg(NTBIN+1)

    Y_s = Ydays * 86400.0_dp
    call seed_spores(nseed, rcell0, zlo, zhi, maxorgs)
    call step3_init()

    if (trait_dump) then                    ! open per-snapshot trait histograms
      call log_edges(r_min, r_max, redg)
      call log_edges(X_min_s, X_max_s, xedg)
      open(newunit=ur, file='traits_r.csv',       status='replace', action='write')
      open(newunit=ux, file='traits_X.csv',       status='replace', action='write')
      open(newunit=us, file='traits_summary.csv', status='replace', action='write')
      write(ur,'(a)',advance='no') 't_day'      ! header = bin CENTERS (geometric), um
      do b = 1, NTBIN; write(ur,'(",",es11.4)',advance='no') sqrt(redg(b)*redg(b+1))*1.0e6_dp; end do
      write(ur,'(a)') ''
      write(ux,'(a)',advance='no') 't_day'      ! header = bin CENTERS (geometric), days
      do b = 1, NTBIN; write(ux,'(",",es11.4)',advance='no') sqrt(xedg(b)*xedg(b+1))/86400.0_dp; end do
      write(ux,'(a)') ''
      write(us,'(a)') 't_day,nACT,nDOR,nDEP,r_mean_um,r_std_um,X_mean_d,X_std_d'
    end if

    write(unit,'(a)') ' '
    write(unit,'(a)') '----------------------------------------------------------'
    write(unit,'(a)') ' EVOLVE TEST (step 3: fission + mutation + carrying capacity + turnover)'
    write(unit,'(a,i0,a,f6.3,a,i0)') '  ', nseed, ' seed spores, r_cell0 = ', &
         rcell0*1.0e6_dp, ' um;  max_orgs = ', maxorgs
    write(unit,'(a,es9.2,a,es10.3,a)') '  A_ref = ', A_ref, &
         ' m2;  total avail droplets = ', sum(avail_lev), ' (population scale)'
    write(unit,'(a,f5.1,a,f5.1,a,f6.2,a,f6.3,a)') '  seed ', zlo/1000.0_dp, '-', &
         zhi/1000.0_dp, ' km;  Y = ', Ydays, ' d;  X0 = ', X_init_s/86400.0_dp, ' d'
    write(unit,'(a,es9.2,a,i0,a,f7.1,a)') '  dt = ', dt, ' s x ', nsteps, &
         ' = ', dt*real(nsteps,dp)/86400.0_dp, ' Earth-days'
    write(unit,'(a,f5.1,a,f5.1,a)') '  cell half-life = ', cell_half_s/86400.0_dp, &
         ' d;  host lifetime = ', host_half_s/86400.0_dp, ' d'
    write(unit,'(a)') '   t[day]  nACT  nDOR   nDEP  hosts colo maxC  <r>um  <X>d  <z>km   births  eject'
    do step = 0, nsteps
      if (mod(step, nout) == 0 .or. step == nsteps) then
        call evolve_stats(nact, ndor, ndep, nhost, ncol, maxocc, rmean_um, Xmean_d, zact)
        write(unit,'(2x,f7.2,1x,i6,1x,i6,1x,i6,1x,i5,1x,i5,1x,i4,1x,f6.3,1x,f5.2,1x,f6.2,1x,i9,1x,i7)') &
             dt*real(step,dp)/86400.0_dp, nact, ndor, ndep, nhost, ncol, maxocc, &
             rmean_um, Xmean_d, zact, n_birth, n_eject
        if (trait_dump) call trait_snapshot(dt*real(step,dp)/86400.0_dp, ur, ux, us, nact, ndor, ndep)
      end if
      if (step < nsteps) call evolve_step(dt, Y_s)
    end do
    write(unit,'(a,i0,a,i0,a,i0)') '  cumulative: germinations = ', n_germ, &
         ',  cell deaths = ', n_dead_cell, ',  host deaths = ', n_host_die
    write(unit,'(a)') '  --- depot return-leg (Seager cycle) feasibility ---'
    write(unit,'(a,i0,a,i0,a,f6.1,a)') '  colony deaths: rain-out to depot = ', n_rainout, &
         ',  desiccate in place = ', n_desicc_place, &
         '   (rain-out = ', 100.0_dp*real(n_rainout,dp)/real(max(n_rainout+n_desicc_place,1),dp), ' %)'
    write(unit,'(a,i0)') '  cells delivered to depot bank = ', n_cells_depot
    write(unit,'(a,f7.3,a,f7.3,a)') '  host-radius envelope: largest that FAILED = ', &
         rmax_fail*1.0e6_dp, ' um;  smallest that REACHED depot = ', rmin_settle*1.0e6_dp, ' um'
    write(unit,'(a,f8.1,a)') '  measured mean rho_eff of rained-out colonies = ', &
         sum_rho_rain / real(max(n_rainout, 1), dp), ' kg/m3'
    write(unit,'(a)') '----------------------------------------------------------'
    if (trait_dump) then; close(ur); close(ux); close(us); end if
  end subroutine evolve_test

end module bio_venus
