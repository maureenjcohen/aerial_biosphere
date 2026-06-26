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
  integer,  allocatable, save :: o_capac(:)        ! host packing capacity (for step-3 division)
  integer,  save              :: n_lost_top = 0, n_lost_bot = 0   ! boundary-loss tallies
  integer,  save              :: n_germ = 0, n_dead_dormant = 0   ! germinations / Y-clock deaths

  ! Per-level total droplet number density [m^-3] = sum over the spectrum.
  ! "Liquid present" (a host can exist) where this exceeds N_LIQ_MIN.
  real(dp), allocatable, save :: n_total(:)
  real(dp), parameter        :: N_LIQ_MIN = 1.0e3_dp   ! [m^-3]

  ! Altitude bands [m] (knobs; defaults from plan section 6 / venus_phase2_plan).
  real(dp), save, public :: z_lethal_floor = 33.0e3_dp   ! below haze -> DEAD
  real(dp), save, public :: z_depot_lo     = 33.0e3_dp   ! depot/haze band bottom
  real(dp), save, public :: z_depot_hi     = 48.0e3_dp   ! depot/haze band top = cloud base
  real(dp), save, public :: z_domain_top   = 85.0e3_dp   ! above cloud top -> lost (DEAD)

  public :: set_rho_spore
  public :: bio_init, seed_spores, transport_step, transport_test, bio_cleanup
  public :: lifecycle_step, lifecycle_test, settling_check

contains

  subroutine set_rho_spore(val)
    real(dp), intent(in) :: val
    if (val > 0.0_dp) rho_spore = val
  end subroutine set_rho_spore

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
  subroutine org_alloc(n)
    integer, intent(in) :: n
    call free_orgs()
    max_orgs = n
    allocate(o_state(n), o_z(n), o_rcell(n), o_age(n), o_rhost(n), o_capac(n))
    o_state = ST_DEAD; o_z = 0.0_dp; o_rcell = 0.0_dp; o_age = 0.0_dp
    o_rhost = 0.0_dp; o_capac = 0
    n_orgs = 0; n_lost_top = 0; n_lost_bot = 0; n_germ = 0; n_dead_dormant = 0
  end subroutine org_alloc

  subroutine free_orgs()
    if (allocated(o_state)) deallocate(o_state)
    if (allocated(o_z))     deallocate(o_z)
    if (allocated(o_rcell)) deallocate(o_rcell)
    if (allocated(o_age))   deallocate(o_age)
    if (allocated(o_rhost)) deallocate(o_rhost)
    if (allocated(o_capac)) deallocate(o_capac)
    max_orgs = 0; n_orgs = 0
  end subroutine free_orgs

  subroutine bio_cleanup()
    call free_orgs()
    if (allocated(n_total)) deallocate(n_total)
  end subroutine bio_cleanup

  ! Seed n DORMANT spores of radius rcell, uniformly in altitude band [zlo,zhi].
  subroutine seed_spores(n, rcell, zlo, zhi)
    integer,  intent(in) :: n
    real(dp), intent(in) :: rcell, zlo, zhi
    integer  :: i
    real(dp) :: u
    call org_alloc(n)
    n_orgs = n
    do i = 1, n
      o_state(i) = ST_DORMANT
      call random_number(u)
      o_z(i)     = zlo + u * (zhi - zlo)
      o_rcell(i) = rcell
      o_age(i)   = 0.0_dp
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

end module bio_venus
