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
  real(dp), allocatable, save :: o_rcell(:)        ! cell radius [m] (heritable; fixed in step 1)
  real(dp), allocatable, save :: o_age(:)          ! time in current dormant spell [s]
  integer,  save              :: n_lost_top = 0, n_lost_bot = 0   ! boundary-loss tallies

  ! Altitude bands [m] (knobs; defaults from plan section 6 / venus_phase2_plan).
  real(dp), save, public :: z_lethal_floor = 33.0e3_dp   ! below haze -> DEAD
  real(dp), save, public :: z_depot_lo     = 33.0e3_dp   ! depot/haze band bottom
  real(dp), save, public :: z_depot_hi     = 48.0e3_dp   ! depot/haze band top = cloud base
  real(dp), save, public :: z_domain_top   = 70.0e3_dp   ! above cloud -> lost (DEAD)

  public :: set_rho_spore
  public :: seed_spores, transport_step, transport_test, bio_cleanup
  public :: settling_check

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

  ! (Re)allocate the organism store to n slots, all DEAD.
  subroutine org_alloc(n)
    integer, intent(in) :: n
    call bio_cleanup()
    max_orgs = n
    allocate(o_state(n), o_z(n), o_rcell(n), o_age(n))
    o_state = ST_DEAD; o_z = 0.0_dp; o_rcell = 0.0_dp; o_age = 0.0_dp
    n_orgs = 0; n_lost_top = 0; n_lost_bot = 0
  end subroutine org_alloc

  subroutine bio_cleanup()
    if (allocated(o_state)) deallocate(o_state)
    if (allocated(o_z))     deallocate(o_z)
    if (allocated(o_rcell)) deallocate(o_rcell)
    if (allocated(o_age))   deallocate(o_age)
    max_orgs = 0; n_orgs = 0
  end subroutine bio_cleanup

  ! Seed n DORMANT spores of radius rcell, uniformly through the depot/haze band.
  subroutine seed_spores(n, rcell)
    integer,  intent(in) :: n
    real(dp), intent(in) :: rcell
    integer  :: i
    real(dp) :: u
    call org_alloc(n)
    n_orgs = n
    do i = 1, n
      o_state(i) = ST_DORMANT
      call random_number(u)
      o_z(i)     = z_depot_lo + u * (z_depot_hi - z_depot_lo)
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

    call seed_spores(n, rcell)
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

end module bio_venus
