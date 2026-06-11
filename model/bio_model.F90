!======================================================================
! bio_model.F90
! Individual-based organism lifecycle model following Yates et al. (2017)
!
! Each organism is a hollow spherical shell characterised by:
!   radius R  [m]       -- outer radius
!   G         [-]       -- growth strategy: inner_radius = G*R
!                          G -> 1 = balloon (thin shell)
!                          G -> 0 = solid sphere
!   rho_org   [kg/m^3]  -- organic skin density
!   m_repr    [kg]      -- reproduction mass threshold
!   age       [days]
!
! Volume of organic material: V_org = (4pi/3) R^3 (1 - G^3)
! Organism mass:              m = V_org * rho_org
!
! Vertical motion: dz/dt = v_conv + v_sed
!   v_conv  [m/s]  -- upward convective velocity (constant, positive)
!   v_sed   [m/s]  -- Stokes sedimentation (negative = downward)
!              v_sed = -g * V_org * (rho_org - rho_gas) / (6 pi nu rho_gas R)
!
! Timestep order (per Yates Sect. 2.2): eat -> move -> AHZ check ->
!                                        reproduce -> age/die
!
! Biomass bookkeeping (the pool is not conserved — it has sources and sinks):
!   - eating moves biomass from a level into the organism;
!   - a half-life (old-age) death returns the organism's mass to the biomass of
!     the level where it died (local return);
!   - an organism convected out the TOP of the AHZ is removed and its mass is lost
!     (carried out of the domain by the updraft — the top sink);
!   - an organism that sinks out the BOTTOM of the AHZ likewise carries its mass
!     out and it is lost;
!   - the biomass field is ADVECTED UPWARD by the convective updraft (v_conv):
!     food enters through the bottom boundary at biomass_flux [kg/day], is carried
!     up through the AHZ, and is lost through the top (see advect_biomass).
!
! Growth is BIOMASS-LIMITED (Yates Sect. 2.2): an organism's growth is limited
! only by the biomass available in its level, not by a fixed maximum rate.  The
! biomass in each level is shared among the organisms there in proportion to
! their mass (bigger organisms eat more).  An optional maximum specific growth
! rate (growth_rate, [/day]) can be re-enabled by setting it > 0; Yates uses no
! such cap, so the default is 0 (disabled).
!
! Units: SI throughout (m, kg, s, K)
!======================================================================
module bio_model
  use bio_params
  use bio_atmosphere
  implicit none

  ! ---- Organism state arrays (index 1..n_slots) ----
  real(8), allocatable, save :: org_z(:)       ! altitude [m]
  real(8), allocatable, save :: org_mass(:)    ! [kg]
  real(8), allocatable, save :: org_radius(:)  ! outer radius R [m]
  real(8), allocatable, save :: org_G(:)       ! growth strategy [-]
  real(8), allocatable, save :: org_rho(:)     ! skin density [kg/m^3]
  real(8), allocatable, save :: org_mrepr(:)   ! reproduction mass [kg]
  real(8), allocatable, save :: org_age(:)     ! age [days]
  logical, allocatable, save :: org_alive(:)

  integer, save :: n_slots    ! allocated array length
  integer, save :: n_orgs     ! current live count
  integer, save :: n_free     ! number of vacant slots (n_slots - n_orgs)

  ! ---- Biomass field ----
  real(8), allocatable, save :: biomass(:)     ! [kg] per atmosphere level

  ! ---- Per-level scratch for biomass-limited growth (frozen at start of step) ----
  real(8), allocatable, save :: mass_lev_start(:) ! organism mass per level, pre-eat [kg]
  real(8), allocatable, save :: biomass_avail(:)  ! biomass per level at step start [kg]

  ! ---- Run parameters (set in bio_init_run) ----
  real(8), save :: v_conv_ms     ! convective velocity [m/s]
  real(8), save :: halflife_days
  real(8), save :: dt_s          ! timestep [s]
  real(8), save :: growth_rate   ! max fractional growth rate [/day]
  ! private: avoids a name clash with the driver's namelist variable of the same name
  real(8), save, private :: biomass_flux  ! tunable food flux into bottom level [kg/day]

contains

  !--------------------------------------------------------------------
  ! Compute outer radius from mass, G, and skin density
  !--------------------------------------------------------------------
  pure function mass_to_radius(m, G, rho) result(R)
    real(8), intent(in) :: m, G, rho
    real(8) :: R
    real(8) :: shell_vol_coeff
    ! m = rho * (4pi/3) R^3 (1 - G^3)  =>  R = (m/(rho*(4pi/3)*(1-G^3)))^(1/3)
    shell_vol_coeff = (4.0d0/3.0d0) * PI * (1.0d0 - G**3)
    R = (m / (rho * shell_vol_coeff))**(1.0d0/3.0d0)
  end function mass_to_radius

  !--------------------------------------------------------------------
  ! Stokes sedimentation velocity (negative = downward)
  ! v_sed = -g * V_org * (rho_org - rho_gas) / (6 pi nu rho_gas R)
  !--------------------------------------------------------------------
  pure function sedimentation_vel(m, R, rho_org, rho_gas, nu) result(v_sed)
    real(8), intent(in) :: m, R, rho_org, rho_gas, nu
    real(8) :: v_sed, V_org
    V_org = m / rho_org
    v_sed = -GRAV * V_org * (rho_org - rho_gas) / &
            (6.0d0 * PI * nu * rho_gas * R)
  end function sedimentation_vel

  !--------------------------------------------------------------------
  ! Clamp x to [lo, hi]
  !--------------------------------------------------------------------
  pure function clamp(x, lo, hi) result(y)
    real(8), intent(in) :: x, lo, hi
    real(8) :: y
    y = max(lo, min(hi, x))
  end function clamp

  !--------------------------------------------------------------------
  ! Box-Muller normal random variate
  !--------------------------------------------------------------------
  subroutine normal_rand(z)
    real(8), intent(out) :: z
    real(8) :: u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, 1.0d-15)
    z  = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * PI * u2)
  end subroutine normal_rand

  !--------------------------------------------------------------------
  ! Place a new organism into the first available (dead) slot
  ! If no slot is free the organism is discarded (population at cap)
  !--------------------------------------------------------------------
  subroutine add_organism(z, mass, G, rho, mrepr)
    real(8), intent(in) :: z, mass, G, rho, mrepr
    integer :: i
    if (n_free == 0) return   ! population at cap; discard offspring O(1)
    do i = 1, n_slots
      if (.not. org_alive(i)) then
        org_z(i)      = z
        org_mass(i)   = mass
        org_radius(i) = mass_to_radius(mass, G, rho)
        org_G(i)      = G
        org_rho(i)    = rho
        org_mrepr(i)  = mrepr
        org_age(i)    = 0.0d0
        org_alive(i)  = .true.
        n_orgs        = n_orgs + 1
        n_free        = n_free - 1
        return
      end if
    end do
  end subroutine add_organism

  !--------------------------------------------------------------------
  ! Kill organism at index i (remove it from the population).
  ! Biomass return is handled by the caller: a half-life death deposits the
  ! organism's mass into its level's biomass; a washout death (left the AHZ)
  ! returns nothing — that biomass is lost from the pool.
  !--------------------------------------------------------------------
  subroutine kill_organism(i)
    integer, intent(in) :: i
    org_alive(i)  = .false.
    n_orgs        = n_orgs - 1
    n_free        = n_free + 1
  end subroutine kill_organism

  !--------------------------------------------------------------------
  ! Neutral-buoyancy mass for a solid sphere (G=0) at given atmosphere
  !--------------------------------------------------------------------
  pure function neutral_buoyancy_mass(v_conv, nu, rho_gas, rho_org) result(m_eq)
    real(8), intent(in) :: v_conv, nu, rho_gas, rho_org
    real(8) :: R_eq, m_eq
    R_eq = sqrt(v_conv * 9.0d0 * nu * rho_gas / (2.0d0 * GRAV * rho_org))
    m_eq = rho_org * (4.0d0/3.0d0) * PI * R_eq**3
  end function neutral_buoyancy_mass

  !--------------------------------------------------------------------
  ! Initialise one ensemble member
  !--------------------------------------------------------------------
  subroutine bio_init_run(n_init, m_init, v_conv, halflife, dt_hrs, &
                          b_ref_kg, b_factor, grwth, mrepr_seed_max, bflux)
    integer, intent(in) :: n_init
    real(8), intent(in) :: m_init, v_conv, halflife, dt_hrs
    real(8), intent(in) :: b_ref_kg, b_factor, grwth, mrepr_seed_max, bflux
    integer  :: i
    real(8)  :: u, z0, G0, rho0, B_total, m0

    v_conv_ms    = v_conv
    halflife_days = halflife
    dt_s         = dt_hrs * 3600.0d0
    growth_rate  = grwth
    biomass_flux = bflux

    ! Allocate / reset organism arrays
    n_slots = MAX_ORGS
    if (allocated(org_z)) &
      deallocate(org_z, org_mass, org_radius, org_G, org_rho, &
                 org_mrepr, org_age, org_alive)
    allocate(org_z(n_slots), org_mass(n_slots), org_radius(n_slots), &
             org_G(n_slots), org_rho(n_slots), org_mrepr(n_slots),   &
             org_age(n_slots), org_alive(n_slots))
    org_alive = .false.
    n_orgs    = 0
    n_free    = n_slots

    ! Biomass pool: initial total = b_ref_kg * b_factor, spread evenly.  The pool
    ! is not conserved — boundary-exit deaths remove biomass, the bottom flux adds
    ! it, and it is advected upward and lost at the top (see bio_step).
    if (allocated(biomass))        deallocate(biomass)
    if (allocated(mass_lev_start)) deallocate(mass_lev_start)
    if (allocated(biomass_avail))  deallocate(biomass_avail)
    allocate(biomass(n_lev), mass_lev_start(n_lev), biomass_avail(n_lev))
    B_total  = b_ref_kg * b_factor
    biomass  = B_total / real(n_lev, 8)
    mass_lev_start = 0.0d0
    biomass_avail  = 0.0d0

    ! Seed initial organisms with random properties throughout the AHZ.
    ! Founder mass / reproduction mass: by default all founders start at m_init
    ! (mrepr_seed_max <= m_init).  If mrepr_seed_max > m_init, the founder mass is
    ! drawn log-uniformly over [m_init, mrepr_seed_max] and m_repr is set to that
    ! mass (Yates: "reproduction mass is close to the birth mass").  This tests
    ! whether seeding some founders near the neutral-buoyancy mass lets a floating
    ! sub-population establish from a cold start.
    do i = 1, n_init
      call random_number(u);  z0   = u * Z_AHZ_TOP
      call random_number(u);  G0   = G_MIN + u * (G_MAX - G_MIN)
      call random_number(u);  rho0 = RHO_ORG_MIN + u * (RHO_ORG_MAX - RHO_ORG_MIN)
      if (mrepr_seed_max > m_init) then
        call random_number(u)
        m0 = m_init * (mrepr_seed_max / m_init)**u   ! log-uniform in [m_init, mrepr_seed_max]
      else
        m0 = m_init
      end if
      call add_organism(z0, m0, G0, rho0, m0)
    end do
  end subroutine bio_init_run

  !--------------------------------------------------------------------
  ! Advance one timestep
  ! Timestep order: eat -> move -> AHZ check -> reproduce -> age/die
  !--------------------------------------------------------------------
  subroutine bio_step(n_born, n_died)
    integer, intent(out) :: n_born, n_died
    integer  :: i, nc, ic, nc_made, lev
    real(8)  :: rho_loc, nu_loc
    real(8)  :: v_sed, dz, dB
    real(8)  :: G_c, rho_c, mrepr_c
    real(8)  :: u, p_death, noise

    n_born = 0
    n_died = 0

    ! Per-timestep death probability from exponential decay with half-life
    p_death = 1.0d0 - 2.0d0**(-dt_s / (halflife_days * SEC_PER_DAY))

    ! Pre-pass: freeze each level's organism mass and biomass for biomass-limited
    ! growth.  Growth shares a level's biomass among its organisms by mass, so we
    ! need the level totals as they stand at the start of the step (before any
    ! eating, moving, births, or deaths).  Using frozen values makes each
    ! organism's share independent of processing order.
    mass_lev_start = 0.0d0
    do i = 1, n_slots
      if (.not. org_alive(i)) cycle
      lev = atm_level(org_z(i))
      mass_lev_start(lev) = mass_lev_start(lev) + org_mass(i)
    end do
    biomass_avail = biomass

    do i = 1, n_slots
      if (.not. org_alive(i)) cycle

      ! Atmosphere at current (pre-move) position
      lev = atm_level(org_z(i))
      call atm_rhovisc(org_z(i), rho_loc, nu_loc)

      ! -- 1. Eat: biomass-limited growth (Yates 2017 Sect. 2.2) --
      !    Growth is limited only by available biomass, not by a fixed rate.
      !    The organism consumes a share of its level's biomass in proportion to
      !    its mass.  The live min(., biomass(lev)) guard keeps biomass >= 0 and
      !    stops organisms born this step (not in mass_lev_start) from consuming
      !    biomass their parent already took.  An optional max specific growth
      !    rate is applied only if growth_rate > 0 (Yates uses none).
      if (mass_lev_start(lev) > 0.0d0 .and. biomass_avail(lev) > 0.0d0) then
        dB = biomass_avail(lev) * org_mass(i) / mass_lev_start(lev)
        if (growth_rate > 0.0d0) &
          dB = min(dB, growth_rate * (dt_s / SEC_PER_DAY) * org_mass(i))
        dB = min(dB, biomass(lev))
        if (dB > 0.0d0) then
          org_mass(i)   = org_mass(i) + dB
          biomass(lev)  = biomass(lev) - dB
          org_radius(i) = mass_to_radius(org_mass(i), org_G(i), org_rho(i))
        end if
      end if

      ! -- 2. Move: constant upward convection + sedimentation (Yates 2017 Sect. 2.2) --
      v_sed    = sedimentation_vel(org_mass(i), org_radius(i), &
                                   org_rho(i), rho_loc, nu_loc)
      dz       = (v_conv_ms + v_sed) * dt_s
      org_z(i) = org_z(i) + dz

      ! -- 3. AHZ boundary check --
      !    An organism that leaves the AHZ at either boundary is removed and its
      !    mass is lost from the domain (top: carried up and out by the updraft;
      !    bottom: sinks out below).
      if (org_z(i) < 0.0d0 .or. org_z(i) > Z_AHZ_TOP) then
        call kill_organism(i)
        n_died = n_died + 1
        cycle
      end if

      ! -- 4. Reproduce --
      !    Progeny are born at the parent's reproduction mass with mutated traits.
      !    Stop if the population is at the cap (n_free == 0) and debit the parent
      !    only for progeny actually created, so biomass stays conserved even when
      !    offspring are discarded at the cap.
      nc = int(org_mass(i) / org_mrepr(i)) - 1
      if (nc >= 1) then
        nc_made = 0
        do ic = 1, nc
          if (n_free == 0) exit
          G_c     = org_G(i)
          rho_c   = org_rho(i)
          mrepr_c = org_mrepr(i)
          call normal_rand(noise)
          G_c     = clamp(G_c + SIG_G * noise, G_MIN, G_MAX)
          call normal_rand(noise)
          rho_c   = clamp(rho_c + SIG_RHO * noise, RHO_ORG_MIN, RHO_ORG_MAX)
          call normal_rand(noise)
          mrepr_c = mrepr_c * exp(SIG_MREPR * noise)
          call add_organism(org_z(i), org_mrepr(i), G_c, rho_c, mrepr_c)
          nc_made = nc_made + 1
          n_born  = n_born + 1
        end do
        if (nc_made > 0) then
          ! Parent retains remainder mass after reproduction
          org_mass(i)   = org_mass(i) - real(nc_made, 8) * org_mrepr(i)
          org_radius(i) = mass_to_radius(org_mass(i), org_G(i), org_rho(i))
        end if
      end if

      ! -- 5. Age and stochastic death --
      !    A half-life death deposits the organism's mass into the biomass of the
      !    level where it died (local return).  Boundary-exit deaths are handled in
      !    step 3 (top exit returned at end of step; bottom exit lost).
      org_age(i) = org_age(i) + dt_s / SEC_PER_DAY
      call random_number(u)
      if (u < p_death) then
        lev = atm_level(org_z(i))
        biomass(lev) = biomass(lev) + org_mass(i)
        call kill_organism(i)
        n_died = n_died + 1
        cycle
      end if

    end do

    ! -- 6. Advect the biomass field upward; inject the bottom flux, lose at top --
    call advect_biomass()

  end subroutine bio_step

  !--------------------------------------------------------------------
  ! Advect the biomass field upward with the convective updraft v_conv.
  ! Finite-volume first-order upwind, sub-stepped to satisfy the CFL
  ! condition (the per-step Courant number is large, ~20, because the
  ! updraft crosses many levels per timestep).  Food enters through the
  ! bottom boundary at biomass_flux [kg/day] and is carried out (lost)
  ! through the top boundary.  Levels run 1 = bottom .. n_lev = top.
  !--------------------------------------------------------------------
  subroutine advect_biomass()
    integer :: k, isub, n_sub
    real(8) :: dz_lev, courant, c_sub, dt_sub, inflow
    real(8) :: flux_up(n_lev)

    if (v_conv_ms <= 0.0d0) then
      ! No updraft: the bottom source (if any) simply accumulates in level 1.
      if (biomass_flux > 0.0d0) &
        biomass(1) = biomass(1) + biomass_flux * (dt_s / SEC_PER_DAY)
      return
    end if

    dz_lev  = Z_AHZ_TOP / real(n_lev - 1, 8)
    courant = v_conv_ms * dt_s / dz_lev
    n_sub   = max(1, ceiling(courant))
    c_sub   = courant / real(n_sub, 8)               ! <= 1 (CFL-stable)
    dt_sub  = dt_s / real(n_sub, 8)
    inflow  = biomass_flux * (dt_sub / SEC_PER_DAY)   ! bottom inflow per sub-step [kg]

    do isub = 1, n_sub
      ! Upwind: fraction c_sub of each cell's biomass crosses its upper face.
      do k = 1, n_lev
        flux_up(k) = c_sub * biomass(k)
      end do
      ! Each cell loses its upward flux and gains the cell-below's flux.
      ! flux_up(n_lev) exits the top and is lost; inflow enters at the bottom.
      do k = n_lev, 2, -1
        biomass(k) = biomass(k) - flux_up(k) + flux_up(k-1)
      end do
      biomass(1) = biomass(1) - flux_up(1) + inflow
    end do
  end subroutine advect_biomass

  !--------------------------------------------------------------------
  ! Write current organism state to an open file unit (one line/organism)
  ! Columns: z[m] mass[kg] radius[m] G rho_org[kg/m3] age[days] skin_width[m]
  !--------------------------------------------------------------------
  subroutine bio_write_state(iunit)
    integer, intent(in) :: iunit
    integer :: i
    do i = 1, n_slots
      if (.not. org_alive(i)) cycle
      write(iunit, '(7(1pe16.8,1x))') &
        org_z(i), org_mass(i), org_radius(i), org_G(i), &
        org_rho(i), org_age(i), org_radius(i) * (1.0d0 - org_G(i))
    end do
  end subroutine bio_write_state

  subroutine bio_cleanup()
    if (allocated(org_z))        deallocate(org_z)
    if (allocated(org_mass))     deallocate(org_mass)
    if (allocated(org_radius))   deallocate(org_radius)
    if (allocated(org_G))        deallocate(org_G)
    if (allocated(org_rho))      deallocate(org_rho)
    if (allocated(org_mrepr))    deallocate(org_mrepr)
    if (allocated(org_age))      deallocate(org_age)
    if (allocated(org_alive))    deallocate(org_alive)
    if (allocated(biomass))        deallocate(biomass)
    if (allocated(mass_lev_start)) deallocate(mass_lev_start)
    if (allocated(biomass_avail))  deallocate(biomass_avail)
  end subroutine bio_cleanup

end module bio_model
