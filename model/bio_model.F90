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

  ! ---- Biomass field ----
  real(8), allocatable, save :: biomass(:)     ! [kg] per atmosphere level

  ! ---- Run parameters (set in bio_init_run) ----
  real(8), save :: v_conv_ms     ! convective velocity [m/s]
  real(8), save :: halflife_days
  real(8), save :: dt_s          ! timestep [s]
  real(8), save :: growth_rate   ! max fractional growth rate [/day]

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
        return
      end if
    end do
  end subroutine add_organism

  !--------------------------------------------------------------------
  ! Kill organism at index i; return its mass to the biomass pool
  !--------------------------------------------------------------------
  subroutine kill_organism(i)
    integer, intent(in) :: i
    integer :: lev
    lev = atm_level(org_z(i))
    biomass(lev)  = biomass(lev) + org_mass(i)
    org_alive(i)  = .false.
    n_orgs        = n_orgs - 1
  end subroutine kill_organism

  !--------------------------------------------------------------------
  ! Neutral-buoyancy mass for a solid sphere (G=0) at given atmosphere
  ! v_conv = GRAV*(2/9)*R^2*rho_org/(nu*rho_gas)  =>
  ! R_eq = sqrt(v_conv*9*nu*rho_gas / (2*GRAV*rho_org))
  ! m_eq = rho_org*(4pi/3)*R_eq^3
  !--------------------------------------------------------------------
  pure function neutral_buoyancy_mass(v_conv, nu, rho_gas, rho_org) result(m_eq)
    real(8), intent(in) :: v_conv, nu, rho_gas, rho_org
    real(8) :: R_eq, m_eq
    R_eq = sqrt(v_conv * 9.0d0 * nu * rho_gas / (2.0d0 * GRAV * rho_org))
    m_eq = rho_org * (4.0d0/3.0d0) * PI * R_eq**3
  end function neutral_buoyancy_mass

  !--------------------------------------------------------------------
  ! Initialise one ensemble member
  ! b_ref_kg: reference biomass pool size [kg] (independent of n_init)
  ! b_factor: multiplier on b_ref_kg (= 1 for control, = 3 for B-sensitivity)
  !--------------------------------------------------------------------
  subroutine bio_init_run(n_init, m_init, v_conv, halflife, dt_hrs, &
                          b_ref_kg, b_factor, grwth)
    integer, intent(in) :: n_init
    real(8), intent(in) :: m_init, v_conv, halflife, dt_hrs
    real(8), intent(in) :: b_ref_kg, b_factor, grwth
    integer  :: i
    real(8)  :: u, z0, G0, rho0, B_total

    v_conv_ms    = v_conv
    halflife_days = halflife
    dt_s         = dt_hrs * 3600.0d0
    growth_rate  = grwth

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

    ! Biomass: total = b_ref_kg * b_factor, evenly distributed over levels.
    ! b_ref_kg is the atmospheric nutrient reservoir (independent of n_init).
    if (allocated(biomass)) deallocate(biomass)
    allocate(biomass(n_lev))
    B_total  = b_ref_kg * b_factor
    biomass  = B_total / real(n_lev, 8)

    ! Seed initial organisms with random properties throughout the AHZ
    do i = 1, n_init
      call random_number(u);  z0   = u * Z_AHZ_TOP
      call random_number(u);  G0   = G_MIN + u * (G_MAX - G_MIN)
      call random_number(u);  rho0 = RHO_ORG_MIN + u * (RHO_ORG_MAX - RHO_ORG_MIN)
      call add_organism(z0, m_init, G0, rho0, m_init)
    end do
  end subroutine bio_init_run

  !--------------------------------------------------------------------
  ! Advance one timestep
  ! Timestep order: eat -> move -> AHZ check -> reproduce -> age/die
  !--------------------------------------------------------------------
  subroutine bio_step(n_born, n_died)
    integer, intent(out) :: n_born, n_died
    integer  :: i, nc, ic, lev
    real(8)  :: rho_loc, nu_loc
    real(8)  :: v_sed, dz, dB
    real(8)  :: G_c, rho_c, mrepr_c
    real(8)  :: u, p_death, noise

    n_born  = 0
    n_died  = 0

    ! Per-timestep death probability from exponential decay with half-life
    p_death = 1.0d0 - 2.0d0**(-dt_s / (halflife_days * SEC_PER_DAY))

    do i = 1, n_slots
      if (.not. org_alive(i)) cycle

      ! Atmosphere at current (pre-move) position
      lev = atm_level(org_z(i))
      call atm_rhovisc(org_z(i), rho_loc, nu_loc)

      ! -- 1. Eat (growth from biomass pool) --
      ! Each organism consumes at most growth_rate*(dt/day)*own_mass,
      ! and at most half the local pool to avoid monopolising a level
      if (biomass(lev) > 0.0d0) then
        dB = min(growth_rate * (dt_s / SEC_PER_DAY) * org_mass(i), &
                 0.5d0 * biomass(lev))
        org_mass(i)   = org_mass(i) + dB
        biomass(lev)  = biomass(lev) - dB
        org_radius(i) = mass_to_radius(org_mass(i), org_G(i), org_rho(i))
      end if

      ! -- 2. Move (convection + Stokes sedimentation) --
      v_sed    = sedimentation_vel(org_mass(i), org_radius(i), &
                                   org_rho(i), rho_loc, nu_loc)
      dz       = (v_conv_ms + v_sed) * dt_s
      org_z(i) = org_z(i) + dz

      ! -- 3. AHZ boundary check (altitude is equivalent to T-check
      !        for the linear profile used here) --
      if (org_z(i) < 0.0d0 .or. org_z(i) > Z_AHZ_TOP) then
        call kill_organism(i)
        n_died = n_died + 1
        cycle
      end if

      ! -- 4. Reproduce --
      nc = int(org_mass(i) / org_mrepr(i)) - 1
      if (nc >= 1) then
        do ic = 1, nc
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
          n_born = n_born + 1
        end do
        ! Parent retains remainder mass after reproduction
        org_mass(i)   = org_mass(i) - real(nc, 8) * org_mrepr(i)
        org_radius(i) = mass_to_radius(org_mass(i), org_G(i), org_rho(i))
      end if

      ! -- 5. Age and stochastic death --
      org_age(i) = org_age(i) + dt_s / SEC_PER_DAY
      call random_number(u)
      if (u < p_death) then
        call kill_organism(i)
        n_died = n_died + 1
      end if

    end do

  end subroutine bio_step

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
    if (allocated(org_z))      deallocate(org_z)
    if (allocated(org_mass))   deallocate(org_mass)
    if (allocated(org_radius)) deallocate(org_radius)
    if (allocated(org_G))      deallocate(org_G)
    if (allocated(org_rho))    deallocate(org_rho)
    if (allocated(org_mrepr))  deallocate(org_mrepr)
    if (allocated(org_age))    deallocate(org_age)
    if (allocated(org_alive))  deallocate(org_alive)
    if (allocated(biomass))    deallocate(biomass)
  end subroutine bio_cleanup

end module bio_model
