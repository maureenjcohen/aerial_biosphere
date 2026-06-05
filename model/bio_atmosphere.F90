!======================================================================
! bio_atmosphere.F90
! Parametric 1-D atmosphere for the Y dwarf AHZ
!
! Uses linear profiles in altitude between the AHZ boundaries,
! consistent with the density and viscosity ranges given in
! Yates et al. (2017) Sect. 3.1.  A more accurate profile
! (e.g. from Morley et al. 2014) can be substituted by replacing
! atm_init with a file-reading version; the interpolation interface
! is unchanged.
!
! Units: SI throughout
!   altitude  [m]
!   temperature [K]
!   density   [kg/m^3]
!   kinematic viscosity [m^2/s]
!======================================================================
module bio_atmosphere
  use bio_params
  implicit none

  integer,  save              :: n_lev
  real(8), allocatable, save  :: z_lev(:)     ! altitude [m]
  real(8), allocatable, save  :: T_lev(:)     ! temperature [K]
  real(8), allocatable, save  :: rho_lev(:)   ! gas density [kg/m^3]
  real(8), allocatable, save  :: nu_lev(:)    ! kinematic viscosity [m^2/s]

contains

  !--------------------------------------------------------------------
  ! Initialise the atmospheric profile on n grid levels
  ! (linearly spaced from 0 to Z_AHZ_TOP)
  !--------------------------------------------------------------------
  subroutine atm_init(n)
    integer, intent(in) :: n
    integer  :: k
    real(8)  :: frac

    n_lev = n
    allocate(z_lev(n_lev), T_lev(n_lev), rho_lev(n_lev), nu_lev(n_lev))

    do k = 1, n_lev
      frac        = real(k-1, 8) / real(n_lev-1, 8)   ! 0 = bottom, 1 = top
      z_lev(k)   = frac * Z_AHZ_TOP
      T_lev(k)   = T_AHZ_BOT + frac * (T_AHZ_TOP - T_AHZ_BOT)
      rho_lev(k) = RHO_GAS_BOT + frac * (RHO_GAS_TOP - RHO_GAS_BOT)
      nu_lev(k)  = NU_BOT      + frac * (NU_TOP      - NU_BOT)
    end do
  end subroutine atm_init

  !--------------------------------------------------------------------
  ! Interpolate atmosphere at altitude z [m]
  ! Returns T [K], rho_gas [kg/m^3], nu [m^2/s]
  ! Clamps to profile limits if z is outside [0, Z_AHZ_TOP]
  !--------------------------------------------------------------------
  subroutine atm_interp(z, T_loc, rho_loc, nu_loc)
    real(8), intent(in)  :: z
    real(8), intent(out) :: T_loc, rho_loc, nu_loc
    real(8) :: frac

    frac    = max(0.0d0, min(1.0d0, z / Z_AHZ_TOP))
    T_loc   = T_AHZ_BOT + frac * (T_AHZ_TOP - T_AHZ_BOT)
    rho_loc = RHO_GAS_BOT + frac * (RHO_GAS_TOP - RHO_GAS_BOT)
    nu_loc  = NU_BOT      + frac * (NU_TOP      - NU_BOT)
  end subroutine atm_interp

  !--------------------------------------------------------------------
  ! Return the level index (1-based) for altitude z [m]
  !--------------------------------------------------------------------
  pure function atm_level(z) result(k)
    real(8), intent(in) :: z
    integer :: k
    k = max(1, min(n_lev, int(z / Z_AHZ_TOP * real(n_lev-1, 8)) + 1))
  end function atm_level

  !--------------------------------------------------------------------
  ! Convenience: return only rho_gas and nu at altitude z [m]
  !--------------------------------------------------------------------
  subroutine atm_rhovisc(z, rho_loc, nu_loc)
    real(8), intent(in)  :: z
    real(8), intent(out) :: rho_loc, nu_loc
    real(8) :: frac
    frac    = max(0.0d0, min(1.0d0, z / Z_AHZ_TOP))
    rho_loc = RHO_GAS_BOT + frac * (RHO_GAS_TOP - RHO_GAS_BOT)
    nu_loc  = NU_BOT      + frac * (NU_TOP      - NU_BOT)
  end subroutine atm_rhovisc

  subroutine atm_cleanup()
    if (allocated(z_lev))   deallocate(z_lev)
    if (allocated(T_lev))   deallocate(T_lev)
    if (allocated(rho_lev)) deallocate(rho_lev)
    if (allocated(nu_lev))  deallocate(nu_lev)
  end subroutine atm_cleanup

end module bio_atmosphere
