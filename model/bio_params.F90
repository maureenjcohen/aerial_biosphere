!======================================================================
! bio_params.F90
! Physical constants and default parameters for the Yates 2017 IBM
! All quantities in SI (m, kg, s, K) unless noted
!======================================================================
module bio_params
  implicit none

  ! ---- Mathematical / physical constants ----
  real(8), parameter :: PI          = 3.14159265358979d0
  real(8), parameter :: GRAV        = 166.70d0          ! gravitational acceleration [m/s^2]
  real(8), parameter :: R_GAS       = 8.314d0           ! universal gas constant [J/(mol K)]
  real(8), parameter :: M_MEAN      = 2.314d-3          ! mean mol. mass, 85% H2 + 15% He [kg/mol]
  real(8), parameter :: SEC_PER_DAY = 86400.0d0         ! [s]

  ! ---- AHZ temperature boundaries ----
  ! Organisms that stray outside this range die immediately
  real(8), parameter :: T_AHZ_BOT  = 395.0d0            ! hot/lower boundary [K]
  real(8), parameter :: T_AHZ_TOP  = 258.0d0            ! cold/upper boundary [K]
  real(8), parameter :: Z_AHZ_TOP  = 1.05d5             ! AHZ depth [m]  (105 km)

  ! ---- Atmospheric profile parameters ----
  ! Linear profiles between AHZ boundaries (Yates 2017 Sect. 3.1 endpoint values;
  ! derived from Morley 2014 Teff=200K, log_g=5.0 atmosphere via ideal gas law)
  real(8), parameter :: RHO_GAS_BOT = 1.2d0             ! gas density at bottom [kg/m^3]
  real(8), parameter :: RHO_GAS_TOP = 0.4d0             ! gas density at top [kg/m^3]
  real(8), parameter :: NU_BOT      = 1.0d-5            ! kinematic viscosity at bottom [m^2/s]
  real(8), parameter :: NU_TOP      = 2.0d-5            ! kinematic viscosity at top [m^2/s]

  ! ---- Organism property bounds ----
  real(8), parameter :: RHO_ORG_MIN = 500.0d0           ! skin density lower bound [kg/m^3]
  real(8), parameter :: RHO_ORG_MAX = 1500.0d0          ! skin density upper bound [kg/m^3]
  real(8), parameter :: G_MIN       = 0.01d0            ! growth strategy lower bound
  real(8), parameter :: G_MAX       = 0.99d0            ! growth strategy upper bound

  ! ---- Mutation standard deviations ----
  real(8), parameter :: SIG_G       = 0.05d0            ! std dev for G mutation
  real(8), parameter :: SIG_RHO     = 50.0d0            ! std dev for rho_org mutation [kg/m^3]
  real(8), parameter :: SIG_MREPR   = 0.10d0            ! fractional std dev for m_repr (log-normal)

  ! ---- Maximum population size ----
  integer, parameter :: MAX_ORGS    = 50000

  ! ---- Default simulation parameters (overridable by namelist) ----
  integer, parameter :: N_ENS_DEF   = 20                ! ensemble members
  integer, parameter :: N_INIT_DEF  = 100               ! initial organism count
  real(8), parameter :: DT_HRS_DEF  = 6.0d0             ! timestep [hours]
  real(8), parameter :: T_YRS_DEF   = 100.0d0           ! simulation length [Earth years]
  real(8), parameter :: VCONV_DEF   = 1.0d0             ! convective velocity [m/s]
  real(8), parameter :: HL_DEF      = 30.0d0            ! organism half-life [days]
  ! Cold start: Yates Table 1 control run uses m_init = 1e-9 g = 1e-12 kg.  Establishment
  ! from this seed requires a finite max growth rate (GRWTH_DEF below); see bio_model.F90.
  real(8), parameter :: MINIT_DEF   = 1.0d-12          ! initial organism mass [kg]
  real(8), parameter :: BFAC_DEF    = 1.0d0             ! biomass scaling factor
  ! Total conserved biomass pool (Yates Sect. 2.2): organisms consume from it,
  ! and return their mass on death.  No external renewal.
  real(8), parameter :: B_REF_KG   = 1.0d-6            ! total biomass pool [kg]
  ! Max specific growth rate (NPZ mu_max; Yates draws from nutrient-phytoplankton models,
  ! Franks 2002, and reproduces "subject to growth rate").  Growth per step is
  ! min(mu_max*dt*mass, available-biomass share) -> Monod-like.  Must be > 0:
  ! mu_max = 0 lets a lone organism eat a whole level in one step and fall out of the AHZ
  ! before reproducing (cold-start extinction).  ~2-3/day establishes from the 1e-12 seed.
  real(8), parameter :: GRWTH_DEF   = 2.5d0
  ! Founder reproduction-mass seeding spread [kg].  Default <= MINIT_DEF means all
  ! founders start at m_init with m_repr = m_init (literal Yates reading).  If set
  ! > m_init, founder mass / m_repr is drawn log-uniform over [m_init, MRSEED_DEF];
  ! used to test whether seeding some founders near the neutral-buoyancy mass
  ! (m_eq ~ 2.4e-11 kg) lets a floating population establish from a cold start.
  real(8), parameter :: MRSEED_DEF  = 0.0d0
  ! Tunable food flux into the bottom level [kg/day], representing biomass supply
  ! from below the AHZ.  Default 0 (no external supply).  Use to probe how much
  ! bottom flux is required to sustain a population against boundary losses.
  real(8), parameter :: BFLUX_DEF   = 0.0d0

end module bio_params
