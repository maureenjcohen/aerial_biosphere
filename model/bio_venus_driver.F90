!======================================================================
! bio_venus_driver.F90
! Run driver for the Venus cloud biosphere model.
! Reads an optional namelist (default bio_cloud.nml, or argv(1)):
!   &cloud_run  indir='../inputs'  src_mode1=1 src_mode2=1 src_mode3=0
!               phi=0.6  do_transport=.true.  n_spore=5000  rcell_um=0.5
!               dt_s=5000.0  nsteps=2000  nout=200 /
!   (source ids: 0 off, 1 vpcm, 2 haus, 3 kh80)
! Default: VPCM modes 1&2, mode 3 off.
!======================================================================
program bio_venus_driver
  use bio_cloud_model
  use bio_venus
  implicit none
  character(len=256) :: indir, logfile, nmlfile
  integer :: ulog, istat
  integer :: src_mode1, src_mode2, src_mode3
  real(8) :: phi
  logical :: do_transport, do_lifecycle, do_evolve, trait_hist
  integer :: n_spore, nsteps, nout, maxorgs
  real(8) :: rcell_um, dt_s, Ydays, seed_lo_km, seed_hi_km
  real(8) :: Aref_m2, X0days, cellhalf_days, hosthalf_days, rhocell_wet
  integer :: rng_seed
  logical :: do_conv
  real(8) :: conv_sigma_ms, conv_tau_min, conv_zlo_km, conv_zhi_km, conv_edge_km
  integer :: conv_nsub_n
  namelist /cloud_run/ indir, src_mode1, src_mode2, src_mode3, phi, &
       do_transport, do_lifecycle, do_evolve, n_spore, rcell_um, dt_s, nsteps, nout, &
       Ydays, seed_lo_km, seed_hi_km, maxorgs, Aref_m2, X0days, cellhalf_days, hosthalf_days, &
       trait_hist, rhocell_wet, rng_seed, &
       do_conv, conv_sigma_ms, conv_tau_min, conv_zlo_km, conv_zhi_km, conv_edge_km, conv_nsub_n

  indir   = '../inputs'
  logfile = 'cloud_read.log'
  src_mode1 = SRC_VPCM
  src_mode2 = SRC_VPCM
  src_mode3 = SRC_OFF
  phi       = PHI_PACK_DEFAULT
  do_transport = .false.        ! step-1 kinematics regression (off by default)
  do_lifecycle = .false.        ! step-2 state machine (no fission)
  do_evolve    = .true.         ! step-3 fission + mutation + carrying capacity
  n_spore   = 5000
  rcell_um  = 0.5d0
  dt_s      = 5000.0d0
  nsteps    = 2000
  nout      = 200
  Ydays     = 10.0d0            ! dormancy survival time [Earth-days]
  seed_lo_km = 33.0d0           ! seed band: the depot/haze (plan IC) by default
  seed_hi_km = 48.0d0
  maxorgs   = 300000            ! organism-store capacity (births reuse slots)
  Aref_m2   = 1.0d-10           ! reference area -> population scale
  X0days    = 1.0d0             ! seed reproduction half-life [Earth-days]
  cellhalf_days = 10.0d0        ! ACTIVE-cell baseline half-life [Earth-days]
  hosthalf_days = 5.0d0         ! host-droplet lifetime [Earth-days]
  trait_hist    = .false.       ! dump per-snapshot ACTIVE trait histograms (traits_*.csv)
  rhocell_wet   = 1050.0d0      ! hydrated in-droplet cell density [kg/m3] (settling verdict)
  rng_seed      = 0             ! >0 seeds the RNG (ensemble member id); 0 = default sequence
  do_conv       = .false.       ! layer Vega convective up/downdrafts on the large-scale wind
  conv_sigma_ms = 0.6d0         ! convective vertical-wind std [m/s] (Vega W_a)
  conv_tau_min  = 20.0d0        ! convective correlation time [minutes]
  conv_zlo_km   = 48.0d0        ! convective layer bottom [km] (= depot top; decisive knob)
  conv_zhi_km   = 55.0d0        ! convective layer top [km]
  conv_edge_km  = 2.0d0         ! taper half-width at each edge [km]
  conv_nsub_n   = 20            ! transport sub-steps per biology step

  nmlfile = 'bio_cloud.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, nmlfile)
  open(newunit=ulog, file=trim(nmlfile), status='old', iostat=istat)
  if (istat == 0) then
    read(ulog, nml=cloud_run, iostat=istat)
    close(ulog)
    if (istat /= 0) write(*,'(a)') 'Warning: error reading bio_cloud.nml; using defaults.'
  end if

  write(*,'(a)') '============================================'
  write(*,'(a)') ' Venus cloud biosphere model'
  write(*,'(a)') '============================================'
  write(*,'(a)') ' input dir : '//trim(indir)
  write(*,'(a,3(1x,a))') ' sources   : mode1=', src_name(src_mode1), &
       ' mode2='//src_name(src_mode2), ' mode3='//src_name(src_mode3)

  call cloud_read(trim(indir))
  call set_source(1, src_mode1)
  call set_source(2, src_mode2)
  call set_source(3, src_mode3)
  call cloud_assemble(trim(indir))
  call cloud_build_spectrum()
  call set_phi(phi)
  call bio_init()
  A_ref     = Aref_m2
  X_init_s  = X0days * 86400.0d0
  cell_half_s = cellhalf_days * 86400.0d0
  host_half_s = hosthalf_days * 86400.0d0
  trait_dump  = trait_hist
  rho_cell_wet = rhocell_wet
  if (rng_seed > 0) call set_rng_seed(rng_seed)
  do_convection = do_conv       ! convection parameterization (bio_cloud_model)
  conv_sigma    = conv_sigma_ms
  conv_tau_s    = conv_tau_min * 60.0d0
  conv_z_lo     = conv_zlo_km  * 1.0d3
  conv_z_hi     = conv_zhi_km  * 1.0d3
  conv_edge     = conv_edge_km * 1.0d3
  conv_nsub     = conv_nsub_n

  open(newunit=ulog, file=trim(logfile), status='replace', action='write')
  call cloud_diagnostics(ulog)
  call cloud_summary(ulog)
  call cloud_spectrum_check(ulog)
  call host_capacity_sweep(ulog)
  call settling_check(ulog)
  if (do_transport) call transport_test(ulog, n_spore, rcell_um*1.0d-6, dt_s, nsteps, nout)
  if (do_lifecycle) call lifecycle_test(ulog, n_spore, rcell_um*1.0d-6, dt_s, nsteps, nout, &
       Ydays, seed_lo_km*1.0d3, seed_hi_km*1.0d3)
  close(ulog)
  ! (evolve_test is heavy: run it once, to stdout only)

  call cloud_summary(6)
  call cloud_spectrum_check(6)
  call host_capacity_sweep(6)
  call settling_check(6)
  if (do_transport) call transport_test(6, n_spore, rcell_um*1.0d-6, dt_s, nsteps, nout)
  if (do_lifecycle) call lifecycle_test(6, n_spore, rcell_um*1.0d-6, dt_s, nsteps, nout, &
       Ydays, seed_lo_km*1.0d3, seed_hi_km*1.0d3)
  if (do_evolve) call evolve_test(6, n_spore, rcell_um*1.0d-6, dt_s, nsteps, nout, &
       Ydays, seed_lo_km*1.0d3, seed_hi_km*1.0d3, maxorgs)
  write(*,'(a)') ' Done.  Full per-level table in '//trim(logfile)//'.'

  call bio_cleanup()
  call cloud_cleanup()
end program bio_venus_driver
