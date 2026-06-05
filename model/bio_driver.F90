!======================================================================
! bio_driver.F90
! Standalone ensemble runner for the Yates 2017 IBM
!
! Reads run parameters from bio_run.nml (if present), otherwise uses
! defaults from bio_params.F90.
!
! Output (written to outdir/):
!   ensemble_NNN_state.dat  -- organism properties sampled over the
!                              final simulated year; use for distributions
!   ensemble_NNN_pop.dat    -- population count sampled every 24 h;
!                              use for steady-state (ADF) test
!   run_summary.txt         -- one line per ensemble member
!
! Compile:  see Makefile
! Run:      ./bio_1d            (uses defaults)
!           ./bio_1d bio_run.nml (reads namelist)
!======================================================================
program bio_driver
  use bio_params
  use bio_atmosphere
  use bio_model
  implicit none

  ! ---- Namelist ----
  integer            :: n_ensemble, n_init_orgs, seed_base, n_levels
  real(8)            :: v_conv, halflife, dt_hrs, t_sim_years
  real(8)            :: m_init, b_factor, b_total_kg, growth_rate_day
  character(len=256) :: outdir

  namelist /bio_run/ n_ensemble, n_init_orgs, seed_base, n_levels,   &
                     v_conv, halflife, dt_hrs, t_sim_years,           &
                     m_init, b_factor, b_total_kg, growth_rate_day,   &
                     outdir

  ! ---- Local variables ----
  integer  :: iens, istep, n_steps, n_steps_per_day, n_steps_last_yr
  integer  :: n_born, n_died, istat

  integer  :: seed_size
  integer, allocatable :: seed(:)
  character(len=256)   :: fname, nml_file

  ! ---- Set defaults ----
  n_ensemble      = N_ENS_DEF
  n_init_orgs     = N_INIT_DEF
  seed_base       = 42
  n_levels        = 106
  v_conv          = VCONV_DEF
  halflife        = HL_DEF
  dt_hrs          = DT_HRS_DEF
  t_sim_years     = T_YRS_DEF
  m_init          = MINIT_DEF
  b_factor        = BFAC_DEF
  b_total_kg      = B_REF_KG
  growth_rate_day = GRWTH_DEF
  outdir          = 'output'

  ! ---- Read namelist if provided ----
  nml_file = 'bio_run.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, nml_file)
  open(10, file=trim(nml_file), status='old', iostat=istat)
  if (istat == 0) then
    read(10, nml=bio_run, iostat=istat)
    close(10)
    if (istat /= 0) then
      write(*,'(a)') 'Warning: error reading namelist; using defaults.'
    end if
  end if

  ! ---- Derived quantities ----
  n_steps          = nint(t_sim_years * 365.25d0 * 24.0d0 / dt_hrs)
  n_steps_per_day  = nint(24.0d0 / dt_hrs)
  n_steps_last_yr  = nint(365.25d0 * 24.0d0 / dt_hrs)

  ! ---- Print configuration ----
  write(*,'(a)')        '============================================'
  write(*,'(a)')        ' Yates 2017 IBM — 1-D standalone runner'
  write(*,'(a)')        '============================================'
  write(*,'(a,i6)')     ' Ensemble members : ', n_ensemble
  write(*,'(a,i6)')     ' Initial organisms : ', n_init_orgs
  write(*,'(a,1pe10.3,a)') ' m_init          : ', m_init,         ' kg'
  write(*,'(a,f8.2,a)') ' v_conv           : ', v_conv,          ' m/s'
  write(*,'(a,f8.2,a)') ' Half-life        : ', halflife,        ' days'
  write(*,'(a,f8.2,a)') ' growth_rate_day  : ', growth_rate_day, ' /day'
  write(*,'(a,1pe10.3,a)') ' biomass pool    : ', b_total_kg*b_factor, ' kg'
  write(*,'(a,f8.2,a)') ' Timestep         : ', dt_hrs,          ' hours'
  write(*,'(a,f8.1,a)') ' Sim length       : ', t_sim_years,     ' years'
  write(*,'(a,i8)')     ' Total steps      : ', n_steps
  write(*,'(a,a)')      ' Output dir       : ', trim(outdir)
  write(*,'(a)')        '============================================'

  call atm_init(n_levels)

  call random_seed(size=seed_size)
  allocate(seed(seed_size))

  ! ---- Summary file ----
  open(99, file=trim(outdir)//'/run_summary.txt', status='replace')
  write(99,'(a)') '# ens  final_n_orgs  run_ok'

  ! ---- Ensemble loop ----
  do iens = 1, n_ensemble

    ! Reproducible per-member seed
    seed = seed_base + iens * 1000 + [(istep, istep = 1, seed_size)]
    call random_seed(put=seed)

    call bio_init_run(n_init_orgs, m_init, v_conv, halflife, &
                      dt_hrs, b_total_kg, b_factor, growth_rate_day)

    ! Open output files for this member
    write(fname,'(a,"/ensemble_",i3.3,"_state.dat")') trim(outdir), iens
    open(20, file=trim(fname), status='replace')
    write(20,'(a)') &
      '# z[m]  mass[kg]  radius[m]  G[-]  rho_org[kg/m3]  age[days]  skin_width[m]'

    write(fname,'(a,"/ensemble_",i3.3,"_pop.dat")') trim(outdir), iens
    open(21, file=trim(fname), status='replace')
    write(21,'(a)') '# day  n_orgs'

    ! ---- Time loop ----
    do istep = 1, n_steps

      call bio_step(n_born, n_died)

      ! Population time series (every 24 h)
      if (mod(istep, n_steps_per_day) == 0) then
        write(21,'(f10.2,i8)') &
          real(istep, 8) * dt_hrs / 24.0d0, n_orgs
      end if

      ! Annual progress to stdout
      if (mod(istep, n_steps_last_yr) == 0) then
        write(*,'(a,i3,a,f6.1,a,i7)') &
          ' [ens', iens, '] year', &
          real(istep, 8) * dt_hrs / (24.0d0 * 365.25d0), &
          '  n_orgs=', n_orgs
      end if

      ! Organism state snapshots (last simulated year only)
      if (istep > n_steps - n_steps_last_yr) then
        call bio_write_state(20)
      end if

    end do

    close(20)
    close(21)

    write(*,'(a,i3,a,i6,a)') &
      ' Member', iens, ' done.  Final population:', n_orgs, ' organisms.'
    write(99,'(2i8,2x,a)') iens, n_orgs, 'ok'

    call bio_cleanup()
  end do

  close(99)
  call atm_cleanup()
  deallocate(seed)

  write(*,'(a)') ' Run complete.'

end program bio_driver
