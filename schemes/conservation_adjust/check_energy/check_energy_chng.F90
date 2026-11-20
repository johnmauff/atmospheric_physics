module check_energy_chng
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public  :: check_energy_chng_init
  public  :: check_energy_chng_timestep_init
  public  :: check_energy_chng_run

  ! Private module options.
  logical :: print_energy_errors = .false.    ! Turn on verbose output identifying columns that fail
                                              ! energy/water checks?

contains

!> \section arg_table_check_energy_chng_init Argument Table
!! \htmlinclude arg_table_check_energy_chng_init.html
  subroutine check_energy_chng_init(print_energy_errors_in)
    ! Input arguments
    logical,            intent(in)    :: print_energy_errors_in

    print_energy_errors = print_energy_errors_in
  end subroutine check_energy_chng_init

  ! Compute initial values of energy and water integrals,
  ! and zero out cumulative boundary tendencies.
!> \section arg_table_check_energy_chng_timestep_init Argument Table
!! \htmlinclude arg_table_check_energy_chng_timestep_init.html
  subroutine check_energy_chng_timestep_init( &
       ncol, pver, pcnst, &
       is_first_timestep, &
       q, pdel, &
       u, v, T, &
       pintdry, phis, zm, &
       cp_phys, &              ! cpairv generally, cpair fixed size for subcolumns code
       cp_or_cv_dycore, &
       te_ini_phys, te_ini_dyn, &
       tw_ini, &
       te_cur_phys, te_cur_dyn, &
       tw_cur, &
       tend_te_tnd, tend_tw_tnd, &
       temp_ini, z_ini, &
       count, &
       teout, &
       energy_formula_physics, energy_formula_dycore, &
       wv_idx, liq_idx, ice_idx, enthalpy_reference_state, &
       rga, latvap, latice, &
       errmsg, errflg)

    use check_energy_hydrostatic_energy, only: hydrostatic_energy_column, enthalpy_reference_code
    use check_energy_hydrostatic_energy, only: ENERGY_FORMULA_DYCORE_SE, ENERGY_FORMULA_DYCORE_MPAS

    ! Input arguments
    integer,            intent(in)    :: ncol           ! number of atmospheric columns
    integer,            intent(in)    :: pver           ! number of vertical layers
    integer,            intent(in)    :: pcnst          ! number of ccpp constituents
    logical,            intent(in)    :: is_first_timestep ! is first step of initial run?
    real(kind_phys),    intent(in)    :: q(:,:,:)       ! constituent mass mixing ratios [kg kg-1]
    real(kind_phys),    intent(in)    :: pdel(:,:)      ! layer thickness [Pa]
    real(kind_phys),    intent(in)    :: u(:,:)         ! zonal wind [m s-1]
    real(kind_phys),    intent(in)    :: v(:,:)         ! meridional wind [m s-1]
    real(kind_phys),    intent(in)    :: T(:,:)         ! temperature [K]
    real(kind_phys),    intent(in)    :: pintdry(:,:)   ! interface pressure dry [Pa]
    real(kind_phys),    intent(in)    :: phis(:)        ! surface geopotential [m2 s-2]
    real(kind_phys),    intent(in)    :: zm(:,:)        ! geopotential height at layer midpoints [m]
    real(kind_phys),    intent(in)    :: cp_phys(:,:)   ! enthalpy (cpairv generally) [J kg-1 K-1]
    real(kind_phys),    intent(in)    :: cp_or_cv_dycore(:,:)  ! enthalpy or heat capacity, dycore dependent [J K-1 kg-1]
    integer,            intent(in)    :: energy_formula_physics! total energy formulation physics
    integer,            intent(in)    :: energy_formula_dycore ! total energy formulation dycore
    integer,            intent(in)    :: wv_idx         ! constituent index of water vapor
    integer,            intent(in)    :: liq_idx(:)     ! constituent indices of thermodynamically active liquid water species
    integer,            intent(in)    :: ice_idx(:)     ! constituent indices of thermodynamically active ice species
    character(len=*),   intent(in)    :: enthalpy_reference_state ! 'ice', 'liq' or 'wv'
    real(kind_phys),    intent(in)    :: rga            ! 1/gravit [s2 m-1]
    real(kind_phys),    intent(in)    :: latvap         ! latent heat of vaporization [J kg-1]
    real(kind_phys),    intent(in)    :: latice         ! latent heat of fusion [J kg-1]

    ! Output arguments
    real(kind_phys),    intent(out)   :: temp_ini(:,:)  ! initial temperature [K]
    real(kind_phys),    intent(out)   :: z_ini(:,:)     ! initial geopotential height [m]
    integer,            intent(out)   :: count          ! count of values with significant energy or water imbalances [1]
    real(kind_phys),    intent(out)   :: teout(:)       ! total energy for global fixer in next timestep [J m-2]
    real(kind_phys),    intent(out)   :: tend_te_tnd(:) ! total energy tendency [J m-2 s-1]
    real(kind_phys),    intent(out)   :: tend_tw_tnd(:) ! total water tendency [kg m-2 s-1]

    ! Input/Output arguments
    real(kind_phys),    intent(inout) :: te_ini_phys(:) ! physics formula: initial total energy [J m-2]
    real(kind_phys),    intent(inout) :: te_ini_dyn (:) ! dycore  formula: initial total energy [J m-2]
    real(kind_phys),    intent(inout) :: tw_ini     (:) ! initial total water [kg m-2]
    real(kind_phys),    intent(inout) :: te_cur_phys(:) ! physics formula: current total energy [J m-2]
    real(kind_phys),    intent(inout) :: te_cur_dyn (:) ! dycore  formula: current total energy [J m-2]
    real(kind_phys),    intent(inout) :: tw_cur     (:) ! current total water [kg m-2]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local variables
    integer :: i, k
    integer :: enthalpy_reference   ! enthalpy reference state code
    integer :: colerr               ! per-column error code (validated on the host)

    errmsg = ''
    errflg = 0

    ! The energy formula and enthalpy reference state are the same for every
    ! column: validate them here (hydrostatic_energy_column's error cases),
    ! so the column kernels cannot fail.
    enthalpy_reference = enthalpy_reference_code(enthalpy_reference_state)
    call check_energy_inputs(energy_formula_physics, enthalpy_reference, &
                             enthalpy_reference_state, errmsg, errflg)
    if (errflg /= 0) return

    ! Zero out current energy unbalances count
    count = 0

    ! OpenACC: one thread per column (levels and species summed in order
    ! inside hydrostatic_energy_column, as before).
    !$acc data copyin(q, pdel, cp_phys, cp_or_cv_dycore, u, v, T, pintdry, phis, zm, &
    !$acc             liq_idx, ice_idx) &
    !$acc      copy(te_ini_phys, te_ini_dyn, tw_ini, te_cur_phys, te_cur_dyn, tw_cur, teout) &
    !$acc      copyout(temp_ini, z_ini, tend_te_tnd, tend_tw_tnd)

    !$acc parallel loop gang vector private(k, colerr)
    do i = 1, ncol
      !------------------------------------------------
      ! Physics total energy.
      !------------------------------------------------
      call hydrostatic_energy_column(i, q, pdel, cp_phys, u, v, T, &
           energy_formula_physics, pintdry, phis, wv_idx, liq_idx, ice_idx, &
           enthalpy_reference, rga, latvap, latice, colerr, &
           te=te_ini_phys, H2O=tw_ini)       ! vertically integrated total energy and water

      ! Save initial state temperature and geopotential height for use in run phase
      !$acc loop seq
      do k = 1, size(temp_ini, 2)
        temp_ini(i, k) = T (i, k)
        z_ini   (i, k) = zm(i, k)
      end do

      !------------------------------------------------
      ! Dynamical core total energy.
      !------------------------------------------------
      if (energy_formula_dycore == ENERGY_FORMULA_DYCORE_SE) then
        ! SE dycore specific hydrostatic energy (enthalpy)
        call hydrostatic_energy_column(i, q, pdel, cp_or_cv_dycore, u, v, T, &
             energy_formula_dycore, pintdry, phis, wv_idx, liq_idx, ice_idx, &
             enthalpy_reference, rga, latvap, latice, colerr, &
             te=te_ini_dyn)                  ! WRITE OPERATION - vertically integrated total energy
      else if (energy_formula_dycore == ENERGY_FORMULA_DYCORE_MPAS) then
        ! MPAS dycore: compute cv if vertical coordinate is height: cv = cp - R (internal energy)
        call hydrostatic_energy_column(i, q, pdel, cp_or_cv_dycore, u, v, T, &
             energy_formula_dycore, pintdry, phis, wv_idx, liq_idx, ice_idx, &
             enthalpy_reference, rga, latvap, latice, colerr, &
             z_mid=z_ini, &                  ! unique for MPAS
             te=te_ini_dyn)                  ! WRITE OPERATION - vertically integrated total energy
      else
        ! FV dycore: dycore energy is the same as physics
        te_ini_dyn(i) = te_ini_phys(i)
      endif

      ! Set current state to be the same as initial
      te_cur_phys(i) = te_ini_phys(i)
      te_cur_dyn (i) = te_ini_dyn (i)
      tw_cur     (i) = tw_ini     (i)

      ! Zero out cumulative boundary fluxes
      tend_te_tnd(i) = 0._kind_phys
      tend_tw_tnd(i) = 0._kind_phys

      ! If first timestep, initialize value of teout
      if (is_first_timestep) then
        teout(i) = te_ini_dyn(i)
      endif
    end do

    !$acc end data

  end subroutine check_energy_chng_timestep_init


  ! Check that energy and water change matches the boundary fluxes.
!> \section arg_table_check_energy_chng_run Argument Table
!! \htmlinclude arg_table_check_energy_chng_run.html
  subroutine check_energy_chng_run( &
       ncol, pver, pcnst, &
       iulog, &
       q, pdel, &
       u, v, T, &
       pintdry, phis, zm, &
       cp_phys, &              ! cpairv generally, cpair fixed size for subcolumns code
       cp_or_cv_dycore, &
       scaling_dycore, &       ! From check_energy_scaling to work around subcolumns code
       te_cur_phys, te_cur_dyn, &
       tw_cur, &
       tend_te_tnd, tend_tw_tnd, &
       temp_ini, z_ini, &
       count, ztodt, &
       latice, latvap, &
       energy_formula_physics, energy_formula_dycore, &
       name, flx_vap, flx_cnd, flx_ice, flx_sen, &
       wv_idx, liq_idx, ice_idx, enthalpy_reference_state, rga, &
       errmsg, errflg)

    use check_energy_hydrostatic_energy, only: hydrostatic_energy_column, enthalpy_reference_code
    use check_energy_hydrostatic_energy, only: ENERGY_FORMULA_DYCORE_SE, ENERGY_FORMULA_DYCORE_MPAS

    ! Input arguments
    integer,            intent(in)    :: ncol           ! number of atmospheric columns
    integer,            intent(in)    :: pver           ! number of vertical layers
    integer,            intent(in)    :: pcnst          ! number of ccpp constituents
    integer,            intent(in)    :: iulog          ! log output unit
    real(kind_phys),    intent(in)    :: q(:,:,:)       ! constituent mass mixing ratios [kg kg-1]
    real(kind_phys),    intent(in)    :: pdel(:,:)      ! layer thickness [Pa]
    real(kind_phys),    intent(in)    :: u(:,:)         ! zonal wind [m s-1]
    real(kind_phys),    intent(in)    :: v(:,:)         ! meridional wind [m s-1]
    real(kind_phys),    intent(in)    :: T(:,:)         ! temperature [K]
    real(kind_phys),    intent(in)    :: pintdry(:,:)   ! interface pressure dry [Pa]
    real(kind_phys),    intent(in)    :: phis(:)        ! surface geopotential [m2 s-2]
    real(kind_phys),    intent(in)    :: zm(:,:)        ! geopotential height at layer midpoints [m]
    real(kind_phys),    intent(in)    :: temp_ini(:,:)  ! initial temperature [K]
    real(kind_phys),    intent(in)    :: z_ini(:,:)     ! initial geopotential height [m]
    real(kind_phys),    intent(in)    :: cp_phys(:,:)   ! enthalpy (cpairv generally) [J kg-1 K-1]
    real(kind_phys),    intent(in)    :: cp_or_cv_dycore(:,:)  ! enthalpy or heat capacity, dycore dependent [J K-1 kg-1]
    real(kind_phys),    intent(in)    :: scaling_dycore(:,:)   ! scaling for conversion of temperature increment [1]
    real(kind_phys),    intent(in)    :: ztodt          ! physics timestep [s]
    real(kind_phys),    intent(in)    :: latice         ! constant, latent heat of fusion of water at 0 C [J kg-1]
    real(kind_phys),    intent(in)    :: latvap         ! constant, latent heat of vaporization of water at 0 C [J kg-1]
    integer,            intent(in)    :: energy_formula_physics! total energy formulation physics
    integer,            intent(in)    :: energy_formula_dycore ! total energy formulation dycore

    ! Input from CCPP-scheme being checked:
    ! parameterization name; surface fluxes of (1) vapor, (2) liquid+ice, (3) ice, (4) sensible heat
    ! to pass in the values to be checked, call check_energy_zero_input_fluxes to reset these values
    ! before a parameterization that is checked, then update these values as-needed
    ! (can be all zero; in fact, most parameterizations calling _chng call with zero arguments)
    !
    ! Original comment from BAB:
    ! Note that the precip and ice fluxes are in precip units (m/s).
    ! I would prefer to have kg/m2/s.
    ! I would also prefer liquid (not total) and ice fluxes
    character(len=*),   intent(in)    :: name           ! parameterization name for fluxes
    real(kind_phys),    intent(in)    :: flx_vap(:)     ! boundary flux of vapor [kg m-2 s-1]
    real(kind_phys),    intent(in)    :: flx_cnd(:)     ! boundary flux of liquid+ice (precip?) [m s-1]
    real(kind_phys),    intent(in)    :: flx_ice(:)     ! boundary flux of ice [m s-1]
    real(kind_phys),    intent(in)    :: flx_sen(:)     ! boundary flux of sensible heat [W m-2]
    integer,            intent(in)    :: wv_idx         ! constituent index of water vapor
    integer,            intent(in)    :: liq_idx(:)     ! constituent indices of thermodynamically active liquid water species
    integer,            intent(in)    :: ice_idx(:)     ! constituent indices of thermodynamically active ice species
    character(len=*),   intent(in)    :: enthalpy_reference_state ! 'ice', 'liq' or 'wv'
    real(kind_phys),    intent(in)    :: rga            ! 1/gravit [s2 m-1]

    ! Input/Output arguments
    real(kind_phys),    intent(inout) :: te_cur_phys(:) ! physics formula: current total energy [J m-2]
    real(kind_phys),    intent(inout) :: te_cur_dyn (:) ! dycore  formula: current total energy [J m-2]
    real(kind_phys),    intent(inout) :: tw_cur     (:) ! current total water [kg m-2]
    integer,            intent(inout) :: count          ! count of columns with significant energy or water imbalances [1]
    real(kind_phys),    intent(inout) :: tend_te_tnd(:) ! total energy tendency [J m-2 s-1]
    real(kind_phys),    intent(inout) :: tend_tw_tnd(:) ! total water tendency [kg m-2 s-1]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local variables
    real(kind_phys) :: te_xpd(ncol)                     ! expected value (f0 + dt*boundary_flux)
    real(kind_phys) :: te_dif(ncol)                     ! energy of input state - original energy
    real(kind_phys) :: te_tnd(ncol)                     ! tendency from last process
    real(kind_phys) :: te_rer(ncol)                     ! relative error in energy column

    real(kind_phys) :: tw_xpd(ncol)                     ! expected value (w0 + dt*boundary_flux)
    real(kind_phys) :: tw_dif(ncol)                     ! tw_inp - original water
    real(kind_phys) :: tw_tnd(ncol)                     ! tendency from last process
    real(kind_phys) :: tw_rer(ncol)                     ! relative error in water column

    real(kind_phys) :: te(ncol)                         ! vertical integral of total energy
    real(kind_phys) :: tw(ncol)                         ! vertical integral of total water
    real(kind_phys) :: temp(ncol,pver)                  ! temperature

    ! Per-column scalars (stored once).
    real(kind_phys) :: te_tnd_i, tw_tnd_i               ! boundary-flux tendencies
    real(kind_phys) :: te_xpd_i, tw_xpd_i               ! expected values
    integer :: k
    integer :: colerr                                   ! per-column error code (validated on the host)

    integer :: i
    integer :: enthalpy_reference                       ! enthalpy reference state code

    errmsg = ''
    errflg = 0

    ! The energy formula and enthalpy reference state are the same for every
    ! column: validate them here, so the column kernel cannot fail.
    enthalpy_reference = enthalpy_reference_code(enthalpy_reference_state)
    call check_energy_inputs(energy_formula_physics, enthalpy_reference, &
                             enthalpy_reference_state, errmsg, errflg)
    if (errflg /= 0) return

    ! OpenACC: one thread per column. The conservation-error report stays on
    ! the host: with print_energy_errors the per-column budget values come
    ! back after the kernel and the original loop prints them.
    ! zm is not used by the run phase; it is listed so the cap treats it as
    ! device data instead of copying it to the host before every call.
    !$acc data copyin(q, pdel, u, v, T, pintdry, phis, zm, cp_phys, cp_or_cv_dycore, &
    !$acc             scaling_dycore, temp_ini, z_ini, flx_vap, flx_cnd, flx_ice, flx_sen, &
    !$acc             liq_idx, ice_idx) &
    !$acc      copy(te_cur_phys, te_cur_dyn, tw_cur, tend_te_tnd, tend_tw_tnd) &
    !$acc      create(te, tw, temp, te_xpd, te_dif, te_tnd, te_rer, tw_xpd, tw_dif, tw_tnd, tw_rer)

    !$acc parallel loop gang vector &
    !$acc   private(k, colerr, te_tnd_i, tw_tnd_i, te_xpd_i, tw_xpd_i)
    do i = 1, ncol
      !------------------------------------------------
      ! Physics total energy.
      !------------------------------------------------
      call hydrostatic_energy_column(i, q, pdel, cp_phys, u, v, T, &
           energy_formula_physics, pintdry, phis, wv_idx, liq_idx, ice_idx, &
           enthalpy_reference, rga, latvap, latice, colerr, &
           te=te, H2O=tw)                    ! vertically integrated total energy and water

      ! compute expected values and tendencies
      ! change in static energy and total water
      te_dif(i) = te(i) - te_cur_phys(i)
      tw_dif(i) = tw(i) - tw_cur     (i)

      ! expected tendencies from boundary fluxes for last process
      te_tnd_i = flx_vap(i)*(latvap+latice) - (flx_cnd(i) - flx_ice(i))*1000._kind_phys*latice + flx_sen(i)
      tw_tnd_i = flx_vap(i) - flx_cnd(i) *1000._kind_phys
      te_tnd(i) = te_tnd_i
      tw_tnd(i) = tw_tnd_i

      ! cummulative tendencies from boundary fluxes
      tend_te_tnd(i) = tend_te_tnd(i) + te_tnd_i
      tend_tw_tnd(i) = tend_tw_tnd(i) + tw_tnd_i

      ! expected new values from previous state plus boundary fluxes
      te_xpd_i = te_cur_phys(i) + te_tnd_i*ztodt
      tw_xpd_i = tw_cur     (i) + tw_tnd_i*ztodt
      te_xpd(i) = te_xpd_i
      tw_xpd(i) = tw_xpd_i

      ! relative error, expected value - input state / previous state
      te_rer(i) = (te_xpd_i - te(i)) / te_cur_phys(i)

      ! relative error for total water (allow for dry atmosphere)
      tw_rer(i) = 0._kind_phys
      if (tw_cur(i) > 0._kind_phys) then
        tw_rer(i) = (tw_xpd_i - tw(i)) / tw_cur(i)
      end if

      ! WRITE OPERATION - copy new value to state, including total water.
      ! the total water operations are consistent regardless of energy formula,
      ! so it only has to be written once.
      te_cur_phys(i) = te(i)
      tw_cur(i)      = tw(i)

      !------------------------------------------------
      ! Dynamical core total energy.
      !------------------------------------------------
      if (energy_formula_dycore == ENERGY_FORMULA_DYCORE_SE .or. &
          energy_formula_dycore == ENERGY_FORMULA_DYCORE_MPAS) then
        ! enthalpy scaling for energy consistency
        ! REMOVECAM: note this scaling is different with subcols off/on which is why it was put into
        ! separate scheme (hplin, 9/5/24)
        !$acc loop seq
        do k = 1, size(temp, 2)
          temp(i, k) = temp_ini(i, k)+scaling_dycore(i, k)*(T(i, k)-temp_ini(i, k))
        end do
        if (energy_formula_dycore == ENERGY_FORMULA_DYCORE_SE) then
          ! SE dycore specific hydrostatic energy
          call hydrostatic_energy_column(i, q, pdel, cp_or_cv_dycore, u, v, temp, &
               energy_formula_dycore, pintdry, phis, wv_idx, liq_idx, ice_idx, &
               enthalpy_reference, rga, latvap, latice, colerr, &
               te=te_cur_dyn)                ! WRITE OPERATION - vertically integrated total energy
        else
          ! MPAS dycore: compute cv if vertical coordinate is height: cv = cp - R
          call hydrostatic_energy_column(i, q, pdel, cp_or_cv_dycore, u, v, temp, &
               energy_formula_dycore, pintdry, phis, wv_idx, liq_idx, ice_idx, &
               enthalpy_reference, rga, latvap, latice, colerr, &
               z_mid=z_ini, &                ! unique for MPAS
               te=te_cur_dyn)                ! WRITE OPERATION - vertically integrated total energy
        end if
      else
        ! FV dycore
        te_cur_dyn(i) = te(i)
      end if
    end do

    ! error checking (on the host, columns in order)
    if (print_energy_errors) then
      !$acc update self(te, te_xpd, te_dif, te_tnd, te_rer, tw, tw_xpd, tw_dif, tw_tnd, tw_rer, &
      !$acc             tend_te_tnd, tend_tw_tnd)
      if (any(abs(te_rer(1:ncol)) > 1.E-14_kind_phys .or. abs(tw_rer(1:ncol)) > 1.E-10_kind_phys)) then
        do i = 1, ncol
          ! the relative error threshold for the water budget has been reduced to 1.e-10
          ! to avoid messages generated by QNEG3 calls
          ! PJR- change to identify if error in energy or water
          if (abs(te_rer(i)) > 1.E-14_kind_phys ) then
            count = count + 1
            write(iulog,*) "significant energy conservation error after ", name,        &
                  " count", count, "col", i
            write(iulog,*) te(i),te_xpd(i),te_dif(i),tend_te_tnd(i)*ztodt,  &
                  te_tnd(i)*ztodt,te_rer(i)
          endif
          if ( abs(tw_rer(i)) > 1.E-10_kind_phys) then
            count = count + 1
            write(iulog,*) "significant water conservation error after ", name,        &
                  " count", count, "col", i
            write(iulog,*) tw(i),tw_xpd(i),tw_dif(i),tend_tw_tnd(i)*ztodt,  &
                  tw_tnd(i)*ztodt,tw_rer(i)
          end if
        end do
      end if
    end if

    !$acc end data

  end subroutine check_energy_chng_run

  ! Check the energy formula and enthalpy reference state used for every
  ! column (hydrostatic_energy_column's error cases; the endrun messages of
  ! cam_thermo's get_hydrostatic_energy).
  subroutine check_energy_inputs(vcoord, enthalpy_reference, enthalpy_reference_state, &
                                 errmsg, errflg)
    use check_energy_hydrostatic_energy, only: ENERGY_FORMULA_DYCORE_FV, &
         ENERGY_FORMULA_DYCORE_SE, ENERGY_FORMULA_DYCORE_MPAS
    integer,            intent(in)  :: vcoord
    integer,            intent(in)  :: enthalpy_reference
    character(len=*),   intent(in)  :: enthalpy_reference_state
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0
    if (vcoord /= ENERGY_FORMULA_DYCORE_FV .and. vcoord /= ENERGY_FORMULA_DYCORE_SE .and. &
        vcoord /= ENERGY_FORMULA_DYCORE_MPAS) then
      errflg = 1
      write(errmsg, '(a,i0)') 'check_energy_chng: energy formula not supported: ', vcoord
    else if (enthalpy_reference == 0) then
      errflg = 1
      errmsg = 'check_energy_chng: enthalpy reference state not supported: ' // &
               trim(enthalpy_reference_state)
    end if
  end subroutine check_energy_inputs

end module check_energy_chng
