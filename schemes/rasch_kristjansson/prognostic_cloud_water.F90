! Prognostic cloud water data and methods (cldwat)
! Original authors as marked in subroutines
! CCPP-ized: Haipeng Lin, January 2025
module prognostic_cloud_water
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: prognostic_cloud_water_init
  public :: prognostic_cloud_water_run

  ! tuning parameters used by prognostic cloud water (RK stratiform)
  real(kind_phys), public :: icritc ! threshold for autoconversion of cold ice [kg kg-1]
  ! REMOVECAM: icritc does not have to be public after CAM is retired (namelist option; CCPP framework can provide it)
  real(kind_phys) :: icritw         ! threshold for autoconversion of warm ice [kg kg-1]
  real(kind_phys) :: conke          ! tunable constant for evaporation of precipitation [kg-0.5 m s-0.5]
  real(kind_phys) :: r3lcrit        ! critical radius where liq conversion begins [m]

  logical         :: do_psrhmin     ! set rh in stratosphere poleward of 50 degrees [flag]
  real(kind_phys) :: psrhmin        ! rh set in stratosphere poleward of 50 degrees [%]

  ! module private variables
  real(kind_phys) :: rhonot         ! air density at surface [g cm-3]
  real(kind_phys) :: rhos           ! assumed snow density [g cm-3]
  real(kind_phys) :: rhow           ! water density [g cm-3]
  real(kind_phys) :: esi            ! Collection efficiency for ice by snow [1]
  real(kind_phys) :: esw            ! Collection efficiency for water by snow [1]
  real(kind_phys) :: t0             ! Approx. freezing temperature [K]
  real(kind_phys) :: cldmin         ! Assumed minimum cloud amount [1]
  real(kind_phys) :: small          ! Small number compared to unity [1]
  real(kind_phys) :: c              ! Constant for graupel-like snow [cm^(1-d)/s]

  real(kind_phys) :: d              ! Constant for graupel-like snow [1]
  real(kind_phys) :: thrpd          ! 3+d
  real(kind_phys) :: gam3pd         ! gamma(3+d)
  real(kind_phys) :: gam4pd         ! gamma(4+d)

  real(kind_phys) :: nos            ! Particles snow concentration [cm^-4]
  real(kind_phys) :: prhonos        ! pi * rhos * nos

  ! cloud microphysics constants
  real(kind_phys) :: mcon01
  real(kind_phys) :: mcon02
  real(kind_phys) :: mcon03
  real(kind_phys) :: mcon04
  real(kind_phys) :: mcon05
  real(kind_phys) :: mcon06
  real(kind_phys) :: mcon07
  real(kind_phys) :: mcon08

  ! Parameters for findmcnew
  real(kind_phys) :: capnw          ! Warm continental cloud particle number concentration [cm-3]
  real(kind_phys) :: capnc          ! Cold/oceanic cloud particle number concentration [cm-3]
  real(kind_phys) :: capnsi         ! Sea ice cloud particle number concentration [cm-3]
  real(kind_phys) :: kconst         ! Terminal velocity constant (Stokes regime) [1]
  real(kind_phys) :: effc           ! Autoconv collection efficiency [1]
  real(kind_phys) :: alpha          ! Ratio of 3rd moment radius to 2nd [1]
  real(kind_phys) :: capc           ! Autoconversion constant [1]
  real(kind_phys) :: convfw         ! Fall velocity calculation constant [1]
  real(kind_phys) :: cracw          ! Rain accreting water constant [1]
  real(kind_phys) :: critpr         ! Critical precip rate for collection efficiency [mm day-1]
  real(kind_phys) :: ciautb         ! Ice autoconversion coefficient [s-1]

  ! module variables used on the device (set in prognostic_cloud_water_init)
  !$acc declare create(icritc, icritw, conke, r3lcrit, do_psrhmin, psrhmin, &
  !$acc&               rhow, esi, esw, t0, cldmin, mcon05, mcon07, mcon08, &
  !$acc&               capnw, capnc, capnsi, capc, critpr, convfw, cracw, ciautb)

  ! errors found on the device in prognostic_cloud_water_run, reported on the host
  integer, parameter :: ERR_NONE               = 0
  integer, parameter :: ERR_FINDSP_CONVERGENCE = 1  ! findsp status 2
  integer, parameter :: ERR_FINDSP_ENTHALPY    = 2  ! findsp status 8
  integer, parameter :: ERR_EMPTY_RH_CLOUD     = 3
  integer, parameter :: ERR_NEGATIVE_RCWN      = 4
  integer, parameter :: ERR_NEGATIVE_RLIQ      = 5
  integer, parameter :: ERR_NEGATIVE_RICE      = 6

contains

  ! Initialize prognostic cloud water module and calculate microphysical constants
!> \section arg_table_prognostic_cloud_water_init Argument Table
!! \htmlinclude arg_table_prognostic_cloud_water_init.html
  subroutine prognostic_cloud_water_init( &
    amIRoot, iulog, &
    tmelt, rhodair, pi, &
    icritc_in, icritw_in, &
    conke_in, r3lcrit_in, &
    do_psrhmin_in, psrhmin_in, &
    one_mom_clouds, &
    errmsg, errflg)

    ! Input arguments
    logical,            intent(in)    :: amIRoot       ! are we on the MPI root task?
    integer,            intent(in)    :: iulog         ! log output unit

    real(kind_phys),    intent(in)    :: tmelt         ! freezing_point_of_water [K]
    real(kind_phys),    intent(in)    :: rhodair       ! density_of_dry_air_at_stp [kg m-3]
    real(kind_phys),    intent(in)    :: pi

    real(kind_phys),    intent(in)    :: icritc_in
    real(kind_phys),    intent(in)    :: icritw_in
    real(kind_phys),    intent(in)    :: conke_in
    real(kind_phys),    intent(in)    :: r3lcrit_in
    logical,            intent(in)    :: do_psrhmin_in
    real(kind_phys),    intent(in)    :: psrhmin_in

    ! Output arguments
    logical,            intent(out)   :: one_mom_clouds ! is_one_moment_cloud_microphysics [flag]
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    errmsg = ''
    errflg = 0

    ! RK prognostic cloud water is a one-moment cloud microphysics: declare
    ! it so downstream consumers (convective_cloud_water emissivity averaging,
    ! cloud_water_paths_diagnostics one-moment/two-moment history split) can
    ! branch without querying which microphysics scheme is active.
    one_mom_clouds = .true.

    ! First populate tuning parameters in-module
    icritc = icritc_in
    icritw = icritw_in
    conke  = conke_in
    r3lcrit = r3lcrit_in
    do_psrhmin = do_psrhmin_in
    psrhmin = psrhmin_in

    if(amIRoot) then
      write(iulog,*) 'tuning parameters prognostic_cloud_water_init: icritw ',icritw,' icritc ',icritc,' conke ',conke,' r3lcrit ',r3lcrit
      write(iulog,*) 'prognostic_cloud_water_init: do_psrhmin = ', do_psrhmin
    endif

    !--------------------------------------------------
    ! Initialize constants used for prognostic condensate (inimc)
    ! Original author: P. Rasch, April 1997
    !--------------------------------------------------

    rhonot = rhodair/1000.0_kind_phys     ! convert from kg m-3 to g cm-3

    ! assumed densities of snow, water [g cm-3]
    rhos   = 0.1_kind_phys
    rhow   = 1._kind_phys

    esi    = 1._kind_phys
    esw    = 0.1_kind_phys
    t0     = tmelt
    cldmin = 0.02_kind_phys
    small  = 1.e-22_kind_phys
    c      = 152.93_kind_phys
    d      = 0.25_kind_phys               ! constant for graupel like snow
                                          ! values other than 0.25 are not supported.
    if(d == 0.25_kind_phys) then
      gam3pd = 2.549256966718531_kind_phys
      gam4pd = 8.285085141835282_kind_phys

      ! on cray machines this can be calculated using gamma function:
      ! call gamma(3._kind_phys+d, signgam, gam3pd)
      ! gam3pd = sign(exp(gam3pd),signgam)
      ! call gamma(4._kind_phys+d, signgam, gam4pd)
      ! gam4pd = sign(exp(gam4pd),signgam)
      ! write(iulog,*) ' d, gamma(3+d), gamma(4+d) =', gam3pd, gam4pd
    else
      errflg = 1
      errmsg = 'prognostic_cloud_water_init: cannot use d /= 0.25'
    endif

    nos    = 3.e-2_kind_phys

    ! note: pi was originally calculated here as 4._kind_phys * atan(1.0_kind_phys)
    ! changed in ccpp-ization to use physconst value.
    prhonos = pi*rhos*nos
    thrpd  = 3._kind_phys + d

    mcon01 = pi*nos*c*gam3pd/4._kind_phys
    mcon02 = 1._kind_phys/(c*gam4pd*sqrt(rhonot)/(6*prhonos**(d/4._kind_phys)))
    mcon03 = -(0.5_kind_phys+d/4._kind_phys)
    mcon04 = 4._kind_phys/(4._kind_phys+d)
    mcon05 = (3+d)/(4+d)
    mcon06 = (3+d)/4._kind_phys
    mcon07 = mcon01*sqrt(rhonot)*mcon02**mcon05/prhonos**mcon06
    mcon08 = -0.5_kind_phys/(4._kind_phys+d)

    ! Initialize parameters used by findmcnew
    ! Cloud particle densities [cm-3] for ...
    capnw  = 400._kind_phys  ! ...warm continental
    capnc  = 150._kind_phys  ! ...cold and oceanic
    capnsi = 75._kind_phys   ! ...sea ice

    kconst = 1.18e6_kind_phys

    ! Autoconv collection efficiencies.
    ! Tripoli and Cotton (default) = 0.55
    ! Boucher 96 = 1.
    ! Baker 93 = 0.55*0.05
    ! Turn off = 0
    effc   = 0.55_kind_phys

    alpha  = 1.1_kind_phys ** 4
    ! constant for autoconversion
    capc = pi**(-.333_kind_phys)*kconst*effc *(0.75_kind_phys)**(1.333_kind_phys)*alpha

    ! critical precip rate at which we assume the collector drops can change the
    ! drop size enough to enhance the auto-conversion process (mm/day)
    critpr = 0.5_kind_phys
    convfw = 1.94_kind_phys*2.13_kind_phys*sqrt(rhow*1000._kind_phys*9.81_kind_phys*2.7e-4_kind_phys)

    ! liquid microphysics - rain accreting water constant
    ! Tripoli and Cotton = default
    ! Beheng = 6.
    cracw = .884_kind_phys*sqrt(9.81_kind_phys/(rhow*1000._kind_phys*2.7e-4_kind_phys)) ! tripoli and cotton

    ! ice microphysics autoconversion coefficient
    ciautb = 5.e-4_kind_phys

    !$acc update device(icritc, icritw, conke, r3lcrit, do_psrhmin, psrhmin, &
    !$acc&              rhow, esi, esw, t0, cldmin, mcon05, mcon07, mcon08, &
    !$acc&              capnw, capnc, capnsi, capc, critpr, convfw, cracw, ciautb)

    if(amIRoot) then
      write(iulog,*) 'tuning parameters prognostic_cloud_water_init: capnw ',capnw,' capnc ',capnc,' capnsi ',capnsi,' kconst ',kconst
      write(iulog,*) 'tuning parameters prognostic_cloud_water_init: effc ',effc,' alpha ',alpha,' capc ',capc
      write(iulog,*) 'tuning parameters prognostic_cloud_water_init: critpr ',critpr,' convfw ',convfw,' cracw ',cracw,' ciautb ',ciautb
    endif

  end subroutine prognostic_cloud_water_init

  ! Calculate prognostic condensate.
  ! Cloud water parameterization, returns tendencies to water vapor, temperature,
  ! and cloud water variables.
  !
  ! Rasch, P. J., and J. E. Kristjánsson, 1998: A Comparison of the CCM3 Model Climate Using Diagnosed and Predicted Condensate Parameterizations. J. Climate, 11, 1587–1614
  ! https://doi.org/10.1175/1520-0442(1998)011<1587:ACOTCM>2.0.CO;2
  !
  ! with modifications to improve the method of determining condensation/evaporation
  ! Zhang, M., W. Lin, C. Bretherton, J. Hack, and P. J. Rasch, 2003, A modified formulation of fractional stratiform condensation rate in the NCAR Community Atmospheric Model (CAM2), J. Geophys. Res., 108(D1), 4035
  ! https://doi.org/10.1029/2002JD002523
  !
  ! Original authors: M. Zhang, W. Lin, P. Rasch and J.E. Kristjansson
  !                   B. A. Boville (latent heat of fusion)
!> \section arg_table_prognostic_cloud_water_run Argument Table
!! \htmlinclude arg_table_prognostic_cloud_water_run.html
  subroutine prognostic_cloud_water_run( &
    ncol, pver, top_lev, deltat, &
    iulog, &
    pi, gravit, rh2o, epsilo, latvap, latice, cpair, &
    dlat, &
    pmid, pdel, &
    zi, &
    troplev, &
    ttend, tn, &
    qtend, qn, &
    ltend, &
    cldice, cldliq, &
    omega, &
    cldn, &
    fice, fsnow, &
    rhdfda, rhu00, &
    landm, seaicef, snowh, &
    qme, &
    prodprec, prodsnow, &
    evapprec, evapsnow, &
    evapheat, prfzheat, meltheat, &
    precip, snowab, &
    ice2pr, liq2pr, liq2snow, &
    lsflxprc, &
    lsflxsnw, &
    pracwo, psacwo, psacio, &
    fwaut, fsaut, fracw, fsacw, fsaci, &
    errmsg, errflg)
    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: top_lev        ! vertical_layer_index_of_troposphere_cloud_physics_top [index]
    real(kind_phys),    intent(in)    :: deltat         ! timestep [s]
    integer,            intent(in)    :: iulog          ! log output unit [1]
    real(kind_phys),    intent(in)    :: pi             ! pi_constant [1]
    real(kind_phys),    intent(in)    :: gravit         ! gravitational acceleration [m s-2]
    real(kind_phys),    intent(in)    :: rh2o           ! gas_constant_of_water_vapor [J K-1 kg-1]
    real(kind_phys),    intent(in)    :: epsilo         ! ratio_of_water_vapor_to_dry_air_molecular_weights [1]
    real(kind_phys),    intent(in)    :: latvap         ! latent_heat_of_vaporization_of_water_at_0c [J kg-1]
    real(kind_phys),    intent(in)    :: latice         ! latent_heat_of_fusion_of_water_at_0c [J kg-1]
    real(kind_phys),    intent(in)    :: cpair          ! specific_heat_of_dry_air_at_constant_pressure [J K-1 kg-1]
    real(kind_phys),    intent(in)    :: dlat(:)        ! latitude_degrees_north [degrees]
    real(kind_phys),    intent(in)    :: pmid(:,:)      ! air_pressure [Pa]
    real(kind_phys),    intent(in)    :: pdel(:,:)      ! air_pressure_thickness [Pa]
    real(kind_phys),    intent(in)    :: zi(:,:)        ! geopotential_height_wrt_surface_at_interface [m]
    integer,            intent(in)    :: tropLev(:)     ! tropopause_vertical_layer_index [index]

    real(kind_phys),    intent(in)    :: ttend(:,:)     ! Temperature tendency [K s-1] -- from non-micro/macrophysics
    real(kind_phys),    intent(in)    :: tn(:,:)        ! New Temperature [K]
    real(kind_phys),    intent(in)    :: qtend(:,:)     ! Water vapor tendency [kg kg-1 s-1] -- from non-micro/macrophysics
    real(kind_phys),    intent(in)    :: qn(:,:)        ! New Water vapor mixing ratio [kg kg-1]
    real(kind_phys),    intent(in)    :: ltend(:,:)     ! Cloud liquid water tendency [kg kg-1 s-1] -- from non-micro/macrophysics
    real(kind_phys),    intent(in)    :: cldice(:,:)    ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: cldliq(:,:)    ! adv: cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]

    real(kind_phys),    intent(in)    :: omega(:,:)     ! lagrangian_tendency_of_air_pressure [Pa s-1]
    real(kind_phys),    intent(in)    :: cldn(:,:)      ! New Cloud fraction [fraction]
    real(kind_phys),    intent(in)    :: fice(:,:)      ! Ice fraction within cwat [fraction]
    real(kind_phys),    intent(in)    :: fsnow(:,:)     ! Fraction of rain that freezes to snow [fraction]

    real(kind_phys),    intent(in)    :: rhdfda(:,:)    ! dG(a)/da, rh=G(a), when rh>u00 [??]
    real(kind_phys),    intent(in)    :: rhu00(:,:)     ! Rhlim for cloud [percent]
    real(kind_phys),    intent(in)    :: landm(:)       ! smoothed_land_area_fraction [fraction]
    real(kind_phys),    intent(in)    :: seaicef(:)     ! Sea ice fraction [1 (fraction)]
    real(kind_phys),    intent(in)    :: snowh(:)       ! Snow depth over land, water equivalent [m]

    ! Output arguments
    real(kind_phys),    intent(out)   :: qme(:,:)       ! Rate of condensation-evaporation of condensate (net_condensation_rate_due_to_microphysics) [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: prodprec(:,:)  ! Conversion rate of condensate to precip (precipitation_production_due_to_microphysics) [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: prodsnow(:,:)  ! Snow production rate (ignored in RK?) [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: evapprec(:,:)  ! Falling precipitation evaporation rate (rate_of_evaporation_of_precipitation_due_to_microphysics) [kg kg-1 s-1] -- & combined to apply q(wv) tendency
    real(kind_phys),    intent(out)   :: evapsnow(:,:)  ! Falling snow evaporation rate [kg kg-1 s-1]
    real(kind_phys),    intent(out)   :: evapheat(:,:)  ! heating rate due to evaporation of precipitation [J kg-1 s-1]
    real(kind_phys),    intent(out)   :: prfzheat(:,:)  ! heating rate due to freezing of precipitation [J kg-1 s-1]
    real(kind_phys),    intent(out)   :: meltheat(:,:)  ! heating rate due to snow melt [J kg-1 s-1]
    ! note -- these are precip units for atmosphere-surface exchange.
    real(kind_phys),    intent(out)   :: precip(:)      ! Precipitation rate (lwe_stratiform_precipitation_rate_at_surface) [m s-1]
    real(kind_phys),    intent(out)   :: snowab(:)      ! Snow rate (lwe_snow_precipitation_rate_at_surface_due_to_microphysics) [m s-1]
    ! / note
    real(kind_phys),    intent(out)   :: ice2pr(:,:)    ! Conversion rate of ice to precip [kg kg-1 s-1] -- apply q(cldice) tendency
    real(kind_phys),    intent(out)   :: liq2pr(:,:)    ! Conversion rate of liquid to precip [kg kg-1 s-1] -- apply q(cldliq) tendency
    real(kind_phys),    intent(out)   :: liq2snow(:,:)  ! Conversion rate of liquid to snow [kg kg-1 s-1] -- discard??
    real(kind_phys),    intent(out)   :: lsflxprc(:,:)  ! Grid-box mean, flux large scale cloud rain+snow at interfaces (stratiform_rain_and_snow_flux_at_interface) [kg m-2 s-1]
    real(kind_phys),    intent(out)   :: lsflxsnw(:,:)  ! Grid-box mean, flux large scale cloud snow at interfaces (stratiform_snow_flux_at_interface) [kg m-2 s-1]
    real(kind_phys),    intent(out)   :: pracwo(:,:)    ! accretion of cloud water by rain [s-1]
    real(kind_phys),    intent(out)   :: psacwo(:,:)    ! accretion of cloud water by snow [s-1]
    real(kind_phys),    intent(out)   :: psacio(:,:)    ! accretion of cloud ice by snow [s-1]
    real(kind_phys),    intent(out)   :: fwaut(:,:)                      ! Relative importance of liquid autoconversion [fraction]
    real(kind_phys),    intent(out)   :: fsaut(:,:)                      ! Relative importance of ice autoconversion [fraction]
    real(kind_phys),    intent(out)   :: fracw(:,:)                      ! Relative importance of liquid collection by rain [fraction]
    real(kind_phys),    intent(out)   :: fsacw(:,:)                      ! Relative importance of liquid collection by snow [fraction]
    real(kind_phys),    intent(out)   :: fsaci(:,:)                      ! Relative importance of ice collection by snow [fraction]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local variables
    real(kind_phys),    parameter     :: rhofw = 1000._kind_phys  ! density of fresh water [kg m-3]
    integer :: i                                        ! Column index [1]
    integer :: iter                                     ! # of iterations for precipitation calculation [1]

    ! Constants and parameters
    real(kind_phys) :: mincld                           ! Minimum cloud fraction [1]
    real(kind_phys) :: cpohl                            ! Ratio of specific heat to latent heat [K-1]
    real(kind_phys) :: hlocp                            ! Ratio of latent heat to specific heat [K]
    real(kind_phys) :: dto2                             ! Half timestep [s]

    ! Errors found on the device, per column (see prognostic_cloud_water_column)
    integer         :: err_code(ncol)                   ! 0: none, else ERR_* [1]
    integer         :: err_k(ncol)                      ! level of the error [index]
    integer         :: err_l(ncol)                      ! iteration of the error (0: findsp) [index]
    real(kind_phys) :: err_val(ncol, 6)                 ! values for the error message
    integer         :: ierr                             ! column of the error reported [index]
    integer         :: key, key_min                     ! order of the errors in the original loops

    errmsg = ''
    errflg = 0

    cpohl  = cpair/latvap
    hlocp  = latvap/cpair
    dto2   = 0.5_kind_phys * deltat

#ifdef PERGRO
    mincld = 1.e-4_kind_phys
    iter = 1           ! number of times to iterate the precipitation calculation
#else
    mincld = 1.e-4_kind_phys
    iter = 2
#endif

    ! One column per thread; levels, iterations and the microphysics of a
    ! column in the original order (columns are independent).
    !$acc data copyin(dlat, pmid, pdel, zi, troplev, ttend, tn, qtend, qn, ltend, &
    !$acc&            cldice, cldliq, omega, cldn, fice, fsnow, rhdfda, rhu00, &
    !$acc&            landm, seaicef, snowh) &
    !$acc&     copyout(qme, prodprec, prodsnow, evapprec, evapsnow, evapheat, &
    !$acc&             prfzheat, meltheat, precip, snowab, ice2pr, liq2pr, liq2snow, &
    !$acc&             lsflxprc, lsflxsnw, pracwo, psacwo, psacio, &
    !$acc&             fwaut, fsaut, fracw, fsacw, fsaci, &
    !$acc&             err_code, err_k, err_l, err_val)
    !$acc parallel loop gang vector
    do i = 1, ncol
      call prognostic_cloud_water_column(i, pver, top_lev, iter, deltat, &
        mincld, cpohl, hlocp, dto2, rhofw, &
        pi, gravit, latvap, latice, cpair, &
        dlat, pmid, pdel, zi, troplev, &
        ttend, tn, qtend, qn, ltend, cldice, cldliq, &
        cldn, fice, fsnow, rhdfda, rhu00, landm, seaicef, snowh, &
        qme, prodprec, prodsnow, evapprec, evapsnow, &
        evapheat, prfzheat, meltheat, precip, snowab, &
        ice2pr, liq2pr, liq2snow, lsflxprc, lsflxsnw, &
        pracwo, psacwo, psacio, fwaut, fsaut, fracw, fsacw, fsaci, &
        err_code, err_k, err_l, err_val)
    end do
    !$acc end data

    ! Report the error the original loops (level, iteration, phase, column)
    ! met first; they stopped there.
    ierr = 0
    key_min = huge(key_min)
    do i = 1, ncol
      if (err_code(i) /= ERR_NONE) then
        key = error_order(err_code(i), err_k(i), err_l(i), iter)
        if (key < key_min) then
          key_min = key
          ierr = i
        end if
      end if
    end do
    if (ierr > 0) then
      i = ierr
      errflg = 1
      select case (err_code(i))
      case (ERR_FINDSP_CONVERGENCE)
        write(errmsg,*) 'wv_saturation::FINDSP -- not converging at i = ', i, &
             '; t, q, p ', err_val(i,1), err_val(i,2), err_val(i,3), &
             '; tsp, qsp ', err_val(i,4), err_val(i,5)
      case (ERR_FINDSP_ENTHALPY)
        write(errmsg,*) 'wv_saturation::FINDSP -- enthalpy is not conserved at i = ', i, &
             '; t, q, p ', err_val(i,1), err_val(i,2), err_val(i,3), &
             '; tsp, qsp ', err_val(i,4), err_val(i,5)
      case (ERR_EMPTY_RH_CLOUD)
        write(iulog,*) 'prognostic_cloud_water: empty rh cloud @ ', i, err_k(i)
        write(iulog,*) 'relhum, iter ', err_val(i,1), err_l(i), err_val(i,2), err_val(i,3), err_val(i,4)
        errmsg = 'prognostic_cloud_water: empty RH cloud'
      case (ERR_NEGATIVE_RCWN)
        write(iulog,*) 'prognostic_cloud_water: prob with neg rcwn1 ', err_val(i,1), err_val(i,2)
        write(iulog,*) ' cwat, qme*deltat, prodprec*deltat ', &
           err_val(i,3), err_val(i,4), err_val(i,5), err_val(i,6)
        errmsg = 'prognostic_cloud_water: negative rcwn1'
      case (ERR_NEGATIVE_RLIQ)
        write(iulog,*) 'prognostic_cloud_water: prob with neg rliq1 ', err_val(i,1)
        errmsg = 'prognostic_cloud_water: negative rliq1'
      case (ERR_NEGATIVE_RICE)
        write(iulog,*) 'prognostic_cloud_water: prob with neg rice ', err_val(i,1)
        errmsg = 'prognostic_cloud_water: negative rice'
      end select
    end if
  end subroutine prognostic_cloud_water_run

  ! Order of an error in the original loops: level, then iteration (findsp
  ! before the iterations), then phase (empty RH cloud before the negative
  ! checks); the column decides between equal orders.
  pure function error_order(code, k, l, iter) result(key)
    integer, intent(in) :: code, k, l, iter
    integer :: key
    integer :: phase

    select case (code)
    case (ERR_EMPTY_RH_CLOUD)
      phase = 1
    case (ERR_NEGATIVE_RCWN, ERR_NEGATIVE_RLIQ, ERR_NEGATIVE_RICE)
      phase = 2
    case default
      phase = 0
    end select
    key = (k*(iter + 1) + l)*3 + phase
  end function error_order

  ! prognostic_cloud_water_run for column i: all levels from top_lev, with
  ! the iterations of each level. Values of the current level are kept in
  ! scalars and stored once per level (avoids storing and reloading the same
  ! array element, which nvfortran's device code has mishandled -- see
  ! cloud_particle_sedimentation's cfdotmc_column). On an error the column
  ! stops and the error is recorded in err_* for the host to report.
  subroutine prognostic_cloud_water_column(i, pver, top_lev, iter, deltat, &
    mincld, cpohl, hlocp, dto2, rhofw, &
    pi, gravit, latvap, latice, cpair, &
    dlat, pmid, pdel, zi, troplev, &
    ttend, tn, qtend, qn, ltend, cldice, cldliq, &
    cldn, fice, fsnow, rhdfda, rhu00, landm, seaicef, snowh, &
    qme, prodprec, prodsnow, evapprec, evapsnow, &
    evapheat, prfzheat, meltheat, precip, snowab, &
    ice2pr, liq2pr, liq2snow, lsflxprc, lsflxsnw, &
    pracwo, psacwo, psacio, fwaut, fsaut, fracw, fsacw, fsaci, &
    err_code, err_k, err_l, err_val)
    !$acc routine seq
    use wv_saturation,  only: qsat, estblf, svp_to_qsat, findsp

    integer,         intent(in)    :: i, pver, top_lev, iter
    real(kind_phys), intent(in)    :: deltat, mincld, cpohl, hlocp, dto2, rhofw
    real(kind_phys), intent(in)    :: pi, gravit, latvap, latice, cpair
    real(kind_phys), intent(in)    :: dlat(:), pmid(:,:), pdel(:,:), zi(:,:)
    integer,         intent(in)    :: troplev(:)
    real(kind_phys), intent(in)    :: ttend(:,:), tn(:,:), qtend(:,:), qn(:,:), ltend(:,:)
    real(kind_phys), intent(in)    :: cldice(:,:), cldliq(:,:)
    real(kind_phys), intent(in)    :: cldn(:,:), fice(:,:), fsnow(:,:), rhdfda(:,:), rhu00(:,:)
    real(kind_phys), intent(in)    :: landm(:), seaicef(:), snowh(:)
    real(kind_phys), intent(inout) :: qme(:,:), prodprec(:,:), prodsnow(:,:)
    real(kind_phys), intent(inout) :: evapprec(:,:), evapsnow(:,:)
    real(kind_phys), intent(inout) :: evapheat(:,:), prfzheat(:,:), meltheat(:,:)
    real(kind_phys), intent(inout) :: precip(:), snowab(:)
    real(kind_phys), intent(inout) :: ice2pr(:,:), liq2pr(:,:), liq2snow(:,:)
    real(kind_phys), intent(inout) :: lsflxprc(:,:), lsflxsnw(:,:)
    real(kind_phys), intent(inout) :: pracwo(:,:), psacwo(:,:), psacio(:,:)
    real(kind_phys), intent(inout) :: fwaut(:,:), fsaut(:,:), fracw(:,:), fsacw(:,:), fsaci(:,:)
    integer,         intent(inout) :: err_code(:), err_k(:), err_l(:)
    real(kind_phys), intent(inout) :: err_val(:,:)

    integer :: k, l                                     ! Level, iteration [index]
    integer :: status                                   ! findsp status [1]

    ! Column values carried down the levels
    real(kind_phys) :: precip_c                         ! Surface precipitation accumulated [kg m-2 s-1]
    real(kind_phys) :: snowab_c                         ! Snow entering the layer [kg m-2 s-1]
    real(kind_phys) :: precab                           ! Rate of precipitation entering layer [kg m-2 s-1]
    real(kind_phys) :: cldmax                           ! Maximum cloud fraction above current level [1]

    ! Values of the current level (the original's per-column arrays)
    real(kind_phys) :: cwat                             ! Cloud water mixing ratio [kg kg-1]
    real(kind_phys) :: tsp, qsp                         ! Saturation point temperature [K], mixing ratio [kg kg-1]
    real(kind_phys) :: es, qs, gamma                    ! qsat outputs at tn
    real(kind_phys) :: relhum, relhum1                  ! Relative humidity [1]
    real(kind_phys) :: rhu_adj                          ! Adjusted relative humidity limit [1]
    real(kind_phys) :: cldm                             ! Mean cloud fraction over timestep [1]
    real(kind_phys) :: icwc                             ! In-cloud water content [kg kg-1]
    real(kind_phys) :: cwm, cwn                         ! Cloud water at midpoint, end of timestep [kg kg-1]
    real(kind_phys) :: coef                             ! Conversion timescale for condensate to rain [s-1]
    real(kind_phys) :: prprov                           ! Provisional precipitation at bottom of layer [kg m-2 s-1]
    real(kind_phys) :: calpha, cbeta, cbetah, cgamma, cgamah, rcgama, csigma
    real(kind_phys) :: qmeres, qmec1, qmec2, qmec3, qmec4
    real(kind_phys) :: qtl                              ! Saturation tendency [kg kg-1 s-1]
    real(kind_phys) :: rcwn, rliq, rice                 ! Cloud water budget checks [kg kg-1]

    ! Outputs at the current level
    real(kind_phys) :: qme_k, prodprec_k, prodsnow_k, evapprec_k, evapsnow_k
    real(kind_phys) :: evapheat_k, prfzheat_k, meltheat_k
    real(kind_phys) :: ice2pr_k, liq2pr_k, liq2snow_k
    real(kind_phys) :: pracwo_k, psacwo_k, psacio_k
    real(kind_phys) :: fwaut_k, fsaut_k, fracw_k, fsacw_k, fsaci_k

    ! Work variables
    real(kind_phys), parameter :: omsm = 0.99999_kind_phys  ! unity minus small number for rounding
    real(kind_phys) :: nice2pr, nliq2pr, nliq2snow
    real(kind_phys) :: dqsdt, qtmp, ttmp, qsn, esn, prtmp, ctmp, cdt, wtthick, pol

    err_code(i) = ERR_NONE
    err_k(i) = 0
    err_l(i) = 0

    ! initialize single level and multi level fields
    precip_c = 0.0_kind_phys
    precab = 0.0_kind_phys
    snowab_c = 0.0_kind_phys
    cldmax = 0.0_kind_phys
    do k = 1, pver
      qme(i,k)      = 0._kind_phys
      evapprec(i,k) = 0._kind_phys
      prodprec(i,k) = 0._kind_phys
      evapsnow(i,k) = 0._kind_phys
      prodsnow(i,k) = 0._kind_phys
      evapheat(i,k) = 0._kind_phys
      meltheat(i,k) = 0._kind_phys
      prfzheat(i,k) = 0._kind_phys
      ice2pr(i,k)   = 0._kind_phys
      liq2pr(i,k)   = 0._kind_phys
      liq2snow(i,k) = 0._kind_phys
      fwaut(i,k)    = 0._kind_phys
      fsaut(i,k)    = 0._kind_phys
      fracw(i,k)    = 0._kind_phys
      fsacw(i,k)    = 0._kind_phys
      fsaci(i,k)    = 0._kind_phys
      pracwo(i,k)   = 0._kind_phys
      psacwo(i,k)   = 0._kind_phys
      psacio(i,k)   = 0._kind_phys
    end do
    do k = 1, pver + 1
      lsflxprc(i,k) = 0._kind_phys
      lsflxsnw(i,k) = 0._kind_phys
    end do
    precip(i) = 0._kind_phys
    snowab(i) = 0._kind_phys

    ! find the wet bulb temp and saturation value
    ! for the provisional t and q without condensation
    mp_level_loop: do k = top_lev, pver
      ! total cloud water [kg kg-1]
      cwat = cldice(i,k) + cldliq(i,k)

      ! "True" means that ice will be taken into account.
      call findsp(qn(i,k), tn(i,k), pmid(i,k), .true., tsp, qsp, status)
      ! Currently, only 2 and 8 are treated as fatal errors (findsp_vc).
      if (status == 2 .or. status == 8) then
        if (status == 2) then
          err_code(i) = ERR_FINDSP_CONVERGENCE
        else
          err_code(i) = ERR_FINDSP_ENTHALPY
        end if
        err_k(i) = k
        err_l(i) = 0
        err_val(i,1) = tn(i,k)
        err_val(i,2) = qn(i,k)
        err_val(i,3) = pmid(i,k)
        err_val(i,4) = tsp
        err_val(i,5) = qsp
        return
      end if

      call qsat(tn(i,k), pmid(i,k), es, qs, gam=gamma)

      relhum = qn(i,k)/qs
      cldm = max(cldn(i,k),mincld)

      ! the max cloud fraction above this level
      cldmax = max(cldmax, cldm)

      ! define the coefficients for C - E calculation
      calpha = 1.0_kind_phys/qs
      cbeta  = qn(i,k)/qs**2*gamma*cpohl
      cbetah = 1.0_kind_phys/qs*gamma*cpohl
      cgamma = calpha+latvap*cbeta/cpair
      cgamah = calpha+latvap*cbetah/cpair
      rcgama = cgamma/cgamah

      if(cldm > mincld) then
        icwc = max(0._kind_phys,cwat/cldm)
      else
        icwc = 0.0_kind_phys
      endif
      ! PJR the above logic give zero icwc with nonzero cwat, dont like it!
      ! PJR generates problems with csigma
      ! PJR set the icwc to a very small number, so we can start from zero cloud cover and make some clouds
      ! icwc = max(1.e-8_kind_phys,cwat/cldm)

      ! initial guess of evaporation, will be updated within iteration
      evapprec_k = conke*(1._kind_phys - cldm)*sqrt(precab) &
                   *(1._kind_phys - min(relhum,1._kind_phys))

      ! zero qmeres before iteration for each level
      qmeres = 0.0_kind_phys

      ! calculate the cooling due to a phase change of the rainwater from above
      if (tn(i,k) >= t0) then
        meltheat_k = -latice * snowab_c * gravit/pdel(i,k)
        snowab_c = 0._kind_phys
      else
        meltheat_k = 0._kind_phys
      endif

      ! level outputs not set before the iterations start from zero
      qme_k      = 0._kind_phys
      prodprec_k = 0._kind_phys
      prodsnow_k = 0._kind_phys
      evapsnow_k = 0._kind_phys
      evapheat_k = 0._kind_phys
      prfzheat_k = 0._kind_phys
      ice2pr_k   = 0._kind_phys
      liq2pr_k   = 0._kind_phys
      liq2snow_k = 0._kind_phys
      pracwo_k   = 0._kind_phys
      psacwo_k   = 0._kind_phys
      psacio_k   = 0._kind_phys
      fwaut_k    = 0._kind_phys
      fsaut_k    = 0._kind_phys
      fracw_k    = 0._kind_phys
      fsacw_k    = 0._kind_phys
      fsaci_k    = 0._kind_phys

      ! calculate qme and formation of precip.
      !
      ! The cloud microphysics is highly nonlinear and coupled with qme
      ! Both rain processes and qme are calculated iteratively.
      qme_iter_loop: do l = 1, iter

        ! Adjust relative humidity above the tropopause in higher
        ! latitude regions (relhum_min_adj):
        ! For special WACCM/CAM-chem cases with CAM4 physics:
        ! Sets rhu to a different value poleward of +/- 50 deg latitude and
        ! levels above the tropopause if polstrat_rhmin is specified
        rhu_adj = rhu00(i,k)
        if (do_psrhmin) then
          ! assumption that up is lower index
          if ((k .lt. tropLev(i)) .and. (abs(dlat(i)) .gt. 50._kind_phys)) then
            rhu_adj = psrhmin
          endif
        endif

        ! calculation of qme has 4 scenarios
        if(relhum > rhu_adj) then
          ! 1. whole grid saturation
          if(relhum >= 0.999_kind_phys .or. cldm >= 0.999_kind_phys) then
            qme_k = (calpha*qtend(i,k)-cbetah*ttend(i,k))/cgamah
          ! 2. fractional saturation
          else
            if (rhdfda(i,k) .eq. 0._kind_phys .and. icwc .eq. 0._kind_phys) then
              err_code(i) = ERR_EMPTY_RH_CLOUD
              err_k(i) = k
              err_l(i) = l
              err_val(i,1) = relhum
              err_val(i,2) = rhu_adj
              err_val(i,3) = cldm
              err_val(i,4) = cldn(i,k)
              return
            endif

            csigma = 1.0_kind_phys/(rhdfda(i,k)+cgamma*icwc)
            qmec1 = (1.0_kind_phys-cldm)*csigma*rhdfda(i,k)
            qmec2 = cldm*calpha/cgamah+(1.0_kind_phys-rcgama*cldm)*   &
                    csigma*calpha*icwc
            qmec3 = cldm*cbetah/cgamah +  &
                    (cbeta-rcgama*cldm*cbetah)*csigma*icwc
            qmec4 = csigma*cgamma*icwc

            ! Q = C-E = -C1*Al + C2*Aq - C3*At + C4*Er
            qme_k = -qmec1*ltend(i,k) + qmec2*qtend(i,k)  &
                    -qmec3*ttend(i,k) + qmec4*evapprec_k
          endif
        ! 3. when rh < rhu00, evaporate existing cloud water
        else if(cwat > 0.0_kind_phys) then
          ! liquid water should be evaporated but not to exceed
          ! saturation point. if qn > qsp, not to evaporate cwat
          qme_k = -min(max(0._kind_phys,qsp-qn(i,k)),cwat)/deltat
        ! 4. no condensation nor evaporation
        else
          qme_k = 0.0_kind_phys
        endif

        ! Because of the finite time step, place a bound here not to exceed wet bulb point
        ! and not to evaporate more than available water
        qtmp = qn(i,k) - qme_k*deltat

        ! possibilities to have qtmp > qsp
        !
        !   1. if qn > qs(tn), it condenses;
        !      if after applying qme,  qtmp > qsp,  more condensation is applied.
        !
        !   2. if qn < qs, evaporation should not exceed qsp,
        if(qtmp > qsp) then
          qme_k = qme_k + (qtmp-qsp)/deltat
        endif

        ! if net evaporation, it should not exceed available cwat
        if(qme_k < -cwat/deltat) then
           qme_k = -cwat/deltat
        endif

        ! addition of residual condensation from previous step of iteration
        qme_k = qme_k + qmeres

        ! limit qme for roundoff errors (multiply by slightly less than unity)
        qme_k = qme_k * omsm

        ! as a safe limit, condensation should not reduce grid mean rh below rhu00
        if(qme_k > 0.0_kind_phys .and. relhum > rhu_adj) then
          qme_k = min(qme_k, (qn(i,k)-qs*rhu_adj)/deltat)
        endif

        ! initial guess for cwm (mean cloud water over time step) if 1st iteration
        if(l < 2) then
          cwm = max(cwat+qme_k*dto2, 0._kind_phys)
        endif

        ! provisional precipitation falling through model layer
        ! prprov =  precab + prodprec_k*pdel(i,k)/gravit
        ! rain produced in this layer not too effective in collection process
        wtthick = max(0._kind_phys,min(0.5_kind_phys,((zi(i,k)-zi(i,k+1))/1000._kind_phys)**2))
        prprov = precab + wtthick*prodprec_k*pdel(i,k)/gravit

        ! calculate conversion of condensate to precipitation by cloud microphysics
        call findmcnew( &
          pi      = pi,   &
          precab  = prprov, &
          snowab  = snowab_c, &
          t       = tn(i,k), &
          p       = pmid(i,k), &
          p_bot   = pmid(i,pver), &
          cwm     = cwm, &
          cldm    = cldm, &
          cldmax  = cldmax, &
          fice    = fice(i,k),  &
          landm   = landm(i),   &
          seaicef = seaicef(i), &
          snowh   = snowh(i),   &
          ! below output
          coef    = coef, &
          fwaut   = fwaut_k, &
          fsaut   = fsaut_k, &
          fracw   = fracw_k, &
          fsacw   = fsacw_k, &
          fsaci   = fsaci_k, &
          pracwo  = pracwo_k,&
          psacwo  = psacwo_k,&
          psacio  = psacio_k)

        ! calculate the precip rate
        if (cldm > 0) then
          ! first predict the cloud water
          cdt = coef*deltat
          if(cdt > 0.01_kind_phys) then
            pol = qme_k/coef ! production over loss
            cwn = max(0._kind_phys,(cwat-pol)*exp(-cdt)+ pol)
          else
            cwn = max(0._kind_phys,(cwat + qme_k*deltat)/(1+cdt))
          endif

          ! now back out the tendency of net rain production
          prodprec_k = max(0._kind_phys,qme_k-(cwn-cwat)/deltat)
        else
          prodprec_k = 0.0_kind_phys
          cwn = 0._kind_phys
        endif

        ! provisional calculation of conversion terms
        ice2pr_k = prodprec_k*(fsaut_k+fsaci_k)
        liq2pr_k = prodprec_k*(fwaut_k+fsacw_k+fracw_k)

        ! revision suggested by Jim McCaa
        ! it controls the amount of snow hitting the sfc
        ! by forcing a lot of conversion of cloud liquid to snow phase
        ! it might be better done later by an explicit representation of
        ! rain accreting ice (and freezing), or by an explicit freezing of raindrops
        !
        ! old:
        ! liq2snow_k = prodprec_k*fsacw_k
        liq2snow_k = max(prodprec_k*fsacw_k, fsnow(i,k)*liq2pr_k)

        ! bounds
        nice2pr = min(ice2pr_k,(cwat+qme_k*deltat)*fice(i,k)/deltat)
        nliq2pr = min(liq2pr_k,(cwat+qme_k*deltat)*(1._kind_phys-fice(i,k))/deltat)

        if (liq2pr_k .ne. 0._kind_phys) then
          nliq2snow = liq2snow_k*nliq2pr/liq2pr_k   ! correction
        else
          nliq2snow = liq2snow_k
        endif

        ! avoid roundoff problems generating negatives
        nliq2snow = nliq2snow*omsm
        nliq2pr   = nliq2pr*omsm
        nice2pr   = nice2pr*omsm

        ! final estimates of conversion to precip and snow
        prodprec_k = (nliq2pr + nice2pr)
        prodsnow_k = (nice2pr + nliq2snow)

        ! compute diagnostics and sanity checks
        rcwn =  cwat + (qme_k - prodprec_k)*deltat
        rliq = (cwat + qme_k*deltat)*(1._kind_phys-fice(i,k)) - nliq2pr*deltat
        rice = (cwat + qme_k*deltat) * fice(i,k) - nice2pr*deltat

        ! Sanity checks
        if(abs(rcwn) < 1.e-300_kind_phys) rcwn = 0._kind_phys
        if(abs(rliq) < 1.e-300_kind_phys) rliq = 0._kind_phys
        if(abs(rice) < 1.e-300_kind_phys) rice = 0._kind_phys
        if(rcwn < 0._kind_phys) then
          err_code(i) = ERR_NEGATIVE_RCWN
          err_val(i,1) = rcwn
          err_val(i,2) = cwn
          err_val(i,3) = cwat
          err_val(i,4) = qme_k*deltat
          err_val(i,5) = prodprec_k*deltat
          err_val(i,6) = (qme_k-prodprec_k)*deltat
        else if (rliq < 0._kind_phys) then
          err_code(i) = ERR_NEGATIVE_RLIQ
          err_val(i,1) = rliq
        else if (rice < 0._kind_phys) then
          err_code(i) = ERR_NEGATIVE_RICE
          err_val(i,1) = rice
        endif
        if (err_code(i) /= ERR_NONE) then
          err_k(i) = k
          err_l(i) = l
          return
        endif

        ! Final version of condensate to precip terms
        liq2pr_k = nliq2pr
        liq2snow_k = nliq2snow
        ice2pr_k = nice2pr
        cwn = rcwn

        ! update any remaining provisional values
        cwm = (cwn + cwat)*0.5_kind_phys

        ! update in cloud water
        if(cldm > mincld) then
          icwc = cwm/cldm
        else
          icwc = 0.0_kind_phys
        endif
        ! PJR the above logic give zero icwc with nonzero cwat, dont like it!
        ! PJR generates problems with csigma
        ! PJR set the icwc to a very small number, so we can start from zero cloud cover and make some clouds
        ! icwc = max(1.e-8_kind_phys,cwm/cldm)

        ! calculate provisional value of cloud water for
        ! evaporation of precipitate (evapprec) calculation
        qtmp = qn(i,k) - qme_k*deltat
        ttmp = tn(i,k) + deltat/cpair * ( meltheat_k + (latvap + latice*fice(i,k)) * qme_k )
        esn = estblf(ttmp)
        qsn = svp_to_qsat(esn, pmid(i,k))
        qtl = max((qsn - qtmp)/deltat,0._kind_phys)
        relhum1 = qtmp/qsn

#ifdef PERGRO
        evapprec_k = conke*(1._kind_phys - max(cldm,mincld))* &
                     sqrt(precab)*(1._kind_phys - min(relhum1,1._kind_phys))
#else
        evapprec_k = conke*(1._kind_phys - cldm)*sqrt(precab) &
                     *(1._kind_phys - min(relhum1,1._kind_phys))
#endif

        ! limit the evaporation to the amount which is entering the box
        ! or saturates the box
        prtmp = precab*gravit/pdel(i,k)
        evapprec_k = min(evapprec_k, prtmp, qtl)*omsm
#ifdef PERGRO
        ! zeroing needed for pert growth
        evapprec_k = 0._kind_phys
#endif

        ! Partition evaporation of precipitate between rain and snow using
        ! the fraction of snow falling into the box. Determine the heating
        ! due to evaporation. Note that evaporation is positive (loss of precip,
        ! gain of vapor) and that heating is negative.
        if (evapprec_k > 0._kind_phys) then
          evapsnow_k = evapprec_k * snowab_c / precab
          evapheat_k = -latvap * evapprec_k - latice * evapsnow_k
        else
          evapsnow_k = 0._kind_phys
          evapheat_k = 0._kind_phys
        end if

        ! Account for the latent heat of fusion for liquid drops collected by falling snow
        prfzheat_k = latice * liq2snow_k

        ! now remove the residual of any over-saturation. Normally,
        ! the oversaturated water vapor should have been removed by
        ! qme formulation plus constraints by wet bulb tsp/qsp
        ! as computed above. However, because of non-linearity,
        ! addition of (qme-evapprec) to update t and q may still cause
        ! a very small amount of over saturation. It is called a
        ! residual of over-saturation because theoretically, qme
        ! should have taken care of all of large scale condensation.
        qtmp = qn(i,k)-(qme_k-evapprec_k)*deltat
        ttmp = tn(i,k) + deltat/cpair * ( meltheat_k + evapheat_k + prfzheat_k      &
               + (latvap + latice*fice(i,k)) * qme_k )

        call qsat(ttmp, pmid(i,k), esn, qsn, dqsdt=dqsdt)

        if(qtmp > qsn) then
          ! now extra condensation to bring air to just saturation
          ctmp = (qtmp-qsn)/(1._kind_phys+hlocp*dqsdt)/deltat
          qme_k = qme_k+ctmp
          ! save residual on qmeres to addtion to qme on entering next iteration
          ! qme exit here contain the residual but overrided if back to iteration
          qmeres = ctmp
        else
          qmeres = 0.0_kind_phys
        endif
      end do qme_iter_loop  ! loop over l (iteration)

      ! precipitation
      precip_c = precip_c + pdel(i,k)/gravit * (prodprec_k - evapprec_k)
      precab = precab + pdel(i,k)/gravit * (prodprec_k - evapprec_k)
      if(precab < 0._kind_phys) then
        precab = 0._kind_phys
      endif

      ! snowab_c = snowab_c + pdel(i,k)/gravit * (prodprec_k*fice(i,k) - evapsnow_k)
      snowab_c = snowab_c + pdel(i,k)/gravit * (prodsnow_k - evapsnow_k)

      ! If temperature above freezing, all precip is rain flux. if temperature below freezing, all precip is snow flux.
      lsflxprc(i,k+1) = precab   !! making this consistent with other precip fluxes.  prc = rain + snow
      ! lsflxprc(i,k+1) = precab - snowab_c
      lsflxsnw(i,k+1) = snowab_c

      ! store the outputs of this level
      qme(i,k)      = qme_k
      prodprec(i,k) = prodprec_k
      prodsnow(i,k) = prodsnow_k
      evapprec(i,k) = evapprec_k
      evapsnow(i,k) = evapsnow_k
      evapheat(i,k) = evapheat_k
      prfzheat(i,k) = prfzheat_k
      meltheat(i,k) = meltheat_k
      ice2pr(i,k)   = ice2pr_k
      liq2pr(i,k)   = liq2pr_k
      liq2snow(i,k) = liq2snow_k
      pracwo(i,k)   = pracwo_k
      psacwo(i,k)   = psacwo_k
      psacio(i,k)   = psacio_k
      fwaut(i,k)    = fwaut_k
      fsaut(i,k)    = fsaut_k
      fracw(i,k)    = fracw_k
      fsacw(i,k)    = fsacw_k
      fsaci(i,k)    = fsaci_k
    end do mp_level_loop    ! loop over k (level)

    ! Convert precip (lwe_stratiform_precipitation_rate_at_surface)
    ! and snowab (lwe_snow_precipitation_rate_at_surface_due_to_microphysics)
    ! from kg m-2 s-1 to m s-1 (precipitation units) for surface exchange.
    !
    ! If this conversion is removed in the future, the metadata needs to
    ! be updated.
    precip(i) = precip_c/rhofw
    snowab(i) = snowab_c/rhofw
  end subroutine prognostic_cloud_water_column

  ! Calculate the conversion of condensate to precipitate at one point
  ! (one column at level k; p_bot is the pressure of the bottom level).
  ! Rasch, P. J., and J. E. Kristjánsson, 1998: A Comparison of the CCM3 Model Climate Using Diagnosed and Predicted Condensate Parameterizations. J. Climate, 11, 1587–1614
  ! https://doi.org/10.1175/1520-0442(1998)011<1587:ACOTCM>2.0.CO;2
  !
  ! Original author: P. Rasch
  subroutine findmcnew( &
    pi, &
    precab, snowab, t, p, p_bot, cwm, cldm, cldmax, fice, &
    landm, seaicef, snowh, &
    coef, fwaut, fsaut, &
    fracw, fsacw, fsaci, &
    pracwo, psacwo, psacio)
    !$acc routine seq

    ! input arguments
    real(kind_phys), intent(in) :: pi               ! pi_constant [1]

    real(kind_phys), intent(in) :: precab           ! rate of precipitation from above [kg m-2 s-1]
    real(kind_phys), intent(in) :: t                ! temperature [K]
    real(kind_phys), intent(in) :: p                ! air_pressure [Pa]
    real(kind_phys), intent(in) :: p_bot            ! air_pressure at the bottom level [Pa]
    real(kind_phys), intent(in) :: cwm              ! condensate mixing ratio [kg kg-1]
    real(kind_phys), intent(in) :: cldm             ! cloud fraction
    real(kind_phys), intent(in) :: cldmax           ! max cloud fraction above this level
    real(kind_phys), intent(in) :: fice             ! fraction of cwat that is ice
    real(kind_phys), intent(in) :: landm            ! Land fraction ramped over water
    real(kind_phys), intent(in) :: seaicef          ! sea ice fraction
    real(kind_phys), intent(in) :: snowab           ! rate of snow from above [kg m-2 s-1]
    real(kind_phys), intent(in) :: snowh            ! Snow depth over land, water equivalent [m]

    ! output arguments
    real(kind_phys), intent(out) :: coef            ! conversion rate [s-1]
    real(kind_phys), intent(out) :: fwaut           ! relative importance of liquid autoconversion (a diagnostic)
    real(kind_phys), intent(out) :: fsaut           ! relative importance of ice autoconversion    (a diagnostic)
    real(kind_phys), intent(out) :: fracw           ! relative importance of rain accreting liquid (a diagnostic)
    real(kind_phys), intent(out) :: fsacw           ! relative importance of snow accreting liquid (a diagnostic)
    real(kind_phys), intent(out) :: fsaci           ! relative importance of snow accreting ice    (a diagnostic)
    ! not set where there is no microphysics (cwm <= 1.e-20), as in the
    ! original
    real(kind_phys), intent(inout) :: pracwo        ! accretion of cloud water by rain [s-1]
    real(kind_phys), intent(inout) :: psacwo        ! accretion of cloud water by snow [s-1]
    real(kind_phys), intent(inout) :: psacio        ! accretion of cloud ice by snow [s-1]

    ! local variables
    real(kind_phys) :: capn                         ! Local cloud particle number concentration [cm-3]
    real(kind_phys) :: cldloc                       ! Non-zero cloud fraction [1]
    real(kind_phys) :: cldpr                        ! Cloud fraction for precipitation [1]
    real(kind_phys) :: totmr                        ! In-cloud total water mixing ratio [kg kg-1]
    real(kind_phys) :: icemr                        ! In-cloud ice mixing ratio [kg kg-1]
    real(kind_phys) :: liqmr                        ! In-cloud liquid water mixing ratio [kg kg-1]
    real(kind_phys) :: rainmr                       ! In-cloud rain mixing ratio [kg kg-1]
    real(kind_phys) :: ciaut                        ! Ice autoconversion coefficient [s-1]
    real(kind_phys) :: pracw                        ! Rate of rain collecting cloud water [s-1]
    real(kind_phys) :: psaci                        ! Rate of snow collecting cloud ice [s-1]
    real(kind_phys) :: psacw                        ! Rate of snow collecting cloud water [s-1]
    real(kind_phys) :: psaut                        ! Rate of ice autoconversion [s-1]
    real(kind_phys) :: pwaut                        ! Rate of liquid water autoconversion [s-1]
    real(kind_phys) :: ptot                         ! Total conversion rate [s-1]
    real(kind_phys) :: prlloc                       ! Local rain flux [mm day-1]
    real(kind_phys) :: prscgs                       ! Local snow amount [g cm-2]
    real(kind_phys) :: snowfr                       ! Snow fraction of precipitation [1]
    real(kind_phys) :: vfallw                       ! Fall speed of liquid precipitation [m s-1]
    real(kind_phys) :: rho                          ! Air density [kg m-3]
    real(kind_phys) :: rhocgs                       ! Air density in CGS units [g cm-3]
    real(kind_phys) :: r3l                          ! Cloud droplet volume radius [m]
    real(kind_phys) :: icrit                        ! Ice autoconversion threshold [kg kg-1]
    real(kind_phys) :: wt                           ! Ice fraction weight [1]
    real(kind_phys) :: con1                         ! Work constant for radius calculation [m]
    real(kind_phys) :: con2                         ! Work constant for density ratios [1]
    real(kind_phys) :: csacx                        ! Constant used for snow accreting liquid or ice [??]
    real(kind_phys) :: rat1                         ! Density ratio work variable [1]
    real(kind_phys) :: rat2                         ! Mixing ratio ratio work variable [1]

    ! set the output variables to zero; the microphysics is done only
    ! where there is condensate
    coef = 0._kind_phys
    fwaut = 0._kind_phys
    fsaut = 0._kind_phys
    fracw = 0._kind_phys
    fsacw = 0._kind_phys
    fsaci = 0._kind_phys
    liqmr = 0._kind_phys
    rainmr = 0._kind_phys
    if (.not. (cwm > 1.e-20_kind_phys)) return

    ! the local cloudiness at this level
    cldloc = max(cldmin,cldm)

    ! a weighted mean between max cloudiness above, and this layer
    cldpr = max(cldmin,(cldmax+cldm)*0.5_kind_phys)

    ! decompose the suspended condensate into
    ! an incloud liquid and ice phase component
    totmr = cwm/cldloc
    icemr = totmr*fice
    liqmr = totmr*(1._kind_phys-fice)

    ! density
    rho = p/(287._kind_phys*t)
    rhocgs = rho*1.e-3_kind_phys     ! density in cgs units

    ! decompose the precipitate into a liquid and ice phase
    if (t > t0) then
       vfallw = convfw/sqrt(rho)
       rainmr = precab/(rho*vfallw*cldpr)
       snowfr = 0
    else
       snowfr = 1
       rainmr = 0._kind_phys
    endif

    ! local snow amount in cgs units
    prscgs = precab/cldpr*0.1_kind_phys*snowfr

    ! local rain amount in mm/day
    prlloc = precab*86400._kind_phys/cldpr

    con1 = 1._kind_phys/(1.333_kind_phys*pi)**0.333_kind_phys * 0.01_kind_phys ! meters

    ! calculate the conversion terms
    rhocgs = rho*1.e-3_kind_phys     ! density in cgs units

    ! some temperature dependent constants
    wt = fice
    icrit = icritc*wt + icritw*(1-wt)

    ! jrm Reworked droplet number concentration algorithm
    ! Start with pressure-dependent value appropriate for continental air
    ! Note: reltab has a temperature dependence here
    capn = capnw + (capnc-capnw) * min(1._kind_phys,max(0._kind_phys,1.0_kind_phys-(p-0.8_kind_phys*p_bot)/(0.2_kind_phys*p_bot)))
    ! Modify for snow depth over land
    capn = capn + (capnc-capn) * min(1.0_kind_phys,max(0.0_kind_phys,snowh*10._kind_phys))
    ! Ramp between polluted value over land to clean value over ocean.
    capn = capn + (capnc-capn) * min(1.0_kind_phys,max(0.0_kind_phys,1.0_kind_phys-landm))
    ! Ramp between the resultant value and a sea ice value in the presence of ice.
    capn = capn + (capnsi-capn) * min(1.0_kind_phys,max(0.0_kind_phys,seaicef))
    ! end jrm

    ! useful terms in following calculations
    rat1 = rhocgs/rhow
    rat2 = liqmr/capn
    con2 = (rat1*rat2)**0.333_kind_phys

    ! volume radius
    r3l = con1*con2

    ! critical threshold for autoconversion if modified for mixed phase
    ! clouds to mimic a bergeron findeisen process
    ! r3lc2 = r3lcrit*(1.-0.5*fice*(1-fice))
    !
    ! autoconversion of liquid
    !
    !        cwaut = 2.e-4
    !        cwaut = 1.e-3
    !        lcrit = 2.e-4
    !        lcrit = 5.e-4
    !        pwaut = max(0._kind_phys,liqmr-lcrit)*cwaut
    !
    ! pwaut is following tripoli and cotton (and many others)
    ! we reduce the autoconversion below critpr, because these are regions where
    ! the drop size distribution is likely to imply much smaller collector drops than
    ! those relevant for a cloud distribution corresponding to the value of effc = 0.55
    ! suggested by cotton (see austin 1995 JAS, baker 1993)

    ! easy to follow form
    !        pwaut = capc*liqmr**2*rhocgs/rhow
    !                *(liqmr*rhocgs/(rhow*capn))**(.333)
    !                *heavy(r3l,r3lcrit)
    !                *max(0.10_kind_phys,min(1._kind_phys,prlloc/critpr))
    ! somewhat faster form

    ! using modified Heaviside function:
    pwaut = capc*liqmr**2*rat1*con2*heavymp(r3l,r3lcrit) * &
            max(0.10_kind_phys,min(1._kind_phys,prlloc/critpr))

    ! not using modified Heaviside function:
    ! pwaut = capc*liqmr**2*rat1*con2*heavym(r3l,r3lcrit)* &
    !         max(0.10_kind_phys,min(1._kind_phys,prlloc/critpr))

    ! autoconversion of ice
    ciaut = ciautb

    ! autoconversion of ice condensate
#ifdef PERGRO
    psaut = heavyp(icemr,icrit)*icemr*ciaut
#else
    psaut = max(0._kind_phys,icemr-icrit)*ciaut
#endif

    ! collection of liquid by rain
    ! pracw = cracw*rho*liqmr*rainmr !(beheng 1994)
    pracw = cracw*rho*sqrt(rho)*liqmr*rainmr !(tripoli and cotton)
    pracwo = pracw

    ! the following lines calculate the slope parameter and snow mixing ratio
    ! from the precip rate using the equations found in lin et al 83
    ! in the most natural form, but it is expensive, so after some tedious
    ! algebraic manipulation you can use the cheaper form found below
    !            vfalls = c*gam4pd/(6*lamdas**d)*sqrt(rhonot/rhocgs)
    !     $               *0.01   ! convert from cm/s to m/s
    !            snowmr = snowfr*precab/(rho*vfalls*cldpr)
    !            snowmr = ( prscgs*mcon02 * (rhocgs**mcon03) )**mcon04
    !            lamdas = (prhonos/max(rhocgs*snowmr,small))**0.25
    !            csacw = mcon01*sqrt(rhonot/rhocgs)/(lamdas**thrpd)
    !
    ! coefficient for collection by snow independent of phase
    csacx = mcon07*rhocgs**mcon08*prscgs**mcon05

    ! collection of liquid by snow (lin et al 1983)
    psacw = csacx*liqmr*esw

#ifdef PERGRO
    ! this is necessary for pergro
    psacw = 0._kind_phys
#endif

    psacwo = psacw

    ! collection of ice by snow (lin et al 1983)
    psaci = csacx*icemr*esi
    psacio = psaci

    ! total conversion of condensate to precipitate
    ptot = pwaut + psaut + pracw + psacw + psaci

    ! the reciprocal of cloud water amount (or zero if no cloud water)
    !  rcwm =  totmr/(max(totmr,small)**2)

    ! turn the tendency back into a loss rate [s-1]
    if (totmr > 0._kind_phys) then
       coef = ptot/totmr
    else
       coef = 0._kind_phys
    endif

    if (ptot.gt.0._kind_phys) then
       fwaut = pwaut/ptot
       fsaut = psaut/ptot
       fracw = pracw/ptot
       fsacw = psacw/ptot
       fsaci = psaci/ptot
    else
       fwaut = 0._kind_phys
       fsaut = 0._kind_phys
       fracw = 0._kind_phys
       fsacw = 0._kind_phys
       fsaci = 0._kind_phys
    endif

  end subroutine findmcnew

  ! Heavyside step functions used in findmcnew
  ! Heavyside step function
  pure function heavy(a1, a2) result(h)
    !$acc routine seq
    ! h: 1 if a1 > a2, 0 otherwise
    real(kind_phys), intent(in) :: a1, a2
    real(kind_phys) :: h
    h = max(0.0_kind_phys, sign(1.0_kind_phys, a1-a2))
  end function heavy

  ! Modified Heavyside function with minimum value of 0.01
  pure function heavym(a1, a2) result(h)
    !$acc routine seq
    ! h: 1 if a1 > a2, 0.01 otherwise
    real(kind_phys), intent(in) :: a1, a2
    real(kind_phys) :: h
    h = max(0.01_kind_phys, sign(1.0_kind_phys, a1-a2))
  end function heavym

  ! Smoothed Heavyside function using ratio formulation
  pure function heavyp(a1, a2) result(h)
    !$acc routine seq
    ! h: Ratio of a1/(a1+a2+eps)
    real(kind_phys), intent(in) :: a1, a2
    real(kind_phys) :: h
    real(kind_phys), parameter :: eps = 1.0e-36_kind_phys
    h = a1/(a2 + a1 + eps)
  end function heavyp

  ! Modified smooth Heavyside with offset term
  pure function heavymp(a1, a2) result(h)
    !$acc routine seq
    ! h: Ratio (a1 + 0.01*a2)/(a1+a2+eps)
    real(kind_phys), intent(in) :: a1, a2
    real(kind_phys) :: h
    real(kind_phys), parameter :: eps = 1.0e-36_kind_phys
    h = (a1 + 0.01_kind_phys*a2)/(a2 + a1 + eps)
  end function heavymp

end module prognostic_cloud_water
