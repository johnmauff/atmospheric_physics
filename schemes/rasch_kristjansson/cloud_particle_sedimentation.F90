! Compute tendencies from sedimentation of cloud liquid and ice particles
! Original authors: Byron Boville, September 2002 from code by P.J. Rasch
! CCPP-ized: Haipeng Lin, January 2025
module cloud_particle_sedimentation

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: cloud_particle_sedimentation_init
  public :: cloud_particle_sedimentation_run

  ! tuning parameters for cloud liquid and ice particle sedimentation
  real(kind_phys), parameter :: vland = 1.5_kind_phys        ! liquid fall velocity over land [cm s-1]
  real(kind_phys), parameter :: vocean = 2.8_kind_phys       ! liquid fall velocity over ocean [cm s-1]
  real(kind_phys), parameter :: mxsedfac = 0.99_kind_phys    ! maximum sedimentation flux factor

  ! use Stokes velocity method
  ! or McFarquhar and Heymsfield (https://doi.org/10.1175/1520-0469(1997)054<2187:POTCIC>2.0.CO;2)
  logical,         parameter :: stokes = .true.              ! use Stokes velocity instead of McFarquhar and Heymsfield

  ! tuning parameters for Stokes terminal velocity
  real(kind_phys)            :: cldsed_ice_stokes_fac        ! factor applied to ice fall vel
                                                             ! from Stokes terminal velocity
  real(kind_phys), parameter :: eta  = 1.7e-5_kind_phys      ! viscosity of air [kg m s-1]
  real(kind_phys), parameter :: r40  = 40._kind_phys         !  40 micron radius
  real(kind_phys), parameter :: r400 = 400._kind_phys        ! 400 micron radius
  real(kind_phys), parameter :: v400 = 1.00_kind_phys        ! fall velocity of 400 micron sphere [m s-1]
  real(kind_phys)            :: v40                          ! Stokes fall velocity of 40 micron sphere [m s-1]
  real(kind_phys)            :: vslope                       ! linear slope for large particles [m s-1 micron-1]

  ! tuning parameters for McFarquhar and Heymsfield terminal velocity
  real(kind_phys), parameter :: vice_small = 1._kind_phys    ! ice fall velocity for small concentration [cm s-1]
  real(kind_phys)            :: lbound                       ! lower bound for iciwc
  !$acc declare create(cldsed_ice_stokes_fac, v40, vslope, lbound)

  ! misc
  real(kind_phys), parameter :: um_to_m = 1.e-12_kind_phys

contains

!> \section arg_table_cloud_particle_sedimentation_init Argument Table
!! \htmlinclude arg_table_cloud_particle_sedimentation_init.html
  subroutine cloud_particle_sedimentation_init(&
    amIRoot, iulog, &
    cldsed_ice_stokes_fac_in, &
    rhoh2o, gravit, &
    errmsg, errflg)

    ! Input arguments
    logical,            intent(in)    :: amIRoot
    integer,            intent(in)    :: iulog          ! log output unit
    real(kind_phys),    intent(in)    :: cldsed_ice_stokes_fac_in
    real(kind_phys),    intent(in)    :: gravit
    real(kind_phys),    intent(in)    :: rhoh2o

    ! Output arguments
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local variables
    real(kind_phys)     :: ac, bc, cc                   ! constants for McFarquhar and Heymsfield method

    errmsg = ''
    errflg = 0
    cldsed_ice_stokes_fac = cldsed_ice_stokes_fac_in

    ! linear ramp variables for fall velocity
    v40 = (2._kind_phys/9._kind_phys)*rhoh2o*gravit/eta*r40**2*1.e-12_kind_phys
    vslope = (v400 - v40)/(r400 - r40)

    ! McFarquhar and Heymsfield lower bound for iciwc
    cc = 128.64_kind_phys
    bc = 53.242_kind_phys
    ac = 5.4795_kind_phys
    lbound = (-bc + sqrt(bc*bc - 4*ac*cc))/(2*ac)
    lbound = 10._kind_phys**lbound
    !$acc update device(cldsed_ice_stokes_fac, v40, vslope, lbound)

    if(amIRoot) then
      write(iulog,*) 'cloud_particle_sedimentation_init: cldsed_ice_stokes_fac = ', cldsed_ice_stokes_fac
    endif

  end subroutine cloud_particle_sedimentation_init

  ! Compute gravitational sedimentation velocities for cloud liquid water
  ! and ice, based on Lawrence and Crutzen (1998).
  ! https://doi.org/10.3402/tellusb.v50i3.16129
  !
  ! and apply cloud particle gravitational sedimentation to condensate.
  !
  ! Original authors: mgl, March 1998 (MATCH-MPIC version 2.0)
  !                   adapted by P.J. Rasch
  !                   B.A. Boville, September 2002
  !                   P.J. Rasch, May 2003
!> \section arg_table_cloud_particle_sedimentation_run Argument Table
!! \htmlinclude arg_table_cloud_particle_sedimentation_run.html
  subroutine cloud_particle_sedimentation_run( &
    ncol, pver, pverp, dtime, &
    tmelt, gravit, latvap, latice, rair, rhoh2o, &
    icritc, & ! from prognostic_cloud_water namelist
    pint, pmid, pdel, t, cloud, &
    icefrac, landfrac, ocnfrac, &
    cldliq, cldice, snowh, landm, &
    ! below output for sedimentation velocity
    pvliq, pvice, &
    ! below output for sedimentation to condensate
    liqtend, icetend, wvtend, htend, sfliq, sfice, &
    errmsg, errflg)
    ! To-be-CCPPized dependencies
    use cloud_optical_properties,     only: rel_point, rei_point

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    integer,            intent(in)    :: pverp
    real(kind_phys),    intent(in)    :: dtime
    real(kind_phys),    intent(in)    :: tmelt
    real(kind_phys),    intent(in)    :: gravit
    real(kind_phys),    intent(in)    :: latvap
    real(kind_phys),    intent(in)    :: latice
    real(kind_phys),    intent(in)    :: rair
    real(kind_phys),    intent(in)    :: rhoh2o
    real(kind_phys),    intent(in)    :: icritc         ! tunable_parameter_for_autoconversion_of_cold_ice_for_rk_microphysics [kg kg-1]
    real(kind_phys),    intent(in)    :: pint(:,:)      ! air_pressure_at_interface [Pa]
    real(kind_phys),    intent(in)    :: pmid(:,:)      ! air_pressure [Pa]
    real(kind_phys),    intent(in)    :: pdel(:,:)      ! air_pressure_thickness [Pa]
    real(kind_phys),    intent(in)    :: t(:,:)         ! air_temperature [K]
    real(kind_phys),    intent(in)    :: cloud(:,:)     ! cloud_area_fraction [fraction]
    real(kind_phys),    intent(in)    :: icefrac(:)     ! sea_ice_area_fraction_from_coupler [fraction]
    real(kind_phys),    intent(in)    :: landfrac(:)    ! land_area_fraction_from_coupler [fraction]
    real(kind_phys),    intent(in)    :: ocnfrac(:)     ! ocean_area_fraction [fraction]
    real(kind_phys),    intent(in)    :: cldliq(:,:)    ! adv: cloud_liquid_water_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: cldice(:,:)    ! adv: cloud_ice_mixing_ratio_wrt_moist_air_and_condensed_water [kg kg-1]
    real(kind_phys),    intent(in)    :: snowh(:)       ! lwe_surface_snow_depth_over_land_from_coupler [m]
    real(kind_phys),    intent(in)    :: landm(:)       ! smoothed_land_area_fraction [fraction]

    ! Output arguments
    ! note: pvliq, pvice are at the interfaces (loss from cell is based on pvel(k+1))
    real(kind_phys),    intent(out)   :: pvliq(:,:)     ! vertical velocity of cloud liquid drops [Pa s-1]
    real(kind_phys),    intent(out)   :: pvice(:,:)     ! vertical velocity of cloud ice particles [Pa s-1]

    real(kind_phys),    intent(out)   :: liqtend(:,:)   ! liquid condensate tendency [kg kg-1 s-1] -- to apply cldliq tendency
    real(kind_phys),    intent(out)   :: icetend(:,:)   ! ice condensate tendency [kg kg-1 s-1] -- to apply cldice tendency
    real(kind_phys),    intent(out)   :: wvtend(:,:)    ! water vapor tendency [kg kg-1 s-1] -- to apply wv tendency
    real(kind_phys),    intent(out)   :: htend(:,:)     ! heating rate [J kg-1 s-1] -- to apply s tendency
    real(kind_phys),    intent(out)   :: sfliq(:)       ! surface flux of liquid (rain) [kg m-2 s-1]
    real(kind_phys),    intent(out)   :: sfice(:)       ! stratiform_lwe_cloud_ice_surface_flux_due_to_sedimentation [m s-1]
    character(len=512), intent(out)   :: errmsg         ! error message
    integer,            intent(out)   :: errflg         ! error flag

    ! Local variables
    real(kind_phys)     :: rho                          ! air density [kg m-3]
    real(kind_phys)     :: vfall                        ! settling velocity of cloud particles [m s-1]
    real(kind_phys)     :: icice                        ! in cloud ice water content [kg kg-1]
    real(kind_phys)     :: iciwc                        ! in cloud ice water content [g m-3]
    real(kind_phys)     :: icefac
    real(kind_phys)     :: logiwc
    real(kind_phys)     :: rei                          ! effective radius of ice particles [um]
    real(kind_phys)     :: rel                          ! effective radius of liquid particles [um]
    real(kind_phys)     :: fxliq(ncol, pverp)           ! fluxes at interface, liquid (positive is down) [Pa]
    real(kind_phys)     :: fxice(ncol, pverp)           ! fluxes at interface, ice (positive is down) [Pa]
    real(kind_phys)     :: cldab                        ! cloud in layer above [fraction]
    real(kind_phys)     :: evapliq                      ! evaporation of cloud liquid into environment [kg kg-1 s-1]
    real(kind_phys)     :: evapice                      ! evaporation of cloud ice into environment [kg kg-1 s-1]
    real(kind_phys)     :: cldovrl                      ! cloud overlap factor
    integer             :: i, k

    ! Work arrays of getflx_column (and cfdotmc_column), one column each
    real(kind_phys)     :: psi(ncol, pverp)             ! integral of phi
    real(kind_phys)     :: fdot(ncol, pverp)            ! derivative of psi at nodes
    real(kind_phys)     :: s(ncol, pverp)               ! first divided differences at nodes
    real(kind_phys)     :: sh(ncol, pverp)              ! first divided differences between nodes
    real(kind_phys)     :: d(ncol, pverp)               ! second divided differences at nodes
    real(kind_phys)     :: dh(ncol, pverp)              ! second divided differences between nodes
    real(kind_phys)     :: e(ncol, pverp)               ! third divided differences at nodes
    real(kind_phys)     :: eh(ncol, pverp)              ! third divided differences between nodes
    real(kind_phys)     :: delxh(ncol, pverp)           ! node spacing

    ! cfint2 errors: the last column whose interval was not found (0 if none)
    logical             :: found
    integer             :: fail_liq, fail_ice

    errmsg = ''
    errflg = 0
    fail_liq = 0
    fail_ice = 0

    ! One column per thread; levels in order within a column. The unused
    ! inputs (landfrac, ocnfrac) are in the data clause so they stay device
    ! arguments.
    !$acc data copyin(pint, pmid, pdel, t, cloud, icefrac, landfrac, ocnfrac, &
    !$acc&            cldliq, cldice, snowh, landm) &
    !$acc&     copyout(pvliq, pvice, liqtend, icetend, wvtend, htend, sfliq, sfice) &
    !$acc&     create(fxliq, fxice, psi, fdot, s, sh, d, dh, e, eh, delxh)
    !$acc parallel loop gang vector reduction(max:fail_liq, fail_ice) &
    !$acc&  private(rho, vfall, icice, iciwc, icefac, logiwc, rei, rel, cldab, &
    !$acc&          evapliq, evapice, cldovrl, found, k)
    do i = 1, ncol

      !--------------------------------------------------
      ! COMPUTE GRAVITATIONAL SEDIMENTATION VELOCITIES
      !--------------------------------------------------
      ! NEED TO BE CAREFUL - VELOCITIES SHOULD BE AT THE *LOWER* INTERFACE
      ! (THAT IS, FOR K+1), FLUXES ARE ALSO AT THE LOWER INTERFACE (K+1),
      ! BUT MIXING RATIOS ARE AT THE MIDPOINTS (K)...
      !
      ! NOTE THAT PVEL IS ON INTERFACES, WHEREAS VFALL IS FOR THE CELL
      ! AVERAGES (I.E., MIDPOINTS); ASSUME THE FALL VELOCITY APPLICABLE TO THE
      ! LOWER INTERFACE (K+1) IS THE SAME AS THAT APPLICABLE FOR THE CELL (V(K))

      ! LIQUID:
      ! The fall velocities assume that droplets have a gamma distribution
      ! with effective radii for land and ocean as assessed by Han et al.;
      ! see Lawrence and Crutzen (1998) for a derivation.
      ! ICE:
      ! The fall velocities are based on data from McFarquhar and Heymsfield
      ! or on Stokes terminal velocity for spheres and the effective radius.
      pvliq(i, 1) = 0._kind_phys
      pvice(i, 1) = 0._kind_phys
      !$acc loop seq
      do k = 1, pver
        rho = pmid(i, k)/(rair*t(i, k))

        pvliq(i, k + 1) = 0._kind_phys
        if (cloud(i, k) > 0._kind_phys .and. cldliq(i, k) > 0._kind_phys) then
          ! effective radius of liquid drop (reltab)
          rel = rel_point(tmelt, t(i, k), snowh(i), landm(i), icefrac(i))
          if (rel < 40._kind_phys) then
            vfall = 2._kind_phys/9._kind_phys*rhoh2o*gravit*rel**2/eta*um_to_m  ! um^2 -> m^2
          else
            vfall = v40 + vslope*(rel - r40)    ! linear above 40 microns
          end if
          ! convert the fall speed to pressure units, but do not apply the traditional
          ! negative convention for pvel.
          pvliq(i, k + 1) = vfall*rho*gravit      ! m s-1 to Pa s-1
        end if

        pvice(i, k + 1) = 0._kind_phys
        if (cloud(i, k) > 0._kind_phys .and. cldice(i, k) > 0._kind_phys) then
          if (stokes) then
            ! stokes terminal velocity < 40 microns
            ! effective radius of ice (reitab)
            rei = rei_point(t(i, k))
            if (rei < 40._kind_phys) then
              vfall = 2._kind_phys/9._kind_phys*rhoh2o*gravit*rei**2/eta*um_to_m  ! microns^2 -> m^2
              vfall = vfall*cldsed_ice_stokes_fac
            else
              vfall = v40 + vslope*(rei - r40)      ! linear above 40 microns
            end if

            ! convert the fall speed to pressure units, but do not apply the traditional
            ! negative convention for pvel.
            pvice(i, k + 1) = vfall*rho*gravit      ! m s-1 to Pa s-1
          else
            ! McFarquhar and Heymsfield > icritc

            ! compute the in-cloud ice concentration (kg/kg)
            icice = cldice(i, k)/cloud(i, k)

            ! compute the ice water content in g/m3
            iciwc = icice*rho*1.e3_kind_phys

            ! set the fall velocity (cm/s) to depend on the ice water content in g/m3,
            if (iciwc > lbound) then ! need this because of log10
              logiwc = log10(iciwc)
              !       Median -
              vfall = 128.64_kind_phys + 53.242_kind_phys*logiwc + 5.4795_kind_phys*logiwc**2
              !       Average -
              ! vfall = 122.63 + 44.111*logiwc + 4.2144*logiwc**2
            else
              vfall = 0._kind_phys
            end if

            ! set ice velocity to 1 cm/s if ice mixing ratio < icritc, ramp to value
            ! calculated above at 2*icritc
            if (icice <= icritc) then
              vfall = vice_small
            else if (icice < 2*icritc) then
              icefac = (icice - icritc)/icritc
              vfall = vice_small*(1._kind_phys - icefac) + vfall*icefac
            end if

            ! bound the terminal velocity of ice particles at high concentration
            vfall = min(100.0_kind_phys, vfall)

            ! convert the fall speed to pressure units, but do not apply the traditional
            ! negative convention for pvel.
            pvice(i, k + 1) = vfall &
                              *0.01_kind_phys &      ! cm to meters
                              *rho*gravit            ! m s-1 to Pa s-1
          end if
        end if
      end do

      !--------------------------------------------------
      ! APPLY CLOUD PARTICLE GRAVITATIONAL SEDIMENTATION TO CONDENSATE
      !--------------------------------------------------
      ! Get fluxes at interior points (getflx sets the boundary fluxes to zero)
      call getflx_column(i, pver, pverp, pint, cldliq, pvliq, dtime, fxliq, &
                         psi, fdot, s, sh, d, dh, e, eh, delxh, found)
      if (.not. found) fail_liq = max(fail_liq, i)
      call getflx_column(i, pver, pverp, pint, cldice, pvice, dtime, fxice, &
                         psi, fdot, s, sh, d, dh, e, eh, delxh, found)
      if (.not. found) fail_ice = max(fail_ice, i)

      ! Calculate fluxes at boundaries
      fxliq(i, 1) = 0._kind_phys
      fxice(i, 1) = 0._kind_phys
      ! surface flux by upstream scheme
      fxliq(i, pverp) = cldliq(i, pver)*pvliq(i, pverp)*dtime
      fxice(i, pverp) = cldice(i, pver)*pvice(i, pverp)*dtime

      ! filter out any negative fluxes from the getflx routine
      ! (typical fluxes are of order > 1e-3 when clouds are present)
      !$acc loop seq
      do k = 2, pver
        fxliq(i, k) = max(0._kind_phys, fxliq(i, k))
        fxice(i, k) = max(0._kind_phys, fxice(i, k))
      end do

      ! Limit the flux out of the bottom of each cell to the water content in each phase.
      ! Apply mxsedfac to prevent generating very small negative cloud water/ice
      ! NOTE, REMOVED CLOUD FACTOR FROM AVAILABLE WATER. ALL CLOUD WATER IS IN CLOUDS.
      ! ***Should we include the flux in the top, to allow for thin surface layers?
      ! ***Requires simple treatment of cloud overlap, already included below.
      !$acc loop seq
      do k = 1, pver
        fxliq(i, k + 1) = min(fxliq(i, k + 1), mxsedfac*cldliq(i, k)*pdel(i, k))
        fxice(i, k + 1) = min(fxice(i, k + 1), mxsedfac*cldice(i, k)*pdel(i, k))
        ! fxliq(i,k+1) = min( fxliq(i,k+1), cldliq(i,k) * pdel(i,k) + fxliq(i,k))
        ! fxice(i,k+1) = min( fxice(i,k+1), cldice(i,k) * pdel(i,k) + fxice(i,k))
        ! fxliq(i,k+1) = min( fxliq(i,k+1), cloud(i,k) * cldliq(i,k) * pdel(i,k) )
        ! fxice(i,k+1) = min( fxice(i,k+1), cloud(i,k) * cldice(i,k) * pdel(i,k) )
      end do

      ! Now calculate the tendencies assuming that condensate evaporates when
      ! it falls into environment, and does not when it falls into cloud.
      ! All flux out of the layer comes from the cloudy part.
      ! Assume maximum overlap for stratiform clouds
      !  if cloud above < cloud, all water falls into cloud below
      !  if cloud above > cloud, water split between cloud and environment
      !$acc loop seq
      do k = 1, pver
        ! As in the original code, the cloud above is reset to zero at
        ! every level, so cldovrl uses cldab = 0.
        cldab = 0._kind_phys
        ! cloud overlap cloud factor
        cldovrl = min(cloud(i, k)/(cldab + .0001_kind_phys), 1._kind_phys)
        cldab = cloud(i, k)
        ! evaporation into environment cause moistening and cooling
        evapliq = fxliq(i, k)*(1._kind_phys - cldovrl)/(dtime*pdel(i, k))  ! into env [kg kg-1 s-1]
        evapice = fxice(i, k)*(1._kind_phys - cldovrl)/(dtime*pdel(i, k))  ! into env [kg kg-1 s-1]
        wvtend(i, k) = evapliq + evapice                                   ! evaporation into environment [kg kg-1 s-1]
        htend(i, k) = -latvap*evapliq - (latvap + latice)*evapice          ! evaporation [W kg-1]
        ! net flux into cloud changes cloud liquid/ice (all flux is out of cloud)
        liqtend(i, k) = (fxliq(i, k)*cldovrl - fxliq(i, k + 1))/(dtime*pdel(i, k))
        icetend(i, k) = (fxice(i, k)*cldovrl - fxice(i, k + 1))/(dtime*pdel(i, k))
      end do

      ! convert flux out the bottom to mass units Pa -> kg/m2/s
      sfliq(i) = fxliq(i, pverp)/(dtime*gravit)
      sfice(i) = fxice(i, pverp)/(dtime*gravit)

      ! Convert stratiform_lwe_cloud_ice_surface_flux_due_to_sedimentation from kg m-2 s-1 to precip units m s-1
      sfice(i) = sfice(i)/1000._kind_phys
    end do
    !$acc end data

    if (fail_liq > 0) then
      write(errmsg,'(a,i0)') 'cloud_particle_sedimentation/cfint2: interval was not found for column ', fail_liq
      errflg = 1
    else if (fail_ice > 0) then
      write(errmsg,'(a,i0)') 'cloud_particle_sedimentation/cfint2: interval was not found for column ', fail_ice
      errflg = 1
    end if

  end subroutine cloud_particle_sedimentation_run

  ! Compute fluxes at interior points of column icol. Arrays are whole,
  ! (ncol, ...); psi .. delxh are work arrays. found is false if cfint2 did
  ! not find the interval at some level (flux is then incomplete).
  subroutine getflx_column(icol, pver, pverp, &
    xw, phi, vel, deltat, flux, &
    psi, fdot, s, sh, d, dh, e, eh, delxh, found)
    !$acc routine seq
    !.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
    !....psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
    !....velw1.....velw2.....velw3.....velw4.....velw5.....velw6
    !.........phi1......phi2.......phi3.....phi4.......phi5.......

    integer,         intent(in)    :: icol
    integer,         intent(in)    :: pver
    integer,         intent(in)    :: pverp
    real(kind_phys), intent(in)    :: xw(:,:)      ! (ncol, pverp)
    real(kind_phys), intent(in)    :: phi(:,:)     ! (ncol, pver)
    real(kind_phys), intent(in)    :: vel(:,:)     ! (ncol, pverp)
    real(kind_phys), intent(in)    :: deltat

    real(kind_phys), intent(inout) :: flux(:,:)    ! (ncol, pverp)
    real(kind_phys), intent(inout) :: psi(:,:), fdot(:,:), s(:,:), sh(:,:), &
                                      d(:,:), dh(:,:), e(:,:), eh(:,:), delxh(:,:)
    logical,         intent(out)   :: found

    integer         :: k
    real(kind_phys) :: xxk
    real(kind_phys) :: psistar

    found = .true.

    ! integral of phi
    psi(icol, 1) = 0._kind_phys
    ! fluxes at boundaries
    flux(icol, 1) = 0._kind_phys
    flux(icol, pverp) = 0._kind_phys

    ! integral function
    do k = 2, pverp
      psi(icol, k) = phi(icol, k - 1)*(xw(icol, k) - xw(icol, k - 1)) + psi(icol, k - 1)
    end do

    ! Calculate the derivatives for the interpolating polynomial
    call cfdotmc_column(icol, pver, pverp, xw, psi, fdot, s, sh, d, dh, e, eh, delxh)

    !  NEW WAY
    !    calculate fluxes at interior pts
    do k = 2, pver
      xxk = xw(icol, k) - vel(icol, k)*deltat
      call cfint2_column(icol, pverp, xw, psi, fdot, xxk, psistar, found)
      if (.not. found) return
      flux(icol, k) = (psi(icol, k) - psistar)
    end do
  end subroutine getflx_column

  ! Vertical level interpolation at xin in column icol. found is false if the
  ! interval was not found (psistar is then zero).
  subroutine cfint2_column(icol, pverp, x, f, fdot, xin, psistar, found)
    !$acc routine seq

    integer,         intent(in)  :: icol
    integer,         intent(in)  :: pverp
    real(kind_phys), intent(in)  :: x(:,:)       ! (ncol, pverp)
    real(kind_phys), intent(in)  :: f(:,:)       ! (ncol, pverp)
    real(kind_phys), intent(in)  :: fdot(:,:)    ! (ncol, pverp)
    real(kind_phys), intent(in)  :: xin

    real(kind_phys), intent(out) :: psistar
    logical,         intent(out) :: found

    integer :: k
    integer :: intz
    real(kind_phys) :: dx
    real(kind_phys) :: s
    real(kind_phys) :: c2
    real(kind_phys) :: c3
    real(kind_phys) :: xx
    real(kind_phys) :: psi1, psi2, psi3, psim
    real(kind_phys) :: cfint
    real(kind_phys) :: cfnew
    real(kind_phys) :: xins

    xins = medan(x(icol, 1), xin, x(icol, pverp))
    intz = 0

    ! first find the interval
    do k = 1, pverp - 1
      if ((xins - x(icol, k))*(x(icol, k + 1) - xins) .ge. 0) then
        intz = k
      end if
    end do

    found = intz /= 0
    if (.not. found) then
      psistar = 0._kind_phys
      return
    end if

    ! now interpolate
    k = intz
    dx = (x(icol, k + 1) - x(icol, k))
    s = (f(icol, k + 1) - f(icol, k))/dx
    c2 = (3*s - 2*fdot(icol, k) - fdot(icol, k + 1))/dx
    c3 = (fdot(icol, k) + fdot(icol, k + 1) - 2*s)/dx**2
    xx = (xins - x(icol, k))
    cfint = ((c3*xx + c2)*xx + fdot(icol, k))*xx + f(icol, k)

    ! limit the interpolant
    psi1 = f(icol, k) + (f(icol, k + 1) - f(icol, k))*xx/dx
    if (k .eq. 1) then
      psi2 = f(icol, 1)
    else
      psi2 = f(icol, k) + (f(icol, k) - f(icol, k - 1))*xx/(x(icol, k) - x(icol, k - 1))
    end if
    if (k + 1 .eq. pverp) then
      psi3 = f(icol, pverp)
    else
      psi3 = f(icol, k + 1) - (f(icol, k + 2) - f(icol, k + 1))*(dx - xx)/(x(icol, k + 2) - x(icol, k + 1))
    end if
    psim = medan(psi1, psi2, psi3)
    cfnew = medan(cfint, psi1, psim)
    !if (abs(cfnew - cfint)/(abs(cfnew) + abs(cfint) + 1.e-36_kind_phys) .gt. .03_kind_phys) then
    !    CHANGE THIS BACK LATER!!!
    !    $      .gt..1) then
    !    UNCOMMENT THIS LATER!!!
    !        write(iulog,*) ' cfint2 limiting important ', cfint, cfnew
    !end if
    psistar = cfnew
  end subroutine cfint2_column

  ! Calculate the derivative for the interpolating polynomial in column icol.
  ! Arrays are whole, (ncol, pverp); s .. delxh are work arrays.
  subroutine cfdotmc_column(icol, pver, pverp, x, f, fdot, s, sh, d, dh, e, eh, delxh)
    !$acc routine seq
    ! assumed variable distribution
    !    x1.......x2.......x3.......x4.......x5.......x6    1,pverp points
    !    f1.......f2.......f3.......f4.......f5.......f6    1,pverp points
    !    ...sh1.......sh2......sh3......sh4......sh5....    1,pver points
    !    .........d2.......d3.......d4.......d5.........    2,pver points
    !    .........s2.......s3.......s4.......s5.........    2,pver points
    !    .............dh2......dh3......dh4.............    2,pver-1 points
    !    .............eh2......eh3......eh4.............    2,pver-1 points
    !    ..................e3.......e4..................    3,pver-1 points
    !    .................ppl3......ppl4................    3,pver-1 points
    !    .................ppr3......ppr4................    3,pver-1 points
    !    .................t3........t4..................    3,pver-1 points
    !    ................fdot3.....fdot4................    3,pver-1 points
    integer,         intent(in)    :: icol
    integer,         intent(in)    :: pver
    integer,         intent(in)    :: pverp
    real(kind_phys), intent(in)    :: x(:,:)
    real(kind_phys), intent(in)    :: f(:,:)

    real(kind_phys), intent(inout) :: fdot(:,:)    ! derivative at nodes
    real(kind_phys), intent(inout) :: s(:,:)       ! first divided differences at nodes
    real(kind_phys), intent(inout) :: sh(:,:)      ! first divided differences between nodes
    real(kind_phys), intent(inout) :: d(:,:)       ! second divided differences at nodes
    real(kind_phys), intent(inout) :: dh(:,:)      ! second divided differences between nodes
    real(kind_phys), intent(inout) :: e(:,:)       ! third divided differences at nodes
    real(kind_phys), intent(inout) :: eh(:,:)      ! third divided differences between nodes
    real(kind_phys), intent(inout) :: delxh(:,:)

    integer :: i, k
    real(kind_phys) :: pp                      ! p prime
    real(kind_phys) :: ppl                     ! p prime on left
    real(kind_phys) :: ppr                     ! p prime on right
    real(kind_phys) :: qpl
    real(kind_phys) :: qpr
    real(kind_phys) :: ttt
    real(kind_phys) :: t
    real(kind_phys) :: tmin
    real(kind_phys) :: tmax
    real(kind_phys) :: fd                      ! fdot(i, k) being computed

    ! fdot(i, k) is computed in the scalar fd and stored once: with the
    ! store-then-reload form (fdot(i, k) = ...; fdot(i, k) = f(fdot(i, k))),
    ! nvfortran 25.9's device code dropped the first store and reloaded a stale
    ! fdot(i, k) (GPU replay mismatch, 2026-11-18).

    i = icol

    do k = 1, pver
      ! first divided differences between nodes
      delxh(i, k) = (x(i, k + 1) - x(i, k))
      sh(i, k) = (f(i, k + 1) - f(i, k))/delxh(i, k)

      ! first and second divided differences at nodes
      if (k .ge. 2) then
        d(i, k) = (sh(i, k) - sh(i, k - 1))/(x(i, k + 1) - x(i, k - 1))
        s(i, k) = minmod(sh(i, k), sh(i, k - 1))
      end if
    end do

    ! second and third divided diffs between nodes
    do k = 2, pver - 1
      eh(i, k) = (d(i, k + 1) - d(i, k))/(x(i, k + 2) - x(i, k - 1))
      dh(i, k) = minmod(d(i, k), d(i, k + 1))
    end do

    ! treat the boundaries
    e(i, 2) = eh(i, 2)
    e(i, pver) = eh(i, pver - 1)
    ! outside level
    fd = sh(i, 1) - d(i, 2)*delxh(i, 1) &
         - eh(i, 2)*delxh(i, 1)*(x(i, 1) - x(i, 3))
    fdot(i, 1) = minmod(fd, 3*sh(i, 1))
    fd = sh(i, pver) + d(i, pver)*delxh(i, pver) &
         + eh(i, pver - 1)*delxh(i, pver)*(x(i, pverp) - x(i, pver - 1))
    fdot(i, pverp) = minmod(fd, 3*sh(i, pver))
    ! one in from boundary
    fd = sh(i, 1) + d(i, 2)*delxh(i, 1) - eh(i, 2)*delxh(i, 1)*delxh(i, 2)
    fdot(i, 2) = minmod(fd, 3*s(i, 2))
    fd = sh(i, pver) - d(i, pver)*delxh(i, pver) &
         - eh(i, pver - 1)*delxh(i, pver)*delxh(i, pver - 1)
    fdot(i, pver) = minmod(fd, 3*s(i, pver))

    do k = 3, pver - 1
      e(i, k) = minmod(eh(i, k), eh(i, k - 1))
    end do

    do k = 3, pver - 1
      ! p prime at k-0.5
      ppl = sh(i, k - 1) + dh(i, k - 1)*delxh(i, k - 1)
      ! p prime at k+0.5
      ppr = sh(i, k) - dh(i, k)*delxh(i, k)

      t = minmod(ppl, ppr)

      ! derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
      pp = sh(i, k - 1) + d(i, k)*delxh(i, k - 1)

      ! quartic estimate of fdot
      fd = pp &
           - delxh(i, k - 1)*delxh(i, k) &
           *(eh(i, k - 1)*(x(i, k + 2) - x(i, k)) &
             + eh(i, k)*(x(i, k) - x(i, k - 2)) &
             )/(x(i, k + 2) - x(i, k - 2))

      ! now limit it
      qpl = sh(i, k - 1) &
            + delxh(i, k - 1)*minmod(d(i, k - 1) + e(i, k - 1)*(x(i, k) - x(i, k - 2)), &
                                     d(i, k) - e(i, k)*delxh(i, k))
      qpr = sh(i, k) &
            + delxh(i, k)*minmod(d(i, k) + e(i, k)*delxh(i, k - 1), &
                                 d(i, k + 1) + e(i, k + 1)*(x(i, k) - x(i, k + 2)))

      fd = medan(fd, qpl, qpr)

      ttt = minmod(qpl, qpr)
      tmin = min(0._kind_phys, 3*s(i, k), 1.5_kind_phys*t, ttt)
      tmax = max(0._kind_phys, 3*s(i, k), 1.5_kind_phys*t, ttt)

      fdot(i, k) = fd + minmod(tmin - fd, tmax - fd)
    end do
  end subroutine cfdotmc_column

  pure function minmod(a, b) result(res)
      !$acc routine seq
      real(kind_phys), intent(in) :: a, b
      real(kind_phys) :: res
      res = 0.5_kind_phys*(sign(1._kind_phys, a) + sign(1._kind_phys, b))*min(abs(a), abs(b))
  end function minmod

  pure function medan(a, b, c) result(res)
      !$acc routine seq
      real(kind_phys), intent(in) :: a, b, c
      real(kind_phys) :: res
      res = a + minmod(b - a, c - a)
  end function medan

end module cloud_particle_sedimentation
