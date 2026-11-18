! Interstitial scheme to prepare inputs for vertical diffusion solver
! common to the Holtslag-Boville (HB) and Bretherton-Park (UW/diag_TKE)
! PBL schemes.
!
! These include:
!   - potential temperature (th)
!   - the p vertical coordinate (coords1d), including destroying this object at
!     end of timestep, in its own scheme (vertical_diffusion_pressure_coordinates),
!     listed only in suites whose schemes read p (e.g. bretherton_park_diff).
module vertical_diffusion_interstitials
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  ! CCPP-compliant public interfaces
  public :: vertical_diffusion_prepare_inputs_run

  public :: vertical_diffusion_pressure_coordinates_run
  public :: vertical_diffusion_pressure_coordinates_timestep_final

  public :: compute_kinematic_fluxes_and_obklen_run

  public :: vertical_diffusion_set_total_surface_stress_run

contains

  ! Interstitial to compute potential temperature.
!> \section arg_table_vertical_diffusion_prepare_inputs_run Argument Table
!! \htmlinclude vertical_diffusion_prepare_inputs_run.html
  subroutine vertical_diffusion_prepare_inputs_run( &
    ncol, pver, &
    t, exner, &
    th, &
    errmsg, errflg)

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    real(kind_phys),    intent(in)  :: t(:,:)    ! temperature [K]
    real(kind_phys),    intent(in)  :: exner(:,:)! exner [1]
    real(kind_phys),    intent(out) :: th(:,:)   ! Potential temperature [K]
    character(len=512), intent(out) :: errmsg    ! Error message
    integer,            intent(out) :: errflg    ! Error flag

    integer :: i, k

    errmsg = ''
    errflg = 0

    ! Calculate potential temperature [K]
    !$acc parallel loop collapse(2) copyin(t, exner) copyout(th)
    do k = 1, pver
      do i = 1, ncol
        th(i,k) = t(i,k) * exner(i,k)
      end do
    end do

  end subroutine vertical_diffusion_prepare_inputs_run

  ! Pressure coordinates (coords1d) of the vertical diffusion solver, for schemes
  ! that use the object (e.g. bretherton_park_diff). The vertical diffusion
  ! schemes themselves compute the coordinates from pint; this host object is
  ! built only in suites that list this scheme.
!> \section arg_table_vertical_diffusion_pressure_coordinates_run Argument Table
!! \htmlinclude vertical_diffusion_pressure_coordinates_run.html
  subroutine vertical_diffusion_pressure_coordinates_run( &
    ncol, pverp, &
    pint, &
    p, &
    errmsg, errflg)

    use coords_1d,  only: coords1d

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pverp
    real(kind_phys),    intent(in)  :: pint(:,:) ! Air pressure at interfaces [Pa]
    type(coords1d),     intent(out) :: p         ! Vertical moist pressure coordinates for vertical diffusion [Pa]
    character(len=512), intent(out) :: errmsg    ! Error message
    integer,            intent(out) :: errflg    ! Error flag

    errmsg = ''
    errflg = 0

    p = coords1d(pint(:ncol,:pverp))

  end subroutine vertical_diffusion_pressure_coordinates_run

!> \section arg_table_vertical_diffusion_pressure_coordinates_timestep_final Argument Table
!! \htmlinclude vertical_diffusion_pressure_coordinates_timestep_final.html
  subroutine vertical_diffusion_pressure_coordinates_timestep_final(p, errmsg, errflg)
    use coords_1d,  only: coords1d

    type(coords1d),     intent(inout) :: p       ! Vertical moist pressure coordinates for vertical diffusion [Pa]
    character(len=512), intent(out)   :: errmsg  ! Error message
    integer,            intent(out)   :: errflg  ! Error flag

    errmsg = ''
    errflg = 0

    call p%finalize()
  end subroutine vertical_diffusion_pressure_coordinates_timestep_final

  ! Compute kinematic fluxes (heat, water vapor, buoyancy), and obukhov length
!> \section arg_table_compute_kinematic_fluxes_and_obklen_run Argument Table
!! \htmlinclude compute_kinematic_fluxes_and_obklen_run.html
  subroutine compute_kinematic_fluxes_and_obklen_run( &
    ncol, pver, pcnst, &
    zvir, cpair, gravit, karman, &
    shf_from_coupler, cflx_from_coupler, &
    q_wv, &
    th, &
    rrho, &
    ustar, &
    ! below output:
    khfs, kqfs, kbfs, &
    obklen, &
    errmsg, errflg)

    ! dependency to get constituent index
    use ccpp_scheme_utils,         only: ccpp_constituent_index

    use atmos_phys_pbl_utils,      only: calc_virtual_temperature, calc_obukhov_length, &
                                         calc_kinematic_heat_flux, calc_kinematic_water_vapor_flux, &
                                         calc_kinematic_buoyancy_flux

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    integer,            intent(in)  :: pcnst
    real(kind_phys),    intent(in)  :: zvir
    real(kind_phys),    intent(in)  :: cpair
    real(kind_phys),    intent(in)  :: gravit
    real(kind_phys),    intent(in)  :: karman
    real(kind_phys),    intent(in)  :: shf_from_coupler(:)     ! Surface upward sensible heat flux from coupler [W m-2]
    real(kind_phys),    intent(in)  :: cflx_from_coupler(:,:)  ! Surface upward constituent fluxes from coupler [kg m-2 s-1]
    real(kind_phys),    intent(in)  :: q_wv(:,:)               ! specific humidity [kg kg-1]

    real(kind_phys),    intent(in)  :: th(:,:)                 ! Potential temperature [K]

    ! rrho and ustar from PBL scheme.
    ! for HB, it is computed in hb_pbl_independent_coefficients_run;
    ! for UW, it is iteratively updated in bretherton_park_diff_run and returned as a final result.
    real(kind_phys),    intent(in)  :: rrho(:)                 ! 1 / bottom level density [m3 kg-1]
    real(kind_phys),    intent(in)  :: ustar(:)                ! surface friction velocity [m s-1]

    ! Output arguments:
    real(kind_phys),    intent(out) :: khfs(:)                 ! kinematic surface heat flux [K m s-1]
    real(kind_phys),    intent(out) :: kqfs(:)                 ! kinematic surface water vapor flux [kg kg-1 m s-1]
    real(kind_phys),    intent(out) :: kbfs(:)                 ! surface kinematic buoyancy flux [m2 s-3]
    real(kind_phys),    intent(out) :: obklen(:)               ! Obukhov length [m]

    character(len=512), intent(out) :: errmsg                  ! Error message
    integer,            intent(out) :: errflg                  ! Error flag

    ! Local variables
    integer :: const_wv_idx                                    ! Water vapor constituent index
    real(kind_phys) :: thvs                                    ! Virtual potential temperature at surface [K]
    integer :: i

    ! Check constituents list and locate water vapor index
    ! (not assumed to be 1)
    call ccpp_constituent_index('water_vapor_mixing_ratio_wrt_moist_air_and_condensed_water', &
         const_wv_idx, errflg, errmsg)
    if (errflg /= 0) return

    !$acc parallel loop private(thvs) &
    !$acc   copyin(th, q_wv, shf_from_coupler, cflx_from_coupler, rrho, ustar) &
    !$acc   copyout(khfs, kqfs, kbfs, obklen)
    do i = 1, ncol
      thvs      = calc_virtual_temperature(th(i,pver), q_wv(i,pver), zvir)
      khfs(i)   = calc_kinematic_heat_flux(shf_from_coupler(i), rrho(i), cpair)
      kqfs(i)   = calc_kinematic_water_vapor_flux(cflx_from_coupler(i,const_wv_idx), rrho(i))
      kbfs(i)   = calc_kinematic_buoyancy_flux(khfs(i), zvir, th(i,pver), kqfs(i))
      obklen(i) = calc_obukhov_length(thvs, ustar(i), gravit, karman, kbfs(i))
    end do

  end subroutine compute_kinematic_fluxes_and_obklen_run

  ! Set total surface stresses for input
  ! 1) into the HB PBL scheme, and
  ! 2) surface zonal/meridional momentum flux diagnostic output (all schemes)
  !
  ! NOTE: Temporarily, TMS and Beljaars are "hard-coupled" into this subroutine because
  ! the current focus is CAM4 and TMS/Beljaars drag have not been CCPPized.
  ! In the future, after TMS/Beljaars is implemented, this subroutine should only initialize
  ! surface stresses from the coupler, then TMS and Beljaars can work on adding the respective
  ! surface stresses to the tautotx/tautoty used by the PBL scheme.
!> \section arg_table_vertical_diffusion_set_total_surface_stress_run Argument Table
!! \htmlinclude vertical_diffusion_set_total_surface_stress_run.html
  subroutine vertical_diffusion_set_total_surface_stress_run( &
    ncol, &
    wsx_from_coupler, wsy_from_coupler, &
    tautmsx, tautmsy, &
    taubljx, taubljy, &
    tautotx, tautoty, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol
    real(kind_phys),    intent(in)  :: wsx_from_coupler(:)      ! Surface eastward wind stress from coupler [Pa]
    real(kind_phys),    intent(in)  :: wsy_from_coupler(:)      ! Surface northward wind stress from coupler [Pa]
    real(kind_phys),    intent(in)  :: tautmsx(:)               ! Eastward turbulent mountain surface stress [Pa]
    real(kind_phys),    intent(in)  :: tautmsy(:)               ! Northward turbulent mountain surface stress [Pa]
    real(kind_phys),    intent(in)  :: taubljx(:)               ! Eastward Beljaars surface stress [Pa]
    real(kind_phys),    intent(in)  :: taubljy(:)               ! Northward Beljaars surface stress [Pa]

    ! Output arguments
    real(kind_phys),    intent(out) :: tautotx(:)               ! Eastward total stress at surface for boundary layer scheme [Pa]
    real(kind_phys),    intent(out) :: tautoty(:)               ! Northward total stress at surface for boundary layer scheme [Pa]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: i

    errmsg = ''
    errflg = 0

    !$acc parallel loop copyin(wsx_from_coupler, wsy_from_coupler, tautmsx, tautmsy, &
    !$acc   taubljx, taubljy) copyout(tautotx, tautoty)
    do i = 1, ncol
      ! Initialize total surface stresses from coupler
      ! These are used for HB diffusion scheme and later PBL diagnostics but
      ! NOT for the vertical diffusion solver, which uses surface stresses from the coupler only,
      ! or just zero (in the case of CLUBB)
      tautotx(i) = wsx_from_coupler(i)
      tautoty(i) = wsy_from_coupler(i)

      ! Add turbulent mountain stress to total surface stress
      tautotx(i) = tautotx(i) + tautmsx(i)
      tautoty(i) = tautoty(i) + tautmsy(i)

      ! Add Beljaars integrated drag to total surface stress
      tautotx(i) = tautotx(i) + taubljx(i)
      tautoty(i) = tautoty(i) + taubljy(i)
    end do

  end subroutine vertical_diffusion_set_total_surface_stress_run

end module vertical_diffusion_interstitials
