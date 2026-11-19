! Interstitial schemes for preparing fields needed for various gravity wave drag schemes
module gravity_wave_drag_interstitials
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  ! public CCPP-compliant subroutines
  public :: gravity_wave_drag_prepare_profiles_run
  public :: gravity_wave_drag_pressure_coordinates_run
  public :: gravity_wave_drag_pressure_coordinates_timestep_final

contains

  ! (originally gw_prof in gw_common)
  ! Compute profiles of background state quantities for multiple
  ! gravity wave drag parameterizations.
  !
  ! The parameterization is assumed to operate only where water vapor
  ! concentrations are negligible in determining the density.
  !
  ! The pressure coordinates are computed per column from pint (the
  ! arithmetic of Coords1D(pint): ifc = pint, mid = 0.5*(ifc(k)+ifc(k+1)),
  ! rdst = 1/(mid(k+1)-mid(k))); schemes that still read the Coords1D object
  ! get it from gravity_wave_drag_pressure_coordinates.
  !
  ! OpenACC: one data region; zero-fills and kvt_gw are pointwise kernels,
  ! the profiles one thread per column with levels in order. Level values
  ! are kept in scalars and stored once.
!> \section arg_table_gravity_wave_drag_prepare_profiles_run Argument Table
!! \htmlinclude gravity_wave_drag_prepare_profiles_run.html
  subroutine gravity_wave_drag_prepare_profiles_run( &
    ncol, pver, &
    cpair, rair, gravit, &
    pint, &
    t, &
    do_molec_diff, nbot_molec, &
    cpairv, kvt, &
    rhoi, &
    nm, ni, &
    egwdffi, &
    tend_q, tend_u, tend_v, tend_s, &
    kvt_gw, &
    flx_heat, &
    scheme_name, &
    errmsg, errflg)

    ! Input arguments
    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver
    real(kind_phys), intent(in)  :: cpair      ! Specific heat of dry air, constant pressure [J K-1 kg-1]
    real(kind_phys), intent(in)  :: rair       ! Dry air gas constant [J K-1 kg-1]
    real(kind_phys), intent(in)  :: gravit     ! standard_gravitational_acceleration [m s-2]
    real(kind_phys), intent(in)  :: pint(:, :) ! Interface pressures [Pa]
    real(kind_phys), intent(in)  :: t(:, :)    ! Midpoint temperatures [K]
    logical,         intent(in)  :: do_molec_diff ! do_molecular_diffusion [flag]
    integer,         intent(in)  :: nbot_molec ! index_of_pressure_at_bottom_of_molecular_diffusion [index]
    real(kind_phys), intent(in)  :: cpairv(:,:)! composition_dependent_specific_heat_of_dry_air_at_constant_pressure [J kg-1 K-1]
    real(kind_phys), intent(in)  :: kvt(:,:)   ! molecular_kinematic_temperature_conductivity_at_interfaces

    ! Output arguments
    real(kind_phys), intent(out) :: rhoi(:, :) ! Density at layer interfaces [kg m-3]
    real(kind_phys), intent(out) :: nm(:, :)   ! Midpoint Brunt-Vaisalla frequencies [s-1]
    real(kind_phys), intent(out) :: ni(:, :)   ! Interface Brunt-Vaisalla frequencies [s-1]
    real(kind_phys), intent(out) :: egwdffi(:, :)   ! Effective diffusion coefficient from gravity waves [m2 s-1]

    ! Initialize tendency accumulators to make framework happy about intent(inout) in the individual gw schemes.
    real(kind_phys), intent(out) :: tend_q(:, :, :) ! Constituent tendencies [kg kg-1 s-1]
    real(kind_phys), intent(out) :: tend_u(:, :)    ! Zonal wind tendency [m s-2]
    real(kind_phys), intent(out) :: tend_v(:, :)    ! Meridional wind tendency [m s-2]
    real(kind_phys), intent(out) :: tend_s(:, :)    ! Dry static energy tendency [J kg-1 s-1]

    real(kind_phys), intent(out) :: kvt_gw(:,:)! scaled_molecular_kinematic_temperature_conductivity_at_interfaces_for_gravity_wave_drag [m2 s-1]

    ! for energy checker:
    real(kind_phys), intent(out) :: flx_heat(:) ! net_sensible_heat_flux_through_top_and_bottom_of_atmosphere_column [W m-2]
    character(len=*), intent(out):: scheme_name

    character(len=512), intent(out) :: errmsg
    integer, intent(out)            :: errflg

    ! Local variables
    integer :: i, k, m

    real(kind_phys) :: dtdp                    ! dt/dp
    real(kind_phys) :: n2                      ! Brunt-Vaisalla frequency squared [s-2]
    real(kind_phys) :: ti                      ! Interface temperature [K]
    real(kind_phys) :: rho                     ! Interface density [kg m-3]
    real(kind_phys) :: ni_k, ni_km1            ! Interface B-V frequency at k, k-1 [s-1]
    real(kind_phys) :: mid_k, mid_km1          ! Midpoint pressure at k, k-1 [Pa]
    real(kind_phys) :: rdst                    ! 1/(mid(k)-mid(k-1)) [Pa-1]

    ! Minimum value of Brunt-Vaisalla frequency squared.
    real(kind_phys), parameter :: n2min = 5.e-5_kind_phys

    errmsg = ''
    errflg = 0
    scheme_name = 'gravity_wave_drag'

    !$acc data copyin(pint, t, cpairv, kvt) &
    !$acc      copyout(rhoi, nm, ni, egwdffi, tend_q, tend_u, tend_v, tend_s, kvt_gw, flx_heat)

    !------------------------------------------------------------------------
    ! Initialization of effective diffusion coefficients,
    ! accummulated tendencies from all gravity wave schemes,
    ! and quantities for energy checker.
    !------------------------------------------------------------------------
    !$acc parallel loop gang vector collapse(2)
    do k = 1, size(egwdffi, 2)
      do i = 1, size(egwdffi, 1)
        egwdffi(i, k) = 0._kind_phys
      end do
    end do
    !$acc parallel loop gang vector
    do i = 1, size(flx_heat)
      flx_heat(i) = 0._kind_phys
    end do
    !$acc parallel loop gang vector collapse(3)
    do m = 1, size(tend_q, 3)
      do k = 1, size(tend_q, 2)
        do i = 1, size(tend_q, 1)
          tend_q(i, k, m) = 0._kind_phys
        end do
      end do
    end do
    !$acc parallel loop gang vector collapse(2)
    do k = 1, size(tend_u, 2)
      do i = 1, size(tend_u, 1)
        tend_u(i, k) = 0._kind_phys
        tend_v(i, k) = 0._kind_phys
        tend_s(i, k) = 0._kind_phys
      end do
    end do

    !------------------------------------------------------------------------
    ! Calculate local molecular diffusivity
    !------------------------------------------------------------------------
    if (do_molec_diff) then
      ! Use linear extrapolation of cpairv to top interface;
      ! interpolate cpairv to other interfaces down to nbot_molec.
      !$acc parallel loop gang vector collapse(2)
      do k = 1, size(kvt_gw, 2)
        do i = 1, size(kvt_gw, 1)
          if (k == 1) then
            kvt_gw(i, k) = kvt(i, k)/ &
                           (1.5_kind_phys*cpairv(i, 1) - &
                            0.5_kind_phys*cpairv(i, 2))
          else if (k <= nbot_molec) then
            kvt_gw(i, k) = kvt(i, k)/ &
                           (cpairv(i, k + 1) + cpairv(i, k))*2._kind_phys
          else
            kvt_gw(i, k) = kvt(i, k)
          end if
        end do
      end do
    else
      !$acc parallel loop gang vector collapse(2)
      do k = 1, size(kvt_gw, 2)
        do i = 1, size(kvt_gw, 1)
          kvt_gw(i, k) = 0._kind_phys
        end do
      end do
    endif

    !------------------------------------------------------------------------
    ! Determine the interface densities and Brunt-Vaisala frequencies,
    ! then the midpoint Brunt-Vaisala frequencies (midpoint_interp of ni).
    !------------------------------------------------------------------------
    !$acc parallel loop gang vector &
    !$acc   private(k, ti, rho, dtdp, n2, ni_k, ni_km1, mid_k, mid_km1, rdst)
    do i = 1, ncol
      ! The top interface values are calculated assuming an isothermal
      ! atmosphere above the top level.
      k = 1
      ti = t(i, k)
      rhoi(i, k) = pint(i, k)/(rair*ti)
      ni_k = sqrt(gravit*gravit/(cpair*ti))
      ni(i, k) = ni_k
      mid_k = 0.5_kind_phys*(pint(i, k) + pint(i, k + 1))

      ! Interior points use centered differences.
      !$acc loop seq
      do k = 2, pver
        ni_km1 = ni_k
        mid_km1 = mid_k
        mid_k = 0.5_kind_phys*(pint(i, k) + pint(i, k + 1))
        rdst = 1._kind_phys/(mid_k - mid_km1)
        ti = 0.5_kind_phys*(t(i, k - 1) + t(i, k))
        rho = pint(i, k)/(rair*ti)
        rhoi(i, k) = rho
        dtdp = (t(i, k) - t(i, k - 1))*rdst
        n2 = gravit*gravit/ti*(1._kind_phys/cpair - rho*dtdp)
        ni_k = sqrt(max(n2min, n2))
        ni(i, k) = ni_k
        nm(i, k - 1) = 0.5_kind_phys*(ni_km1 + ni_k)
      end do

      ! Bottom interface uses bottom level temperature, density; next interface
      ! B-V frequency.
      k = pver + 1
      ti = t(i, k - 1)
      rhoi(i, k) = pint(i, k)/(rair*ti)
      ni(i, k) = ni_k
      nm(i, pver) = 0.5_kind_phys*(ni_k + ni_k)
    end do

    !$acc end data

  end subroutine gravity_wave_drag_prepare_profiles_run

  ! Pressure coordinates (coords1d) of the gravity wave schemes that still
  ! read the object (the CAM7 sources). prepare_profiles and the orographic
  ! source compute the coordinates from pint; this host object is built only
  ! in suites that list this scheme.
!> \section arg_table_gravity_wave_drag_pressure_coordinates_run Argument Table
!! \htmlinclude gravity_wave_drag_pressure_coordinates_run.html
  subroutine gravity_wave_drag_pressure_coordinates_run( &
    ncol, pver, &
    pint, &
    p, &
    errmsg, errflg)

    use coords_1d, only: Coords1D

    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    real(kind_phys),    intent(in)  :: pint(:, :) ! Interface pressures [Pa]
    type(coords1d),     intent(out) :: p          ! Vertical coordinate for gravity waves [Coords1D, Pa]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    p = Coords1D(pint(:ncol, :pver+1))

  end subroutine gravity_wave_drag_pressure_coordinates_run

!> \section arg_table_gravity_wave_drag_pressure_coordinates_timestep_final Argument Table
!! \htmlinclude gravity_wave_drag_pressure_coordinates_timestep_final.html
  subroutine gravity_wave_drag_pressure_coordinates_timestep_final(p, errmsg, errflg)
    use coords_1d, only: Coords1D

    type(coords1d),     intent(inout)  :: p          ! Vertical coordinate for gravity waves [Coords1D, Pa]
    character(len=512), intent(out)    :: errmsg
    integer, intent(out)               :: errflg

    errmsg = ''
    errflg = 0

    call p%finalize()

  end subroutine gravity_wave_drag_pressure_coordinates_timestep_final

end module gravity_wave_drag_interstitials
