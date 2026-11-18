! Vertical diffusion module.
! Solves vertical diffusion equations using a tri-diagonal solver
! The module will also apply countergradient fluxes
! Includes interstitials for applying vertical diffusion tendencies and
! other intermediate calculations formerly in vertical_diffusion_tend.
!
! This is a reduced functionality version;
!  molecular diffusion is unsupported until WACCM-X is ported.
!
! Original authors: B. Boville and others, 1991-2004
! Modularized: J. McCaa, September 2004.
! Updates: Sungsu Park, August 2006 - January 2010
! Reduced/CCPPized: Haipeng Lin, May 2025
module diffusion_solver
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! CCPP-compliant public interfaces

  public :: vertical_diffusion_set_temperature_at_toa_default_run
  public :: vertical_diffusion_interpolate_to_interfaces_run   ! Interpolate t, rho, rair to interfaces
  public :: implicit_surface_stress_add_drag_coefficient_run   ! Add implicit turbulent surface stress (if do_iss)
  public :: vertical_diffusion_wind_damping_rate_run           ! Compute wind damping rate from drag coeff
  ! Beljaars will add wind damping rates at this point

  public :: vertical_diffusion_diffuse_horizontal_momentum_timestep_init
  public :: vertical_diffusion_diffuse_horizontal_momentum_run ! Vertical diffusion solver (horizontal momentum)

  public :: vertical_diffusion_set_dry_static_energy_at_toa_zero_run
  public :: vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run
  public :: vertical_diffusion_diffuse_dry_static_energy_run   ! Vertical diffusion solver (dry static energy)

  public :: vertical_diffusion_diffuse_tracers_init
  public :: vertical_diffusion_diffuse_tracers_run             ! Vertical diffusion solver (tracers)

  public :: vertical_diffusion_tendencies_run                  ! Using output from solver compute phys. tend

  ! Parameters for implicit surface stress treatment
  real(kind_phys), parameter :: wsmin = 1._kind_phys         ! Minimum sfc wind speed for estimating frictional
                                                             ! transfer velocity ksrf. [m s-1]
  real(kind_phys), parameter :: ksrfmin = 1.e-4_kind_phys    ! Minimum surface drag coefficient [kg s-1 m-2]
  real(kind_phys), parameter :: timeres = 7200._kind_phys    ! Relaxation time scale of residual stress (>= dt) [s]

contains

  ! Set (default) temperature used at top of atmosphere for interpolation.
  ! In the future, this can use the upper boundary condition provided for temperature,
  ! or a different extrapolation in WACCM-X mode.
!> \section arg_table_vertical_diffusion_set_temperature_at_toa_default_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_set_temperature_at_toa_default_run.html
  pure subroutine vertical_diffusion_set_temperature_at_toa_default_run( &
    ncol, pver, &
    t, &
    t_toai, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    real(kind_phys),    intent(in)  :: t(:, :)     ! Temperature [K]

    ! Output arguments
    real(kind_phys),    intent(out) :: t_toai(:)   ! Temperature at interface above TOA for vertical diffusion [K]
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: i

    errmsg = ''
    errflg = 0

    ! Default
    !$acc parallel loop copyin(t) copyout(t_toai)
    do i = 1, ncol
      t_toai(i) = t(i, 1)  ! Use temperature at top model level [K]
    end do

    ! Eventually, if ubc_fixed_temp is available, use ubc_t.

    ! Eventually, in WACCM-X, use extrapolated version:
    ! if (waccmx_mode) then
    !    ! For WACCM-X, set ubc temperature to extrapolate from next two lower interface level temperatures
    !    ! the original formulation is:
    !    ! t_toai(:ncol) = 1.5_r8*tint(:ncol,2)-.5_r8*tint(:ncol,3)
    !    ! this appears to be:
    !    !               = tint(:ncol,2) + 0.5_r8*(tint(:ncol,2) - tint(:ncol,3))
    !    ! assuming that the extrapolated gradient is 1/2 of the temperature gradient between the lower interfaces
    !    ! because the interpolation will be done later in the CCPP-ized scheme, formulate this in terms of
    !    ! the temperature (at midpoints):
    !    t_toai(:ncol) = 1.5_r8*(state%t(:ncol,2)+state%t(:ncol,1))/2._r8-.5_r8*(state%t(:ncol,3)+state%t(:ncol,2))/2._r8

  end subroutine vertical_diffusion_set_temperature_at_toa_default_run

  ! Interpolates temperature, air density (moist and dry),
  ! and sets gas constant (not constituent dependent in non-WACCM-X mode)
  ! at interfaces for vertical diffusion
!> \section arg_table_vertical_diffusion_interpolate_to_interfaces_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_interpolate_to_interfaces_run.html
  pure subroutine vertical_diffusion_interpolate_to_interfaces_run( &
    ncol, pver, pverp, &
    gravit, &
    rair, rairv, &
    flag_for_constituent_dependent_gas_constant, &
    t, &
    t_toai, &
    pint, pintdry, &
    ! below output
    ti, rairi, rhoi, rhoi_dry, &
    dpidz_sq, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)  :: ncol       ! Number of columns
    integer,            intent(in)  :: pver       ! Number of vertical layers
    integer,            intent(in)  :: pverp      ! Number of vertical interfaces (pver+1)
    real(kind_phys),    intent(in)  :: gravit
    real(kind_phys),    intent(in)  :: rair       ! Gas constant for dry air [J kg-1 K-1]
    real(kind_phys),    intent(in)  :: rairv(:,:) ! composition_dependent_gas_constant_of_dry_air [J kg-1 K-1]
    logical,            intent(in)  :: flag_for_constituent_dependent_gas_constant
    real(kind_phys),    intent(in)  :: t(:,:)     ! Air temperature at midpoints [K]
    real(kind_phys),    intent(in)  :: t_toai(:)  ! Air temperature to use at interface above TOA [K]
    real(kind_phys),    intent(in)  :: pint(:,:)  ! Air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)  :: pintdry(:,:) ! Dry air pressure at interfaces [Pa]

    ! Output variables
    real(kind_phys),    intent(out) :: ti(:,:)    ! Air temperature at interfaces [K]
    real(kind_phys),    intent(out) :: rairi(:,:) ! Gas constant for dry air at interfaces [J kg-1 K-1]
    real(kind_phys),    intent(out) :: rhoi(:,:)  ! Air density at interfaces [kg m-3]
    real(kind_phys),    intent(out) :: rhoi_dry(:,:) ! Dry air density at interfaces [kg m-3]
    real(kind_phys),    intent(out) :: dpidz_sq(:,:) ! Square of derivative of pressure with height (moist) [kg2 m-4 s-4]
    character(len=512), intent(out) :: errmsg     ! Error message
    integer,            intent(out) :: errflg     ! Error flag

    ! Local variables
    integer :: i, k

    errmsg = ''
    errflg = 0

    !$acc data copyin(rairv, t, t_toai, pint, pintdry) &
    !$acc      copyout(ti, rairi, rhoi, rhoi_dry, dpidz_sq)

    ! Interpolate temperature to interfaces
    !$acc parallel loop
    do i = 1, ncol
       ti(i,1) = t_toai(i)
       ! Assume temperature at surface interface is the same as the midpoint value in the lowest layer:
       ti(i,pverp) = t(i,pver)
    end do
    !$acc parallel loop collapse(2)
    do k = 2, pver
       do i = 1, ncol
          ti(i,k) = 0.5_kind_phys * (t(i,k) + t(i,k-1))
       end do
    end do

    ! Set gas constant at interfaces
    if(flag_for_constituent_dependent_gas_constant) then
      !$acc parallel loop
      do i = 1, ncol
        rairi(i,1) = rairv(i,1)
        ! Assume gas constant at surface interface is the same as the midpoint value in the lowest layer:
        rairi(i,pverp) = rairv(i,pver)
      end do
      !$acc parallel loop collapse(2)
      do k = 2, pver
        do i = 1, ncol
          rairi(i,k) = 0.5_kind_phys * (rairv(i,k)+rairv(i,k-1))
        end do
      end do
    else
      !$acc parallel loop collapse(2)
      do k = 1, pverp
        do i = 1, ncol
          rairi(i,k) = rair
        end do
      end do
    endif

    ! Compute rho and rho_dry at interfaces, and the square of the derivative of
    ! pressure with height at interfaces.
    ! Note that this equation assumes hydrostatic balance, i.e. the *derivative* dp/dz is g*rho
    !$acc parallel loop collapse(2)
    do k = 1, pverp
       do i = 1, ncol
          rhoi(i,k)  = pint(i,k) / (rairi(i,k)*ti(i,k))
          rhoi_dry(i,k)  = pintdry(i,k) / (rairi(i,k)*ti(i,k))
          dpidz_sq(i,k) = gravit*rhoi(i,k)
          dpidz_sq(i,k) = dpidz_sq(i,k)*dpidz_sq(i,k)
       end do
    end do

    ! FIXME: The following two lines are kept in only to preserve answers;
    !        they really should be taken out completely.
    !$acc parallel loop
    do i = 1, ncol
       dpidz_sq(i, 1) = gravit*(pint(i, 1)/(rairv(i, 1)*t(i, 1)))
       dpidz_sq(i, 1) = dpidz_sq(i, 1)*dpidz_sq(i, 1)
    end do

    !$acc end data

  end subroutine vertical_diffusion_interpolate_to_interfaces_run

  ! Add turbulent surface stress to total surface drag coefficient, if implicit
!> \section arg_table_implicit_surface_stress_add_drag_coefficient_run Argument Table
!! \htmlinclude arg_table_implicit_surface_stress_add_drag_coefficient_run.html
  pure subroutine implicit_surface_stress_add_drag_coefficient_run( &
             ncol, pver, &
             do_iss, &
             taux, tauy, &
             u0, v0, &
             ! below output
             ksrf, &
             errmsg, errflg)

    ! Input arguments
    integer,         intent(in)       :: ncol
    integer,         intent(in)       :: pver
    logical,         intent(in)       :: do_iss             ! Use implicit turbulent surface stress computation
    real(kind_phys), intent(in)       :: taux(:)            ! Surface zonal stress [N m-2]
                                                            ! Input u-momentum per unit time per unit area into the atmosphere.
    real(kind_phys), intent(in)       :: tauy(:)            ! Surface meridional stress [N m-2]
                                                            ! Input v-momentum per unit time per unit area into the atmosphere.

    real(kind_phys), intent(in)       :: u0(:,:)            ! Input u-wind [m s-1]
    real(kind_phys), intent(in)       :: v0(:,:)            ! Input v-wind [m s-1]

    ! Input/output arguments
    real(kind_phys),    intent(inout) :: ksrf(:)            ! total surface drag coefficient [kg m-2 s-1]

    ! Output arguments
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    ! Local variables
    integer :: i
    real(kind_phys) :: ws                        ! Lowest-level wind speed [m s-1]
    real(kind_phys) :: ksrfturb                  ! Surface drag coefficient of 'normal' stress. > 0.
                                                 ! Virtual mass input per unit time per unit area [kg m-2 s-1]
    real(kind_phys) :: tau                       ! Turbulent surface stress (not including mountain stress)

    errmsg = ''
    errflg = 0

    ! Do 'Implicit Surface Stress' treatment for numerical stability
    ! in the lowest model layer.
    if (do_iss) then
      ! Add surface drag coefficient for implicit diffusion
      !$acc parallel loop private(ws, tau, ksrfturb) copyin(taux, tauy, u0, v0) copy(ksrf)
      do i = 1, ncol
        ws = max(sqrt(u0(i, pver)**2 + v0(i, pver)**2), wsmin)
        tau = sqrt(taux(i)**2 + tauy(i)**2)
        ksrfturb = max(tau/ws, ksrfmin)
        ksrf(i) = ksrf(i) + ksrfturb  ! Do all surface stress implicitly
      end do
    endif

  end subroutine implicit_surface_stress_add_drag_coefficient_run

  ! Computes rate at which wind is exponentially damped by surface stress
  ! based on total surface stress (after all surface stress is added)
!> \section arg_table_vertical_diffusion_wind_damping_rate_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_wind_damping_rate_run.html
  pure subroutine vertical_diffusion_wind_damping_rate_run( &
    ncol, pver, &
    gravit, &
    pint, &
    ksrf, &
    tau_damp_rate, &
    errmsg, errflg)

    ! Input arguments
    integer,            intent(in)    :: ncol
    integer,            intent(in)    :: pver
    real(kind_phys),    intent(in)    :: gravit
    real(kind_phys),    intent(in)    :: pint(:,:)          ! Air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)    :: ksrf(:)            ! total surface drag coefficient [kg m-2 s-1]

    ! Output arguments
    real(kind_phys),    intent(out)   :: tau_damp_rate(:,:) ! Rate at which external (surface) stress damps wind speeds [s-1]
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    integer :: i, k
    real(kind_phys) :: rdel_bot  ! 1/(pressure thickness) of the bottom layer, as Coords1D's rdel

    errmsg = ''
    errflg = 0

    !$acc data copyin(ksrf, pint) copyout(tau_damp_rate)

    ! In most layers, no damping at all.
    !$acc parallel loop collapse(2)
    do k = 1, pver
      do i = 1, ncol
        tau_damp_rate(i,k) = 0._kind_phys
      end do
    end do

    ! Compute damp rate based on total surface drag coeff
    !
    ! Physical interpretation:
    ! ksrf is stress per unit wind speed.
    ! p%del / gravit is approximately the mass in the layer per unit of
    ! surface area.
    ! Therefore, gravit*ksrf/p%del is the acceleration of wind per unit
    ! wind speed, i.e. the rate at which wind is exponentially damped by
    ! surface stress.
    !$acc parallel loop private(rdel_bot)
    do i = 1, ncol
      rdel_bot = 1._kind_phys/(pint(i,pver+1) - pint(i,pver))
      tau_damp_rate(i, pver) = -gravit*ksrf(i)*rdel_bot
    end do

    !$acc end data

  end subroutine vertical_diffusion_wind_damping_rate_run

  ! Zero surface stresses at timestep initialization
  ! Any surface stresses, e.g., Turbulent Mountain Stress
  ! will be added to this later as they are computed.
!> \section arg_table_vertical_diffusion_diffuse_horizontal_momentum_timestep_init Argument Table
!! \htmlinclude arg_table_vertical_diffusion_diffuse_horizontal_momentum_timestep_init.html
  pure subroutine vertical_diffusion_diffuse_horizontal_momentum_timestep_init( &
    ncol, pver, &
    ksrf, &
    tau_damp_rate, &
    errmsg, errflg)

    integer,         intent(in)       :: ncol
    integer,         intent(in)       :: pver

    ! Output (initial zeroed) arguments
    real(kind_phys),    intent(out)   :: ksrf(:)            ! total surface drag coefficient [kg m-2 s-1]
    real(kind_phys),    intent(out)   :: tau_damp_rate(:,:) ! Rate at which external (surface) stress damps wind speeds [s-1]

    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    errmsg = ''
    errflg = 0

    ksrf(:ncol) = 0._kind_phys

    ! In most layers, no damping at all.
    tau_damp_rate(:ncol,:pver) = 0._kind_phys

  end subroutine vertical_diffusion_diffuse_horizontal_momentum_timestep_init

  ! Driver routine to compute vertical diffusion of horizontal momentum,
  ! including the loss of momentum due to frictional (KE to heat) dissipation.
  ! Turbulent diffusivities and boundary layer nonlocal transport terms are
  ! obtained from the turbulence module.
!> \section arg_table_vertical_diffusion_diffuse_horizontal_momentum_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_diffuse_horizontal_momentum_run.html
  subroutine vertical_diffusion_diffuse_horizontal_momentum_run( &
             ncol, pver, pverp, &
             dt, &
             rair, gravit, &
             do_iss, &
             am_correction, &
             itaures, &
             do_beljaars, &
             pint, t, rhoi, &
             taux, tauy, &
             tau_damp_rate, &
             ksrftms, &
             dragblj, &
             kvm, &
             dpidz_sq, &
             u0, v0, dse0, &
             tauresx, tauresy, &
             dtk, &
             tautmsx, tautmsy, &
             u1, v1, dse1, &
             errmsg, errflg)

    use atmos_phys_fin_vol_solver, only: fin_vol_coords_column, fin_vol_decomp_column, &
                                         fin_vol_left_div_column

    ! Input Arguments
    integer,         intent(in)       :: ncol
    integer,         intent(in)       :: pver
    integer,         intent(in)       :: pverp
    real(kind_phys), intent(in)       :: dt               ! physics timestep [s]
    real(kind_phys), intent(in)       :: rair
    real(kind_phys), intent(in)       :: gravit
    logical,         intent(in)       :: do_iss           ! Use implicit turbulent surface stress computation
    logical,         intent(in)       :: am_correction    ! Do angular momentum conservation correction
    logical,         intent(in)       :: itaures          ! Flag for updating tauresx tauresy in this subroutine.
    logical,         intent(in)       :: do_beljaars      ! Flag indicating Beljaars drag
    real(kind_phys), intent(in)       :: pint(:, :)       ! Air pressure at interfaces [Pa]
    real(kind_phys), intent(in)       :: t(:, :)          ! Temperature [K]
    real(kind_phys), intent(in)       :: rhoi(:, :)       ! Density of air at interfaces [kg m-3]
    real(kind_phys), intent(in)       :: taux(:)          ! Surface zonal stress [N m-2]
                                                          ! Input u-momentum per unit time per unit area into the atmosphere.
    real(kind_phys), intent(in)       :: tauy(:)          ! Surface meridional stress [N m-2]
                                                          ! Input v-momentum per unit time per unit area into the atmosphere.
    real(kind_phys), intent(in)       :: tau_damp_rate(:,:) ! Rate at which external (surface) stress damps wind speeds [s-1]
    real(kind_phys), intent(in)       :: ksrftms(:)       ! Surface drag coefficient for turbulent mountain stress. > 0. [kg m-2 s-1]
    real(kind_phys), intent(in)       :: dragblj(:, :)    ! Drag profile from Beljaars SGO form drag > 0. [s-1]
    real(kind_phys), intent(in)       :: kvm(:, :)        ! Eddy viscosity (Eddy diffusivity for momentum) [m^2 s-1], interfaces
    real(kind_phys), intent(in)       :: dpidz_sq(:,:)    ! Square of derivative of pressure with height (moist) [kg2 m-4 s-4], interfaces

    real(kind_phys), intent(in)       :: u0(:,:)          ! Input u-wind [m s-1]
    real(kind_phys), intent(in)       :: v0(:,:)          ! Input v-wind [m s-1]
    real(kind_phys), intent(in)       :: dse0(:,:)        ! Input Dry static energy [J kg-1]

    ! Input:  Reserved surface stress at previous time step
    ! Output: Reserved surface stress at current  time step
    real(kind_phys), intent(inout)    :: tauresx(:)
    real(kind_phys), intent(inout)    :: tauresy(:)

    ! Output Arguments
    real(kind_phys), intent(out)      :: dtk(:, :)        ! T tendency from KE dissipation [J kg-1]

    ! Recomputed TMS stresses using provisionally updated winds:
    real(kind_phys), intent(out)      :: tautmsx(:)       ! Implicit zonal turbulent mountain surface stress [N m-2]
    real(kind_phys), intent(out)      :: tautmsy(:)       ! Implicit meridional turbulent mountain surface stress [N m-2]

    ! Outputs from vertical diffusion will be converted to physics tendencies in separate scheme
    ! These provisional outputs are used to compute diagnostics, and thus need to be kept in non-tend form
    real(kind_phys), intent(out)      :: u1(:,:)          ! After vertical diffusion u-wind [m s-1]
    real(kind_phys), intent(out)      :: v1(:,:)          ! After vertical diffusion v-wind [m s-1]
    real(kind_phys), intent(out)      :: dse1(:,:)        ! After vertical diffusion Dry static energy [J kg-1]

    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    ! Local variables
    integer  :: i, k

    real(kind_phys) :: flux_to_state_conversion_factor(ncol)
    real(kind_phys) :: tmpi1(ncol, pverp)        ! Interface KE dissipation
    real(kind_phys) :: tmpi2(ncol, pverp)        ! dt*(g*rho)**2/dp at interfaces
    real(kind_phys) :: ramda                     ! dt/timeres [unitless]

    real(kind_phys) :: keg_in(ncol, pver)        ! KE on entry to subroutine
    real(kind_phys) :: keg_out(ncol, pver)       ! KE after U and V dissipation/diffusion

    real(kind_phys) :: tautotx(ncol)             ! Total surface stress (zonal)
    real(kind_phys) :: tautoty(ncol)             ! Total surface stress (meridional)

    real(kind_phys) :: dinp_u(ncol, pverp)       ! Vertical difference at interfaces, input u
    real(kind_phys) :: dinp_v(ncol, pverp)       ! Vertical difference at interfaces, input v
    real(kind_phys) :: dout_u                    ! Vertical difference at interfaces, output u
    real(kind_phys) :: dout_v                    ! Vertical difference at interfaces, output v

    real(kind_phys) :: usum_in(ncol)             ! Vertical integral of input u-momentum. Total zonal
                                                 ! momentum per unit area in column (sum of u*dp/g) [kg m s-1 m-2]
    real(kind_phys) :: vsum_in(ncol)             ! Vertical integral of input v-momentum. Total meridional
                                                 ! momentum per unit area in column (sum of v*dp/g) [kg m s-1 m-2]
    real(kind_phys) :: usum_out(ncol)            ! Vertical integral of u-momentum after doing implicit diffusion
    real(kind_phys) :: vsum_out(ncol)            ! Vertical integral of v-momentum after doing implicit diffusion

    real(kind_phys) :: tauimpx(ncol)             ! Actual net stress added at the current step other than mountain stress
    real(kind_phys) :: tauimpy(ncol)             ! Actual net stress added at the current step other than mountain stress
    real(kind_phys) :: taubljx(ncol)             ! recomputed explicit/residual beljaars stress
    real(kind_phys) :: taubljy(ncol)             ! recomputed explicit/residual beljaars stress

    ! Pressure coordinates (as Coords1D(pint)) and the tridiagonal decomposition
    real(kind_phys) :: mid(ncol, pver), del(ncol, pver), rdel(ncol, pver)
    real(kind_phys) :: rdst(ncol, pver-1)
    real(kind_phys) :: ca(ncol, pver), ze(ncol, pver), dnom(ncol, pver)

    errmsg = ''
    errflg = 0

    if (am_correction) then
      ! preserve time-mean torque
      ramda = 1._kind_phys
    else
      ramda = dt/timeres
    end if

    ! OpenACC: one column per thread, levels in order (each column's operations
    ! are those of the former whole-array code, in the same order). t is not used
    ! here; it is in the clause so that it is not taken for host data.
    !$acc data copyin(pint, t, rhoi, taux, tauy, tau_damp_rate, ksrftms, dragblj, kvm, &
    !$acc             dpidz_sq, u0, v0, dse0) &
    !$acc      copy(tauresx, tauresy) &
    !$acc      copyout(dtk, tautmsx, tautmsy, u1, v1, dse1) &
    !$acc      create(flux_to_state_conversion_factor, tmpi1, tmpi2, keg_in, keg_out, &
    !$acc             tautotx, tautoty, dinp_u, dinp_v, usum_in, vsum_in, usum_out, &
    !$acc             vsum_out, tauimpx, tauimpy, taubljx, taubljy, &
    !$acc             mid, del, rdel, rdst, ca, ze, dnom)

    !$acc parallel loop gang vector private(k, dout_u, dout_v)
    do i = 1, ncol
      call fin_vol_coords_column(i, pver, pint, mid, del, rdel, rdst)

      ! Set initial iterative outputs to input values
      !$acc loop seq
      do k = 1, pver
        u1(i,k) = u0(i,k)
        v1(i,k) = v0(i,k)
        dse1(i,k) = dse0(i,k)
      end do

      ! necessary temporaries used in computation
      flux_to_state_conversion_factor(i) = dt*gravit*rdel(i, pver)

      ! second term here uses dpidz_sq(:,1) as defined because
      ! the actual dpidz_sq term is "kept in only to preserve answers"
      ! but it is actually not used here ...
      tmpi2(i, 1) = dt*(gravit*rhoi(i,1)*gravit*rhoi(i,1))/(mid(i, 1) - pint(i, 1))
      !$acc loop seq
      do k = 2, pver
        tmpi2(i, k) = dt*dpidz_sq(i, k)*rdst(i, k-1)
      end do

      !---------------------------- !
      ! Diffuse Horizontal Momentum !
      !---------------------------- !
      ! Compute the vertical upward differences of the input u,v for KE dissipation
      ! at each interface.
      ! Velocity = 0 at surface, so difference at the bottom interface is -u,v(pver)
      ! These 'dinp_u, dinp_v' are computed using the non-diffused input wind.
      dinp_u(i, 1) = 0._kind_phys
      dinp_v(i, 1) = 0._kind_phys
      dinp_u(i, pver + 1) = -u0(i, pver)
      dinp_v(i, pver + 1) = -v0(i, pver)
      !$acc loop seq
      do k = 2, pver
        dinp_u(i, k) = u0(i, k) - u0(i, k - 1)
        dinp_v(i, k) = v0(i, k) - v0(i, k - 1)
      end do

      ! -------------------------------------------------------------- !
      ! Do 'Implicit Surface Stress' treatment for numerical stability !
      ! in the lowest model layer.                                     !
      ! -------------------------------------------------------------- !
      if (do_iss) then
        ! Vertical integration of input momentum.
        ! This is total horizontal momentum per unit area [ kg*m/s/m2 ] in each column.
        ! Note (u,v) are the raw input to the PBL scheme, not the
        ! provisionally-marched ones within the iteration loop of the PBL scheme.
        usum_in(i) = 0._kind_phys
        vsum_in(i) = 0._kind_phys
        !$acc loop seq
        do k = 1, pver
          usum_in(i) = usum_in(i) + (1._kind_phys/gravit)*u0(i, k)*del(i, k)
          vsum_in(i) = vsum_in(i) + (1._kind_phys/gravit)*v0(i, k)*del(i, k)
        end do

        ! Add residual stress of previous time step explicitly into the lowest
        ! model layer with a relaxation time scale of 'timeres' (ramda).
        u1(i, pver) = u1(i, pver) + flux_to_state_conversion_factor(i)*tauresx(i)*ramda
        v1(i, pver) = v1(i, pver) + flux_to_state_conversion_factor(i)*tauresy(i)*ramda

      else ! .not. do_iss
        ! In this case, do 'turbulent mountain stress' implicitly,
        ! but do 'normal turbulent stress' explicitly.
        ! In this case, there is no 'residual stress' as long as 'tms' is
        ! treated in a fully implicit way, which is true.

        ! Do 'normal stress' explicitly
        u1(i, pver) = u1(i, pver) + flux_to_state_conversion_factor(i)*taux(i)
        v1(i, pver) = v1(i, pver) + flux_to_state_conversion_factor(i)*tauy(i)
      end if  ! End of 'do iss' (implicit surface stress)

      ! ------------------------------------------------------------------------------------- !
      ! Diffuse horizontal momentum implicitly using tri-diagnonal matrix.                    !
      ! The 'u,v' are input-output: the output 'u,v' are implicitly diffused winds.           !
      !    For implicit 'normal' stress : ksrf = ksrftms + ksrfturb,                          !
      !                                   u(pver) : explicitly include 'residual normal' stress !
      !    For explicit 'normal' stress : ksrf = ksrftms                                      !
      !                                   u(pver) : explicitly include 'normal' stress        !
      ! Note that in all the two cases above, 'tms' is fully implicitly treated.              !
      ! ------------------------------------------------------------------------------------- !

      ! (fin_vol_solve for v1 and then u1 with coef_q=tau_damp_rate and
      ! coef_q_diff=kvm*dpidz_sq; both use the same operator, so it is
      ! decomposed once.)
      call fin_vol_decomp_column(i, pver, dt, del, rdel, rdst, kvm, dpidz_sq, &
                                 ca, ze, dnom, coef_q=tau_damp_rate)
      call fin_vol_left_div_column(i, pver, ca, ze, dnom, v1)
      call fin_vol_left_div_column(i, pver, ca, ze, dnom, u1)

      ! ---------------------------------------------------------------------- !
      ! Calculate 'total' ( tautotx ) and 'tms' ( tautmsx ) stresses that      !
      ! have been actually added into the atmosphere at the current time step. !
      ! Also, update residual stress, if required.                             !
      ! ---------------------------------------------------------------------- !

      ! Compute the implicit 'tms' using the updated winds.
      ! Below 'tautmsx(i),tautmsy(i)' are pure implicit mountain stresses
      ! that has been actually added into the atmosphere both for explicit
      ! and implicit approach.

      tautmsx(i) = -ksrftms(i)*u1(i, pver)
      tautmsy(i) = -ksrftms(i)*v1(i, pver)

      ! We want to add vertically-integrated Beljaars drag to residual stress.
      ! So this has to be calculated locally.
      ! We may want to rethink the residual drag calculation performed here on. (jtb)
      taubljx(i) = 0._kind_phys
      taubljy(i) = 0._kind_phys
      !$acc loop seq
      do k = 1, pver
        taubljx(i) = taubljx(i) + (1._kind_phys/gravit)*dragblj(i, k)*u1(i, k)*del(i, k)
        taubljy(i) = taubljy(i) + (1._kind_phys/gravit)*dragblj(i, k)*v1(i, k)*del(i, k)
      end do

      if (do_iss) then
        ! Compute vertical integration of final horizontal momentum
        usum_out(i) = 0._kind_phys
        vsum_out(i) = 0._kind_phys
        !$acc loop seq
        do k = 1, pver
          usum_out(i) = usum_out(i) + (1._kind_phys/gravit)*u1(i, k)*del(i, k)
          vsum_out(i) = vsum_out(i) + (1._kind_phys/gravit)*v1(i, k)*del(i, k)
        end do

        ! Compute net stress added into the atmosphere at the current time step.
        ! Note that the difference between 'usum_in' and 'usum_out' are induced
        ! by 'explicit residual stress + implicit total stress' for implicit case, while
        ! by 'explicit normal   stress + implicit tms   stress' for explicit case.
        ! Here, 'tautotx(i)' is net stress added into the air at the current time step.
        tauimpx(i) = (usum_out(i) - usum_in(i))/dt
        tauimpy(i) = (vsum_out(i) - vsum_in(i))/dt

        tautotx(i) = tauimpx(i)
        tautoty(i) = tauimpy(i)

        ! Compute residual stress and update if required.
        ! Note that the total stress we should have added at the current step is
        ! the sum of 'taux(i) - ksrftms(i)*u(i,pver) + tauresx(i)'.
        if (itaures) then
          tauresx(i) = taux(i) + tauresx(i) - tauimpx(i) + tautmsx(i) + taubljx(i)
          tauresy(i) = tauy(i) + tauresy(i) - tauimpy(i) + tautmsy(i) + taubljy(i)
        end if
      else
        ! .not. do_iss
        tautotx(i) = taux(i) + tautmsx(i)
        tautoty(i) = tauy(i) + tautmsy(i)
        tauresx(i) = 0._kind_phys
        tauresy(i) = 0._kind_phys
      end if  ! End of 'do_iss' if

      ! ------------------------------------ !
      ! Calculate kinetic energy dissipation !
      ! ------------------------------------ !

      ! Modification : In future, this should be set exactly same as
      !                the ones in the convection schemes

      if (do_beljaars) then
        ! 2. Add Kinetic Energy change across dissipation to Static Energy
        !$acc loop seq
        do k = 1, pver
          keg_in(i, k) = 0.5_kind_phys*(u0(i, k)*u0(i, k) + v0(i, k)*v0(i, k))
          keg_out(i, k) = 0.5_kind_phys*(u1(i, k)*u1(i, k) + v1(i, k)*v1(i, k))
          dtk(i, k) = keg_in(i, k) - keg_out(i, k)
          dse1(i, k) = dse1(i, k) + dtk(i, k) ! + dkeblj(i,k)
        end do
      else
        ! .not. do_beljaars

        ! 1. Compute dissipation term at interfaces
        !    Note that 'u,v' are already diffused wind, and 'tautotx,tautoty' are
        !    implicit stress that has been actually added. On the other hand,
        !    'dinp_u, dinp_v' were computed using non-diffused input wind.

        ! Modification : I should check whether non-consistency between 'u' and 'dinp_u'
        !                is correctly intended approach. I think so.
        tmpi1(i, 1) = 0._kind_phys
        tmpi1(i, pver + 1) = 0.5_kind_phys*dt*gravit* &
                      ((-u1(i, pver) + dinp_u(i, pver + 1))*tautotx(i) &
                       + (-v1(i, pver) + dinp_v(i, pver + 1))*tautoty(i))

        !$acc loop seq
        do k = 2, pver
          dout_u = u1(i, k) - u1(i, k - 1)
          dout_v = v1(i, k) - v1(i, k - 1)
          tmpi1(i, k) = 0.25_kind_phys*tmpi2(i, k)*kvm(i, k)* &
                        (dout_u**2 + dout_v**2 + dout_u*dinp_u(i, k) + dout_v*dinp_v(i, k))
        end do

        ! 2. Compute dissipation term at midpoints, add to dry static energy
        !$acc loop seq
        do k = 1, pver
          dtk(i, k) = (tmpi1(i, k + 1) + tmpi1(i, k))*rdel(i, k)
          dse1(i, k) = dse1(i, k) + dtk(i, k)
        end do
      end if
    end do

    !$acc end data

    ! End of kinetic energy dissipation due to horizontal momentum diffusion
  end subroutine vertical_diffusion_diffuse_horizontal_momentum_run

  ! Set dry static energy top boundary condition - zero
  ! used when (1) no molecular diffusion is active, or (2) molecular diffusion is active and not WACCM-X.
  !
  ! Note: this is an artifact to the way the vertical diffusion calculations are done;
  ! it does not represent physical reality.
!> \section arg_table_vertical_diffusion_set_dry_static_energy_at_toa_zero_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_set_dry_static_energy_at_toa_zero_run.html
  pure subroutine vertical_diffusion_set_dry_static_energy_at_toa_zero_run( &
    ncol, &
    dse_top, &
    errmsg, errflg)

    ! Input arguments
    integer,         intent(in)       :: ncol

    ! Output arguments
    real(kind_phys),    intent(out)   :: dse_top(:)       ! Dry static energy top boundary condition.
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    integer :: i

    errmsg = ''
    errflg = 0

    !$acc parallel loop copyout(dse_top)
    do i = 1, ncol
      dse_top(i) = 0._kind_phys
    end do

  end subroutine vertical_diffusion_set_dry_static_energy_at_toa_zero_run

  ! Set dry static energy top boundary condition - molecular diffusion non-WACCM-X version
!> \section arg_table_vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run.html
  pure subroutine vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run( &
    ncol, &
    gravit, &
    cpairv, &
    zi, &
    tint, &
    dse_top, &
    errmsg, errflg)

    ! Input arguments
    integer,         intent(in)       :: ncol
    real(kind_phys), intent(in)       :: gravit
    real(kind_phys), intent(in)       :: cpairv(:,:)    ! specific heat of dry air at constant pressure [J kg-1 K-1]
    real(kind_phys), intent(in)       :: zi(:,:)        ! geopotential_height_wrt_surface_at_interface [m]
    real(kind_phys), intent(in)       :: tint(:,:)      ! Air temperature at interfaces [K]

    ! Output arguments
    real(kind_phys),    intent(out)   :: dse_top(:)     ! Dry static energy top boundary condition.
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    errmsg = ''
    errflg = 0

    dse_top(:ncol) = cpairv(:ncol,1) * tint(:ncol,1) + gravit * zi(:ncol,1)

  end subroutine vertical_diffusion_set_dry_static_energy_at_toa_molecdiff_run

  ! Driver routine to compute vertical diffusion of dry static energy.
  ! The new temperature is computed from the diffused dry static energy
  ! after tendencies are applied.
  ! Turbulent diffusivities and boundary layer nonlocal transport terms are
  ! obtained from the turbulence module.
!> \section arg_table_vertical_diffusion_diffuse_dry_static_energy_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_diffuse_dry_static_energy_run.html
  subroutine vertical_diffusion_diffuse_dry_static_energy_run( &
    ncol, pver, &
    dt, &
    gravit, &
    pint, rhoi, &
    shflx, &
    dse_top, &
    kvh, cgh, &
    dpidz_sq, &
    ! input/output (may have been updated by horizontal momentum diffusion)
    dse, &
    errmsg, errflg)

    use atmos_phys_fin_vol_solver, only: fin_vol_coords_column, fin_vol_decomp_column, &
                                         fin_vol_left_div_column

    ! Input Arguments
    integer,         intent(in)       :: ncol
    integer,         intent(in)       :: pver
    real(kind_phys), intent(in)       :: dt               ! physics timestep [s]

    real(kind_phys), intent(in)       :: gravit
    real(kind_phys), intent(in)       :: pint(:, :)       ! Air pressure at interfaces [Pa]
    real(kind_phys), intent(in)       :: rhoi(:, :)       ! Density of air at interfaces [kg m-3]
    real(kind_phys), intent(in)       :: shflx(:)         ! Surface sensible heat flux [W m-2]
    real(kind_phys), intent(in)       :: dse_top(:)       ! Dry static energy top boundary condition.
    real(kind_phys), intent(in)       :: kvh(:, :)        ! Eddy diffusivity for heat [m^2 s-1], interfaces
    real(kind_phys), intent(in)       :: cgh(:, :)        ! Counter-gradient term for heat [J kg-1 m-1], interfaces
    real(kind_phys), intent(in)       :: dpidz_sq(:,:)    ! Square of derivative of pressure with height (moist) [kg2 m-4 s-4]

    ! Input/output arguments
    real(kind_phys), intent(inout)    :: dse(:,:)         ! After vertical diffusion Dry static energy [J kg-1]
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    ! Local variables
    integer :: i, k

    real(kind_phys) :: flux_to_state_conversion_factor(ncol)

    ! Pressure coordinates (as Coords1D(pint)) and the tridiagonal decomposition
    real(kind_phys) :: mid(ncol, pver), del(ncol, pver), rdel(ncol, pver)
    real(kind_phys) :: rdst(ncol, pver-1)
    real(kind_phys) :: ca(ncol, pver), ze(ncol, pver), dnom(ncol, pver)

    errmsg = ''
    errflg = 0

    ! OpenACC: one column per thread, levels in order.
    !$acc data copyin(pint, rhoi, shflx, dse_top, kvh, cgh, dpidz_sq) copy(dse) &
    !$acc      create(flux_to_state_conversion_factor, mid, del, rdel, rdst, ca, ze, dnom)

    !$acc parallel loop gang vector private(k)
    do i = 1, ncol
      call fin_vol_coords_column(i, pver, pint, mid, del, rdel, rdst)

      flux_to_state_conversion_factor(i) = dt*gravit*rdel(i, pver)

      ! Modification : In future, we should diffuse the fully conservative
      !                moist static energy, not the dry static energy.

      ! Add counter-gradient to input static energy profiles
      !$acc loop seq
      do k = 1, pver
        dse(i, k) = dse(i, k) + dt*rdel(i, k)*gravit* &
                    (rhoi(i, k + 1)*kvh(i, k + 1)*cgh(i, k + 1) &
                     - rhoi(i, k)*kvh(i, k)*cgh(i, k))
      end do

      ! Add the explicit surface fluxes to the lowest layer
      dse(i, pver) = dse(i, pver) + flux_to_state_conversion_factor(i)*shflx(i)

      !---------------------------------------------------
      ! Solve for temperature using thermal conductivity
      !---------------------------------------------------

      ! Boundary layer thickness of "0._kind_phys" signifies that the boundary
      ! condition is defined directly on the top interface.
      ! (fin_vol_solve with coef_q_diff=kvh*dpidz_sq, upper_bndry a fixed layer
      ! of width 0 (BoundaryFixedLayer) and l_cond=BoundaryData(dse_top).)
      call fin_vol_decomp_column(i, pver, dt, del, rdel, rdst, kvh, dpidz_sq, &
                                 ca, ze, dnom, top_width=0._kind_phys)
      call fin_vol_left_div_column(i, pver, ca, ze, dnom, dse, top_value=dse_top(i))
    end do

    !$acc end data

  end subroutine vertical_diffusion_diffuse_dry_static_energy_run

  ! Set tracers to diffuse using vertical diffusion.
!> \section arg_table_vertical_diffusion_diffuse_tracers_init Argument Table
!! \htmlinclude arg_table_vertical_diffusion_diffuse_tracers_init.html
  subroutine vertical_diffusion_diffuse_tracers_init( &
    ncnst, &
    const_props, &
    do_diffusion_const, &
    errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! Input arguments
    integer,            intent(in)    :: ncnst                 ! # of constituents
    type(ccpp_constituent_prop_ptr_t), &
                        intent(in)    :: const_props(:)        ! CCPP constituent properties pointer

    ! Output arguments
    logical,            intent(out)   :: do_diffusion_const(:) ! diffuse constituents (size ncnst) [flag]
    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    integer :: m
    character(len=256) :: cnst_name

    errmsg = ''
    errflg = 0

    ! Check if logical array is of the expected size
    if (size(do_diffusion_const) .ne. ncnst) then
      errflg = 1
      write(errmsg, *) 'compute_vdiff: do_diffusion_const size ', &
        size(do_diffusion_const), ' is not equal to ncnst ', ncnst
      return
    end if

    ! Diffuse all constituents using moist coordinates.
    do_diffusion_const(:ncnst) = .true.

    ! FIXME: when dropmixnuc is CCPPized, it should flag if a constituent
    ! is treated there (i.e., vertically mixed by ndrop activation process)
    ! and use this flag to skip handling of vertical diffusion in this module.
    !
    ! For now, we replicate the chemistry.F90 code that skips modal aerosol
    ! species (which undergo ndrop activation mixing) but we do not want to
    ! do this specific handling here...
    do m = 1, ncnst
      call const_props(m)%standard_name(cnst_name, errflg, errmsg)
      if(errflg /= 0) return

      ! cannot match just _a because it would match something like
      ! moist_air_and_condensed_water ...
      if(index(cnst_name, '_a1') > 0 .or. &     ! MAM modal species.
         index(cnst_name, '_a2') > 0 .or. &     ! MAM modal species.
         index(cnst_name, '_a3') > 0 .or. &     ! MAM modal species.
         index(cnst_name, '_a4') > 0 .or. &     ! MAM modal species.
         index(cnst_name, '_c1') > 0 .or. &     ! MAM cloud-borne species.
         index(cnst_name, '_c2') > 0 .or. &     ! MAM cloud-borne species.
         index(cnst_name, '_c3') > 0 .or. &     ! MAM cloud-borne species.
         index(cnst_name, '_c4') > 0 .or. &     ! MAM cloud-borne species.
         index(cnst_name, 'NUMLIQ') > 0 .or. &  ! micro_pumas_cam if modal/trop_strat carma
         index(cnst_name, 'mass_number_concentration_of_cloud_liquid_water_droplets_in_moist_air_and_condensed_water') > 0) then
         ! micro_pumas_cam (including MG1) if modal/trop_strat carma

        ! is a modal aerosol species (or numliq) undergoing ndrop activation mixing.
        ! in this case, vertical diffusion does not diffuse it.
        ! cloud-borne (_c) modal species are excluded as well: in CAM they are
        ! pbuf fields (qqcw), not constituents, so vertical diffusion never
        ! sees them; in CAM-SIMA they are registered as constituents.
        do_diffusion_const(m) = .false.
      end if
    end do

  end subroutine vertical_diffusion_diffuse_tracers_init

  ! Diffuse tracers (no molecular diffusion)
!> \section arg_table_vertical_diffusion_diffuse_tracers_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_diffuse_tracers_run.html
  subroutine vertical_diffusion_diffuse_tracers_run( &
    ncol, pver, &
    ncnst, &
    dt, &
    rair, &
    gravit, &
    do_diffusion_const, &
    pint, &
    t, &
    rhoi, &
    cflx, &
    kvh, kvq, cgs, &
    qmincg, &
    dpidz_sq, &
    ubc_mmr, &
    cnst_fixed_ubc, &
    q0, &
    ! below output
    q1, &
    errmsg, errflg)

    use atmos_phys_fin_vol_solver, only: fin_vol_coords_column, fin_vol_decomp_column, &
                                         fin_vol_left_div_column_3d

    ! Input arguments
    integer,         intent(in)       :: ncol             ! Number of atmospheric columns
    integer,         intent(in)       :: pver
    integer,         intent(in)       :: ncnst            ! # of constituents to diffuse. In eddy_diff, only wv. Others, all advected CCPP constituents.
    real(kind_phys), intent(in)       :: dt               ! physics timestep [s]
    real(kind_phys), intent(in)       :: rair
    real(kind_phys), intent(in)       :: gravit
    logical,         intent(in)       :: do_diffusion_const(:) ! diffuse constituents (size ncnst) [flag]

    real(kind_phys), intent(in)       :: pint(:, :)       ! Air pressure at interfaces [Pa]
    real(kind_phys), intent(in)       :: t(:, :)          ! Temperature [K]
    real(kind_phys), intent(in)       :: rhoi(:, :)       ! Density of air at interfaces [kg m-3]

    real(kind_phys), intent(in)       :: cflx(:, :)       ! Surface constituent flux [kg m-2 s-1] (ncol,ncnst)

    real(kind_phys), intent(in)       :: kvh(:, :)        ! Eddy diffusivity for heat [m^2 s-1], interfaces
    real(kind_phys), intent(in)       :: kvq(:, :)        ! Eddy diffusivity for constituents [m^2 s-1], interfaces
    real(kind_phys), intent(in)       :: cgs(:, :)        ! Counter-gradient star [cg/flux] [s m-2], interfaces

    real(kind_phys), intent(in)       :: qmincg(:)        ! Minimum constituent mixing ratios from cg fluxes (ncnst)
    real(kind_phys), intent(in)       :: dpidz_sq(:,:)    ! Square of derivative of pressure with height (moist) [kg2 m-4 s-4]

    ! Upper boundary properties (for when molecular diffusion is off)
    real(kind_phys), intent(in)       :: ubc_mmr(:, :)    ! Upper boundary condition of constituents [none]
    logical,         intent(in)       :: cnst_fixed_ubc(:)! Whether upper boundary condition is fixed

    real(kind_phys), intent(in)       :: q0(:,:,:)        ! Input Constituents [kg kg-1]

    ! Output arguments
    real(kind_phys), intent(out)      :: q1(:,:,:)        ! After vertical diffusion Constituents [kg kg-1]


    character(len=512), intent(out)   :: errmsg  ! error message
    integer,            intent(out)   :: errflg  ! error flag

    ! Local variables
    integer :: i, m, k
    real(kind_phys) :: flux_to_state_conversion_factor(ncol)
    logical  :: lqtst(ncol)                      ! Adjust vertical profiles
    real(kind_phys) :: qtm(ncol, pver)           ! Temporary copy of q
    real(kind_phys) :: rrho(ncol)                ! 1./bottom level density

    ! Pressure coordinates (as Coords1D(pint)) and the tridiagonal decomposition
    real(kind_phys) :: mid(ncol, pver), del(ncol, pver), rdel(ncol, pver)
    real(kind_phys) :: rdst(ncol, pver-1)
    real(kind_phys) :: ca(ncol, pver), ze(ncol, pver), dnom(ncol, pver)

    errmsg = ''
    errflg = 0

    ! OpenACC: one column per thread, levels and constituents in order (each
    ! column's operations are those of the former whole-array code, in the same
    ! order).
    !$acc data copyin(pint, t, rhoi, cflx, kvh, kvq, cgs, qmincg, dpidz_sq, ubc_mmr, &
    !$acc             do_diffusion_const, cnst_fixed_ubc, q0) &
    !$acc      copyout(q1) &
    !$acc      create(flux_to_state_conversion_factor, lqtst, qtm, rrho, &
    !$acc             mid, del, rdel, rdst, ca, ze, dnom)

    !$acc parallel loop gang vector private(k, m)
    do i = 1, ncol
      call fin_vol_coords_column(i, pver, pint, mid, del, rdel, rdst)

      ! Set initial iterative outputs to input values
      !$acc loop seq
      do m = 1, ncnst
        !$acc loop seq
        do k = 1, pver
          q1(i,k,m) = q0(i,k,m)
        end do
      end do

      ! necessary temporaries used in computation
      rrho(i) = rair*t(i, pver)/mid(i, pver)
      flux_to_state_conversion_factor(i) = dt*gravit*rdel(i, pver)

      ! LU decomposition (fin_vol_lu_decomp with coef_q_diff=kvq*dpidz_sq)
      call fin_vol_decomp_column(i, pver, dt, del, rdel, rdst, kvq, dpidz_sq, ca, ze, dnom)

      ! Loop through constituents
      !$acc loop seq
      do m = 1, ncnst
        if (do_diffusion_const(m)) then
          ! Add the nonlocal transport terms to constituents in the PBL.
          ! Check for neg q's in each constituent and put the original vertical
          ! profile back if a neg value is found. A neg value implies that the
          ! quasi-equilibrium conditions assumed for the countergradient term are
          ! strongly violated.
          !$acc loop seq
          do k = 1, pver
            qtm(i, k) = q1(i, k, m)
          end do

          !$acc loop seq
          do k = 1, pver
            q1(i, k, m) = q1(i, k, m) + &
                          dt*rdel(i, k)*gravit*(cflx(i, m)*rrho(i))* &
                          (rhoi(i, k + 1)*kvh(i, k + 1)*cgs(i, k + 1) &
                           - rhoi(i, k)*kvh(i, k)*cgs(i, k))
          end do
          ! lqtst = all(q1(i,1:pver,m) >= qmincg(m)); else the original profile
          lqtst(i) = .true.
          !$acc loop seq
          do k = 1, pver
            if (.not. (q1(i, k, m) >= qmincg(m))) lqtst(i) = .false.
          end do
          if (.not. lqtst(i)) then
            !$acc loop seq
            do k = 1, pver
              q1(i, k, m) = qtm(i, k)
            end do
          end if

          ! Add the explicit surface fluxes to the lowest layer
          q1(i, pver, m) = q1(i, pver, m) + flux_to_state_conversion_factor(i)*cflx(i, m)

          ! not doing molecular diffusion
          ! explicitly set mmr in top layer for cases where molecular diffusion is not active
          if (cnst_fixed_ubc(m)) then
            ! assumes TOA is at vertical level 1 here.
            q1(i, 1, m) = ubc_mmr(i, m)
          end if
          call fin_vol_left_div_column_3d(i, pver, m, ca, ze, dnom, q1)
        end if
      end do
    end do

    !$acc end data

  end subroutine vertical_diffusion_diffuse_tracers_run

  ! Based on provisional output from vertical diffusion,
  ! calculate final physics tendencies to be applied by the tendency updater.
!> \section arg_table_vertical_diffusion_tendencies_run Argument Table
!! \htmlinclude arg_table_vertical_diffusion_tendencies_run.html
  subroutine vertical_diffusion_tendencies_run( &
    ncol, pver, pcnst, &
    const_props, &
    dt, &
    pdel, pdeldry, &
    u0, v0, s0, q0, &  ! actual values at beginning of vdiff
    u1, v1, s1, q1, &  ! provisional values after vdiff, not actual
    ! below output
    tend_s, tend_u, tend_v, tend_q, &
    scheme_name, &
    errmsg, errflg)

    ! framework dependency for const_props
    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! Input arguments
    integer,         intent(in)     :: ncol            ! Number of atmospheric columns [count]
    integer,         intent(in)     :: pver            ! Number of vertical levels [count]
    integer,         intent(in)     :: pcnst           ! Number of ccpp constituents [count]
    type(ccpp_constituent_prop_ptr_t), &
                     intent(in)     :: const_props(:)  ! CCPP constituent properties pointer
    real(kind_phys), intent(in)     :: dt              ! Timestep [s]
    real(kind_phys), intent(in)     :: pdel(:,:)       ! Pressure thickness of layers [Pa]
    real(kind_phys), intent(in)     :: pdeldry(:,:)    ! Dry pressure thickness of layers [Pa]
    real(kind_phys), intent(in)     :: u0(:,:)         ! Initial eastward wind [m s-1]
    real(kind_phys), intent(in)     :: v0(:,:)         ! Initial northward wind [m s-1]
    real(kind_phys), intent(in)     :: s0(:,:)         ! Initial dry static energy [J kg-1]
    real(kind_phys), intent(in)     :: q0(:,:,:)       ! Initial constituent mixing ratios [kg kg-1]
    real(kind_phys), intent(in)     :: u1(:,:)          ! Provisional eastward wind after diffusion [m s-1]
    real(kind_phys), intent(in)     :: v1(:,:)          ! Provisional northward wind after diffusion [m s-1]
    real(kind_phys), intent(in)     :: s1(:,:)          ! Provisional dry static energy after diffusion [J kg-1]
    real(kind_phys), intent(in)     :: q1(:,:,:)        ! Provisional constituent mixing ratios after diffusion [kg kg-1]

    ! Output arguments
    real(kind_phys), intent(out)    :: tend_s(:,:)     ! Tendency of dry air enthalpy at constant pressure [J kg-1 s-1]
    real(kind_phys), intent(out)    :: tend_u(:,:)     ! Eastward wind tendency [m s-2]
    real(kind_phys), intent(out)    :: tend_v(:,:)     ! Northward wind tendency [m s-2]
    real(kind_phys), intent(out)    :: tend_q(:,:,:)   ! Constituent mixing ratio tendencies [kg kg-1 s-1]
    character(len=64),  intent(out) :: scheme_name    ! scheme name
    character(len=512), intent(out) :: errmsg     ! Error message
    integer,            intent(out) :: errflg     ! Error flag

    integer :: i, k, m
    logical :: const_is_dry(pcnst)

    ! for bit-to-bitness with CAM.
    real(kind_phys) :: rdt

    scheme_name = "vertical_diffusion_tendencies"

    rdt = 1._kind_phys/dt

    errmsg = ''
    errflg = 0

    ! Which constituents are dry (on the host: const_props is a host object)
    do m = 1, pcnst
      call const_props(m)%is_dry(const_is_dry(m), errflg, errmsg)
    enddo

    !$acc data copyin(pdel, pdeldry, u0, v0, s0, q0, u1, v1, s1, q1, const_is_dry) &
    !$acc      copyout(tend_s, tend_u, tend_v, tend_q)

    ! calculate physics tendencies
    !$acc parallel loop collapse(2)
    do k = 1, pver
      do i = 1, ncol
        tend_s(i,k) = (s1(i,k) - s0(i,k)) * rdt
        tend_u(i,k) = (u1(i,k) - u0(i,k)) * rdt
        tend_v(i,k) = (v1(i,k) - v0(i,k)) * rdt
      end do
    end do
    !$acc parallel loop collapse(3)
    do m = 1, size(tend_q, 3)
      do k = 1, pver
        do i = 1, ncol
          tend_q(i,k,m) = (q1(i,k,m) - q0(i,k,m)) * rdt
        end do
      end do
    end do

    ! convert tendencies of dry constituents to dry basis
    !$acc parallel loop collapse(3)
    do m = 1, pcnst
      do k = 1, pver
        do i = 1, ncol
          if (const_is_dry(m)) then
            tend_q(i,k,m) = tend_q(i,k,m) * pdel(i,k) / pdeldry(i,k)
          end if
        end do
      end do
    end do

    !$acc end data

  end subroutine vertical_diffusion_tendencies_run

end module diffusion_solver
