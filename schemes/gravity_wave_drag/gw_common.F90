! This module contains code common to different gravity wave
! parameterizations.
module gw_common

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  ! Public CCPP-compliant initialization interface.
  ! This scheme also reads in all namelist parameters for all gravity wave
  ! parameterizations.
  public :: gravity_wave_drag_common_init

  ! Public interfaces for use by multiple gravity wave parameterizations.
  public :: GWBand

  public :: gw_drag_prof
  public :: qbo_hdepth_scaling
  public :: calc_taucd, momentum_flux, momentum_fixer
  public :: energy_change, energy_fixer
  public :: gw_drag_prof_column, energy_change_column
  public :: coriolis_speed, adjust_inertial
  public :: find_bin

  public :: west, east, north, south
  public :: wavelength_mid, wavelength_long

  real(kind_phys), public, parameter :: unset_kind_phys = huge(1._kind_phys)

  ! Number of levels and interfaces in the atmosphere.
  integer, protected :: pver = 0
  integer, protected :: pverp = 0

  ! Whether or not to enforce an upper boundary condition of tau = 0.
  logical :: tau_0_ubc = .false.

  ! Index the cardinal directions.
  integer, parameter :: west = 1
  integer, parameter :: east = 2
  integer, parameter :: south = 3
  integer, parameter :: north = 4

  ! Number of seconds per day [s]
  real(kind_phys), parameter :: seconds_per_day = 86400._kind_phys

  ! Scaling factor for generating QBO
  real(kind_phys), protected :: qbo_hdepth_scaling

  ! Physical constants from host model.
  real(kind_phys), protected :: pi     = huge(1._kind_phys)
  real(kind_phys), protected :: gravit = huge(1._kind_phys)
  real(kind_phys), protected :: rair   = huge(1._kind_phys)

  real(kind_phys), protected :: omega_earth = huge(1._kind_phys)

  ! Horizontal wavelengths [m].
  real(kind_phys), parameter :: wavelength_mid = 1.e5_kind_phys
  real(kind_phys), parameter :: wavelength_long = 3.e5_kind_phys

  ! Definition of the bin boundaries.
  real(kind_phys), parameter :: bounds(4) = [-40._kind_phys, -15._kind_phys, &
       15._kind_phys, 40._kind_phys]

  ! Private variables

  ! Top level for gravity wave sources.
  integer,         parameter :: ktop = 1

  ! Background diffusivity.
  real(kind_phys), parameter :: dback = 0.05_kind_phys

  ! rair/gravit
  real(kind_phys) :: rog = huge(1._kind_phys)

  ! Newtonian cooling coefficients.
  real(kind_phys), allocatable :: alpha(:)

  ! Inverse Prandtl number.
  real(kind_phys) :: prndl

  ! Module data read by gw_drag_prof_column/energy_change_column on the
  ! device; set in gravity_wave_drag_common_init (update device at its end).
  !$acc declare create(gravit, rog, alpha, prndl, tau_0_ubc)

  !
  ! Limits to keep values reasonable.
  !

  ! Minimum non-zero stress.
  real(kind_phys), parameter :: taumin = 1.e-10_kind_phys
  ! Maximum wind tendency from stress divergence (before efficiency applied).
  ! 400 m/s/day
  real(kind_phys), parameter :: tndmax = 400._kind_phys/seconds_per_day
  ! Maximum allowed change in u-c (before efficiency applied).
  real(kind_phys), parameter :: umcfac = 0.5_kind_phys
  ! Minimum value of (u-c)**2.
  real(kind_phys), parameter :: ubmc2mn = 0.01_kind_phys

  ! Type describing a band of wavelengths into which gravity waves can be emitted.
  ! Currently this has to have uniform spacing (i.e. adjacent elements of cref are exactly dc apart).
  type :: GWBand
    ! Dimension of the spectrum.
    integer :: ngwv
    ! Delta between nearest phase speeds [m/s].
    real(kind_phys) :: dc
    ! Reference speeds [m/s].
    real(kind_phys), allocatable :: cref(:)
    ! Critical Froude number, squared (usually 1, but CAM3 used 0.5).
    real(kind_phys) :: fcrit2
    ! Horizontal wave number [1/m].
    real(kind_phys) :: kwv
    ! Effective horizontal wave number [1/m] (fcrit2*kwv).
    real(kind_phys) :: effkwv
  end type GWBand

  interface GWBand
    module procedure new_GWBand
  end interface

contains

!==========================================================================

  ! Constructor for a GWBand that calculates derived components.
  function new_GWBand(ngwv, dc, fcrit2, wavelength) result(band)
    ! Used directly to set the type's components.
    integer, intent(in) :: ngwv
    real(kind_phys), intent(in) :: dc
    real(kind_phys), intent(in) :: fcrit2
    ! Wavelength in meters.
    real(kind_phys), intent(in) :: wavelength

    ! Output.
    type(GWBand) :: band

    ! Wavenumber index.
    integer :: l

    ! Simple assignments.
    band%ngwv = ngwv
    band%dc = dc
    ! For now just ensure fcrit is always set to 1
    band%fcrit2 = 1.0_kind_phys ! fcrit2

    ! Uniform phase speed reference grid.
    allocate (band%cref(-ngwv:ngwv))
    band%cref = [(dc*l, l=-ngwv, ngwv)]

    ! Wavenumber and effective wavenumber come from the wavelength.
    band%kwv = 2._kind_phys*pi/wavelength
    band%effkwv = band%fcrit2*band%kwv

  end function new_GWBand

  ! Common initialization.
!> \section arg_table_gravity_wave_drag_common_init Argument Table
!! \htmlinclude arg_table_gravity_wave_drag_common_init.html
  subroutine gravity_wave_drag_common_init( &
             pver_in, &
             pverp_in, &
             amIRoot, iulog, &
             pref_edge, &
             tau_0_ubc_in, &
             pi_in, gravit_in, rair_in, &
             prndl_in, qbo_hdepth_scaling_in, &
             errmsg, errflg)

    use atmos_phys_interp_utils, only: lininterp

    integer,          intent(in)  :: pver_in
    integer,          intent(in)  :: pverp_in
    logical,          intent(in)  :: amIRoot
    integer,          intent(in)  :: iulog
    real(kind_phys),  intent(in)  :: pref_edge(:)
    logical,          intent(in)  :: tau_0_ubc_in
    real(kind_phys),  intent(in)  :: pi_in
    real(kind_phys),  intent(in)  :: gravit_in
    real(kind_phys),  intent(in)  :: rair_in
    real(kind_phys),  intent(in)  :: prndl_in
    real(kind_phys),  intent(in)  :: qbo_hdepth_scaling_in

    character(len=512), intent(out) :: errmsg
    integer, intent(out)            :: errflg

    integer :: k

    ! Levels of pre-calculated Newtonian cooling (1/day).
    ! The following profile is digitized from:
    ! Wehrbein and Leovy (JAS, 39, 1532-1544, 1982) Figure 5
    ! https://doi.org/10.1175/1520-0469(1982)039<1532:AARHAC>2.0.CO;2
    integer, parameter :: nalph = 71
    real(kind_phys) :: alpha0(nalph) = [ &
                       0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, &
                       0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, 0.1_kind_phys, &
                       0.1_kind_phys, 0.1_kind_phys, 0.10133333_kind_phys, 0.104_kind_phys, &
                       0.108_kind_phys, 0.112_kind_phys, 0.116_kind_phys, 0.12066667_kind_phys, &
                       0.126_kind_phys, 0.132_kind_phys, 0.138_kind_phys, 0.144_kind_phys, &
                       0.15133333_kind_phys, 0.16_kind_phys, 0.17_kind_phys, 0.18_kind_phys, &
                       0.19_kind_phys, 0.19933333_kind_phys, 0.208_kind_phys, 0.216_kind_phys, &
                       0.224_kind_phys, 0.232_kind_phys, 0.23466667_kind_phys, 0.232_kind_phys, &
                       0.224_kind_phys, 0.216_kind_phys, 0.208_kind_phys, 0.20133333_kind_phys, &
                       0.196_kind_phys, 0.192_kind_phys, 0.188_kind_phys, 0.184_kind_phys, &
                       0.18266667_kind_phys, 0.184_kind_phys, 0.188_kind_phys, 0.192_kind_phys, &
                       0.196_kind_phys, 0.19333333_kind_phys, 0.184_kind_phys, 0.168_kind_phys, &
                       0.152_kind_phys, 0.136_kind_phys, 0.12133333_kind_phys, 0.108_kind_phys, &
                       0.096_kind_phys, 0.084_kind_phys, 0.072_kind_phys, 0.061_kind_phys, &
                       0.051_kind_phys, 0.042_kind_phys, 0.033_kind_phys, 0.024_kind_phys, &
                       0.017666667_kind_phys, 0.014_kind_phys, 0.013_kind_phys, 0.012_kind_phys, &
                       0.011_kind_phys, 0.010333333_kind_phys, 0.01_kind_phys, 0.01_kind_phys, &
                       0.01_kind_phys, 0.01_kind_phys, 0.01_kind_phys]

    ! Pressure levels that were used to calculate alpha0 (hPa).
    real(kind_phys) :: palph(nalph) = [ &
                       2.06115E-06_kind_phys, 2.74280E-06_kind_phys, 3.64988E-06_kind_phys, 4.85694E-06_kind_phys, &
                       6.46319E-06_kind_phys, 8.60065E-06_kind_phys, 1.14450E-05_kind_phys, 1.52300E-05_kind_phys, &
                       2.02667E-05_kind_phys, 2.69692E-05_kind_phys, 3.58882E-05_kind_phys, 4.77568E-05_kind_phys, &
                       6.35507E-05_kind_phys, 8.45676E-05_kind_phys, 0.000112535_kind_phys, 0.000149752_kind_phys, &
                       0.000199277_kind_phys, 0.000265180_kind_phys, 0.000352878_kind_phys, 0.000469579_kind_phys, &
                       0.000624875_kind_phys, 0.000831529_kind_phys, 0.00110653_kind_phys, 0.00147247_kind_phys, &
                       0.00195943_kind_phys, 0.00260744_kind_phys, 0.00346975_kind_phys, 0.00461724_kind_phys, &
                       0.00614421_kind_phys, 0.00817618_kind_phys, 0.0108801_kind_phys, 0.0144783_kind_phys, &
                       0.0192665_kind_phys, 0.0256382_kind_phys, 0.0341170_kind_phys, 0.0453999_kind_phys, &
                       0.0604142_kind_phys, 0.0803939_kind_phys, 0.106981_kind_phys, 0.142361_kind_phys, &
                       0.189442_kind_phys, 0.252093_kind_phys, 0.335463_kind_phys, 0.446404_kind_phys, &
                       0.594036_kind_phys, 0.790490_kind_phys, 1.05192_kind_phys, 1.39980_kind_phys, &
                       1.86273_kind_phys, 2.47875_kind_phys, 3.29851_kind_phys, 4.38936_kind_phys, &
                       5.84098_kind_phys, 7.77266_kind_phys, 10.3432_kind_phys, 13.7638_kind_phys, &
                       18.3156_kind_phys, 24.3728_kind_phys, 32.4332_kind_phys, 43.1593_kind_phys, &
                       57.4326_kind_phys, 76.4263_kind_phys, 101.701_kind_phys, 135.335_kind_phys, &
                       180.092_kind_phys, 239.651_kind_phys, 318.907_kind_phys, 424.373_kind_phys, &
                       564.718_kind_phys, 751.477_kind_phys, 1000._kind_phys]

    errmsg = ''
    errflg = 0

    pver = pver_in
    pverp = pverp_in
    tau_0_ubc = tau_0_ubc_in

    pi = pi_in
    gravit = gravit_in
    rair = rair_in

    omega_earth = 2._kind_phys*pi/seconds_per_day

    ! Interpolate Newtonian cooling to model interface levels
    allocate (alpha(pver + 1), stat=errflg, errmsg=errmsg)
    if (errflg /= 0) return

    ! pre-calculated newtonian damping:
    !     * convert to s-1
    !     * ensure it is not smaller than 1e-6
    !     * convert palph from hpa to pa

    do k = 1, nalph
      alpha0(k) = alpha0(k)/seconds_per_day
      alpha0(k) = max(alpha0(k), 1.e-6_kind_phys)
      palph(k) = palph(k)*1.e2_kind_phys
    end do

    call lininterp(alpha0, palph, nalph, alpha, pref_edge, pver + 1, errmsg, errflg)
    if (errflg /= 0) then
      errmsg = 'gravity_wave_drag_common_init: ' // trim(errmsg)
      return
    end if

    if (amIRoot) then
      write (iulog, *) 'gravity_wave_drag_common_init: newtonian damping (1/day):'
      write (iulog, fmt='(a4,a12,a10)') ' k  ', '  pref_edge      ', '  alpha   '
      do k = 1, pver + 1
        write (iulog, fmt='(i4,1e12.5,1f10.2)') k, pref_edge(k), &
          alpha(k)*seconds_per_day
      end do

      write (iulog, *) 'gravity_wave_drag_common_init: ktop = ', ktop
    end if

    prndl = prndl_in
    qbo_hdepth_scaling = qbo_hdepth_scaling_in

    rog = rair/gravit

    !$acc update device(gravit, rog, alpha, prndl, tau_0_ubc)

  end subroutine gravity_wave_drag_common_init

!==========================================================================

  ! Solve for the drag profile from the multiple gravity wave drag
  ! parameterization, for all columns: a host loop over
  ! gw_drag_prof_column (see there). Kept for the schemes that still pass
  ! a Coords1D and a GWBand.
  subroutine gw_drag_prof(ncol, band, p, src_level, tend_level, dt, &
                          t, vramp, &
                          piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                          effgw, c, kvtt, q, dse, tau, utgw, vtgw, &
                          ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                          ro_adjust, kwvrdg, satfac_in, lapply_effgw_in, lapply_vdiff, tau_diag)

    use coords_1d, only: coords1d

    !------------------------------Arguments--------------------------------
    ! Column dimension.
    integer, intent(in) :: ncol
    ! Wavelengths.
    type(GWBand), intent(in) :: band
    ! Pressure coordinates.
    type(coords1d), intent(in) :: p
    ! Level from which gravity waves are propagated upward.
    integer, intent(in) :: src_level(:)
    ! Lowest level where wind tendencies are calculated.
    integer, intent(in) :: tend_level(:)
    ! Time step.
    real(kind_phys), intent(in) :: dt
    ! Midpoint temperatures.
    real(kind_phys), intent(in) :: t(:, :)
    ! Log of interface pressures.
    real(kind_phys), intent(in) :: piln(:, :)
    ! Interface densities.
    real(kind_phys), intent(in) :: rhoi(:, :)
    ! Midpoint and interface Brunt-Vaisalla frequencies.
    real(kind_phys), intent(in) :: nm(:, :), ni(:, :)
    ! Projection of wind at midpoints and interfaces.
    real(kind_phys), intent(in) :: ubm(:, :), ubi(:, :)
    ! Unit vectors of source wind (zonal and meridional components).
    real(kind_phys), intent(in) :: xv(:), yv(:)
    ! Tendency efficiency.
    real(kind_phys), intent(in) :: effgw(:)
    ! Wave phase speeds for each column.
    real(kind_phys), intent(in) :: c(:, -band%ngwv:)
    ! Molecular thermal diffusivity.
    real(kind_phys), intent(in) :: kvtt(:, :)
    ! Constituent array.
    real(kind_phys), intent(in) :: q(:, :, :)
    ! Dry static energy.
    real(kind_phys), intent(in) :: dse(:, :)
    ! Coefficient to ramp down diffusion coeff.
    real(kind_phys), intent(in) :: vramp(:)
    ! Wave Reynolds stress.
    real(kind_phys), intent(inout) :: tau(:, -band%ngwv:, :)
    ! Zonal/meridional wind tendencies.
    real(kind_phys), intent(out) :: utgw(:, :), vtgw(:, :)
    ! Gravity wave heating tendency.
    real(kind_phys), intent(out) :: ttgw(:, :)
    ! Gravity wave constituent tendency.
    real(kind_phys), intent(out) :: qtgw(:, :, :)
    ! Effective gravity wave diffusivity at interfaces.
    real(kind_phys), intent(out) :: egwdffi(:, :)
    ! Gravity wave wind tendency for each wave.
    real(kind_phys), intent(out) :: gwut(:, :, -band%ngwv:)
    ! Temperature tendencies from diffusion and kinetic energy.
    real(kind_phys), intent(out) :: dttdf(:, :)
    real(kind_phys), intent(out) :: dttke(:, :)
    ! Adjustment parameter for IGWs.
    real(kind_phys), intent(in), optional :: ro_adjust(:, -band%ngwv:, :)
    ! Diagnosed horizontal wavenumber for ridges.
    real(kind_phys), intent(in), optional :: kwvrdg(:)
    ! Factor for saturation calculation.
    real(kind_phys), intent(in), optional :: satfac_in
    logical, intent(in), optional :: lapply_effgw_in, lapply_vdiff
    ! Provisional Wave Reynolds stress.
    real(kind_phys), intent(out), optional :: tau_diag(:, :)

    !---------------------------Local storage-------------------------------
    integer :: i
    ! Lowest tendency and source levels (over all columns).
    integer :: kbot_tend, kbot_src
    ! Work arrays for gw_drag_prof_column.
    real(kind_phys) :: egwdffm(ncol, pver), dpidz_sq(ncol, pver + 1)
    real(kind_phys) :: ca(ncol, pver), ze(ncol, pver), dnom(ncol, pver)
    real(kind_phys) :: qwork(ncol, pver)

    ! Lowest levels that loops need to iterate over.
    kbot_tend = maxval(tend_level)
    kbot_src = maxval(src_level)

    do i = 1, ncol
      call gw_drag_prof_column(i, band%ngwv, band%kwv, band%effkwv, &
                               kbot_src, kbot_tend, src_level, tend_level, dt, &
                               t, vramp, piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                               effgw, c, kvtt, q, dse, p%del, p%rdel, p%rdst, &
                               tau, utgw, vtgw, ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                               egwdffm, dpidz_sq, ca, ze, dnom, qwork, &
                               ro_adjust=ro_adjust, kwvrdg=kwvrdg, satfac_in=satfac_in, &
                               lapply_effgw_in=lapply_effgw_in, lapply_vdiff=lapply_vdiff, &
                               tau_diag=tau_diag)
    end do

  end subroutine gw_drag_prof

!==========================================================================

  subroutine gw_drag_prof_column(icol, ngwv, kwv, effkwv, kbot_src, kbot_tend, &
                                 src_level, tend_level, dt, &
                                 t, vramp, piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                                 effgw, c, kvtt, q, dse, del, rdel, rdst, &
                                 tau, utgw, vtgw, ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                                 egwdffm, dpidz_sq, ca, ze, dnom, qwork, &
                                 ro_adjust, kwvrdg, satfac_in, lapply_effgw_in, lapply_vdiff, &
                                 tau_diag)
    !$acc routine seq

    !-----------------------------------------------------------------------
    ! Solve for the drag profile from the multiple gravity wave drag
    ! parameterization, for column icol.
    ! 1. scan up from the wave source to determine the stress profile
    ! 2. scan down the stress profile to determine the tendencies
    !     => apply bounds to the tendency
    !          a. from wkb solution
    !          b. from computational stability constraints
    !     => adjust stress on interface below to reflect actual bounded
    !        tendency
    !
    ! Arrays are passed whole, (ncol, ...), and only column icol is read or
    ! written; the arithmetic per column is that of the former array code
    ! (whole-column where blocks written out per point, same statement
    ! order). Level values are kept in scalars and stored once.
    ! kbot_src and kbot_tend are maxval(src_level) and maxval(tend_level)
    ! over all columns, so every column loops over the same levels as before.
    !-----------------------------------------------------------------------

    use gw_diffusion, only: gw_ediff_column, gw_diff_tend_column, gw_diff_tend_column_3d

    !------------------------------Arguments--------------------------------
    ! Column index.
    integer, intent(in) :: icol
    ! Band: number of waves either side of zero, wavenumber, effective
    ! wavenumber (GWBand's ngwv, kwv, effkwv).
    integer, intent(in) :: ngwv
    real(kind_phys), intent(in) :: kwv, effkwv
    ! Lowest source and tendency levels over all columns.
    integer, intent(in) :: kbot_src, kbot_tend
    ! Level from which gravity waves are propagated upward.
    integer, intent(in) :: src_level(:)
    ! Lowest level where wind tendencies are calculated.
    ! Using tend_level > src_level allows the orographic waves to prescribe
    ! wave propagation up to a certain level, but then allow wind tendencies
    ! and adjustments to tau below that level.
    integer, intent(in) :: tend_level(:)
    ! Time step.
    real(kind_phys), intent(in) :: dt
    ! Midpoint temperatures.
    real(kind_phys), intent(in) :: t(:, :)
    ! Coefficient to ramp down diffusion coeff.
    real(kind_phys), intent(in) :: vramp(:)
    ! Log of interface pressures.
    real(kind_phys), intent(in) :: piln(:, :)
    ! Interface densities.
    real(kind_phys), intent(in) :: rhoi(:, :)
    ! Midpoint and interface Brunt-Vaisalla frequencies.
    real(kind_phys), intent(in) :: nm(:, :), ni(:, :)
    ! Projection of wind at midpoints and interfaces.
    real(kind_phys), intent(in) :: ubm(:, :), ubi(:, :)
    ! Unit vectors of source wind (zonal and meridional components).
    real(kind_phys), intent(in) :: xv(:), yv(:)
    ! Tendency efficiency.
    real(kind_phys), intent(in) :: effgw(:)
    ! Wave phase speeds for each column.
    real(kind_phys), intent(in) :: c(:, -ngwv:)
    ! Molecular thermal diffusivity.
    real(kind_phys), intent(in) :: kvtt(:, :)
    ! Constituent array.
    real(kind_phys), intent(in) :: q(:, :, :)
    ! Dry static energy.
    real(kind_phys), intent(in) :: dse(:, :)
    ! Pressure coordinates: layer thickness, its inverse, and the inverse
    ! distance between midpoints (Coords1D's del, rdel, rdst).
    real(kind_phys), intent(in) :: del(:, :), rdel(:, :), rdst(:, :)
    ! Wave Reynolds stress.
    real(kind_phys), intent(inout) :: tau(:, -ngwv:, :)
    ! Zonal/meridional wind tendencies.
    real(kind_phys), intent(inout) :: utgw(:, :), vtgw(:, :)
    ! Gravity wave heating tendency.
    real(kind_phys), intent(inout) :: ttgw(:, :)
    ! Gravity wave constituent tendency.
    real(kind_phys), intent(inout) :: qtgw(:, :, :)
    ! Effective gravity wave diffusivity at interfaces.
    real(kind_phys), intent(inout) :: egwdffi(:, :)
    ! Gravity wave wind tendency for each wave.
    real(kind_phys), intent(inout) :: gwut(:, :, -ngwv:)
    ! Temperature tendencies from diffusion and kinetic energy.
    real(kind_phys), intent(inout) :: dttdf(:, :)
    real(kind_phys), intent(inout) :: dttke(:, :)
    ! Work arrays: (ncol, pver) except dpidz_sq (ncol, pver+1).
    real(kind_phys), intent(inout) :: egwdffm(:, :), dpidz_sq(:, :)
    real(kind_phys), intent(inout) :: ca(:, :), ze(:, :), dnom(:, :)
    real(kind_phys), intent(inout) :: qwork(:, :)
    ! Adjustment parameter for IGWs.
    real(kind_phys), intent(in), optional :: ro_adjust(:, -ngwv:, :)
    ! Diagnosed horizontal wavenumber for ridges.
    real(kind_phys), intent(in), optional :: kwvrdg(:)
    ! Factor for saturation calculation. Here backwards
    ! compatibility. I believe it should be 1.0 (jtb).
    ! Looks like it has been 2.0 for a while in CAM.
    real(kind_phys), intent(in), optional :: satfac_in
    logical, intent(in), optional :: lapply_effgw_in, lapply_vdiff
    ! Provisional Wave Reynolds stress.
    real(kind_phys), intent(inout), optional :: tau_diag(:, :)

    !---------------------------Local storage-------------------------------

    ! Level, wavenumber and constituent loop indices.
    integer :: k, l, m

    ! "Total" and saturation diffusivity.
    real(kind_phys) :: d
    ! Imaginary part of vertical wavenumber.
    real(kind_phys) :: mi
    ! Stress after damping.
    real(kind_phys) :: taudmp
    ! Saturation stress.
    real(kind_phys) :: tausat
    ! (ub-c) and (ub-c)**2
    real(kind_phys) :: ubmc, ubmc2
    ! Temporary ubar tendencies (overall, and at wave l).
    real(kind_phys) :: ubt, ubtl
    real(kind_phys) :: wrk
    ! Ratio used for ubt tndmax limiting.
    real(kind_phys) :: ubt_lim_ratio
    ! Wind tendency of one wave, and of the mean wind projected on the
    ! components, at one level.
    real(kind_phys) :: gwut_lev, utgw_lev, vtgw_lev
    ! Kinetic energy heating at one level, summed over waves.
    real(kind_phys) :: dttke_lev

    ! saturation factor. Defaults to 2.0
    ! unless overidden by satfac_in
    real(kind_phys) :: satfac

    logical :: lapply_effgw, do_vertical_diffusion

    !------------------------------------------------------------------------

    if (present(satfac_in)) then
      satfac = satfac_in
    else
      satfac = 2._kind_phys
    end if

    ! Default behavior is to apply vertical diffusion.
    ! The user has the option to turn off vert diffusion
    do_vertical_diffusion = .true.
    if (present(lapply_vdiff)) then
      do_vertical_diffusion = lapply_vdiff
    end if

    ! Default behavior is to apply effgw and
    ! tendency limiters as designed by Sean
    ! Santos (lapply_effgw=.TRUE.). However,
    ! WACCM non-oro GW need to be retuned before
    ! this can done to them. --jtb 03/02/16
    if (present(lapply_effgw_in)) then
      lapply_effgw = lapply_effgw_in
    else
      lapply_effgw = .TRUE.
    end if

    ! Initialize gravity wave drag tendencies to zero.
    do k = 1, size(utgw, 2)
      utgw(icol, k) = 0._kind_phys
      vtgw(icol, k) = 0._kind_phys
    end do
    do l = -ngwv, ngwv
      do k = 1, size(gwut, 2)
        gwut(icol, k, l) = 0._kind_phys
      end do
    end do
    do k = 1, size(dttke, 2)
      dttke(icol, k) = 0._kind_phys
      ttgw(icol, k) = 0._kind_phys
      dttdf(icol, k) = 0._kind_phys
    end do
    do m = 1, size(qtgw, 3)
      do k = 1, size(qtgw, 2)
        qtgw(icol, k, m) = 0._kind_phys
      end do
    end do

    ! Workaround floating point exception issues on Intel by initializing
    ! everything that's first set in a where block.
    mi = 0._kind_phys
    taudmp = 0._kind_phys
    tausat = 0._kind_phys
    ubmc = 0._kind_phys
    ubmc2 = 0._kind_phys
    wrk = 0._kind_phys

    !------------------------------------------------------------------------
    ! Compute the stress profiles and diffusivities
    !------------------------------------------------------------------------

    ! Loop from bottom to top to get stress profiles.
    ! do k = kbot_src-1, ktop, -1 !++jtb I think this is right
    do k = kbot_src, ktop, -1  !++ but this is in model now

      ! Determine the diffusivity for each column.

      d = dback + kvtt(icol, k)

      do l = -ngwv, ngwv

        ! Determine the absolute value of the saturation stress.
        ! Define critical levels where the sign of (u-c) changes between
        ! interfaces.
        ubmc = ubi(icol, k) - c(icol, l)

        tausat = 0.0_kind_phys

        if (src_level(icol) >= k) then
          ! Test to see if u-c has the same sign here as the level below.
          if (ubmc > 0.0_kind_phys .eqv. ubi(icol, k + 1) > c(icol, l)) then
            if (present(kwvrdg)) then
              tausat = abs(kwvrdg(icol)*rhoi(icol, k)*ubmc**3/ &
                           (satfac*ni(icol, k)))
            else
              tausat = abs(effkwv*rhoi(icol, k)*ubmc**3/ &
                           (satfac*ni(icol, k)))
            end if
          end if

          if (present(ro_adjust)) then
            tausat = tausat*sqrt(ro_adjust(icol, l, k))
          end if

          ! Compute stress for each wave. The stress at this level is the
          ! min of the saturation stress and the stress at the level below
          ! reduced by damping. The sign of the stress must be the same as
          ! at the level below.

          ubmc2 = max(ubmc**2, ubmc2mn)
          if (present(kwvrdg)) then
            mi = ni(icol, k)/(2._kind_phys*kwvrdg(icol)*ubmc2)* &  ! Is this 2._kind_phys related to satfac?
                 (alpha(k) + ni(icol, k)**2/ubmc2*d)
            wrk = -2._kind_phys*mi*rog*t(icol, k)*(piln(icol, k + 1) - piln(icol, k))

            taudmp = tau(icol, l, k + 1)
          else
            mi = ni(icol, k)/(2._kind_phys*kwv*ubmc2)* &
                 (alpha(k) + ni(icol, k)**2/ubmc2*d)
            wrk = -2._kind_phys*mi*rog*t(icol, k)*(piln(icol, k + 1) - piln(icol, k))

            taudmp = tau(icol, l, k + 1)*exp(wrk)
          end if

          ! For some reason, PGI 14.1 loses bit-for-bit reproducibility if
          ! we limit tau, so instead limit the arrays used to set it.
          if (tausat <= taumin) tausat = 0._kind_phys
          if (taudmp <= taumin) taudmp = 0._kind_phys

          tau(icol, l, k) = min(taudmp, tausat)
        end if

      end do
    end do

    ! Force tau at the top of the model to zero, if requested.
    if (tau_0_ubc) then
      do l = -ngwv, ngwv
        tau(icol, l, ktop) = 0._kind_phys
      end do
    end if

    ! Write out pre-adjustment tau profile for diagnostc purposes.
    ! Current implementation only makes sense for orographic waves.
    ! Fix later.
    if (present(tau_diag)) then
      do k = 1, size(tau_diag, 2)
        tau_diag(icol, k) = tau(icol, 0, k)
      end do
    end if

    ! Apply efficiency to completed stress profile.
    if (lapply_effgw) then
      do k = ktop, kbot_tend + 1
        do l = -ngwv, ngwv
          if (k - 1 <= tend_level(icol)) then
            tau(icol, l, k) = tau(icol, l, k)*effgw(icol)
          end if
        end do
      end do
    end if

    !------------------------------------------------------------------------
    ! Compute the tendencies from the stress divergence.
    !------------------------------------------------------------------------

    ! Loop over levels from top to bottom
    do k = ktop, kbot_tend

      ! Accumulate the mean wind tendency over wavenumber.
      ubt = 0.0_kind_phys

      do l = -ngwv, ngwv    ! loop over wave

        ! Determine the wind tendency, including excess stress carried down
        ! from above.
        ubtl = gravit*(tau(icol, l, k + 1) - tau(icol, l, k))*rdel(icol, k)

        ! Apply first tendency limit to maintain numerical stability.
        ! Enforce du/dt < |c-u|/dt  so u-c cannot change sign
        !    (u^n+1 = u^n + du/dt * dt)
        ! The limiter is somewhat stricter, so that we don't come anywhere
        ! near reversing c-u.
        ubtl = min(ubtl, umcfac*abs(c(icol, l) - ubm(icol, k))/dt)

        if (.not. lapply_effgw) ubtl = min(ubtl, tndmax)

        if (k <= tend_level(icol)) then

          ! Save tendency for each wave (for later computation of kzz).
          ! sign function returns magnitude of ubtl with sign of c-ubm
          ! Renders ubt/ubm check for mountain waves unecessary
          gwut_lev = sign(ubtl, c(icol, l) - ubm(icol, k))
          gwut(icol, k, l) = gwut_lev
          ubt = ubt + gwut_lev

        end if

      end do

      if (lapply_effgw) then
        ! Apply second tendency limit to maintain numerical stability.
        ! Enforce du/dt < tndmax so that ridicuously large tendencies are not
        ! permitted.
        ! This can only happen above tend_level, so don't bother checking the
        ! level explicitly.
        if (abs(ubt) > tndmax) then
          ubt_lim_ratio = tndmax/abs(ubt)
          ubt = ubt_lim_ratio*ubt
        else
          ubt_lim_ratio = 1._kind_phys
        end if
      else
        ubt_lim_ratio = 1._kind_phys
      end if

      do l = -ngwv, ngwv
        gwut_lev = ubt_lim_ratio*gwut(icol, k, l)
        ! Redetermine the effective stress on the interface below from the
        ! wind tendency. If the wind tendency was limited above, then the
        ! new stress will be smaller than the old stress, causing stress
        ! divergence in the next layer down. This smoothes large stress
        ! divergences downward while conserving total stress.

        ! Protection on SMALL gwut to prevent floating point
        ! issues.
        !--------------------------------------------------
        if (abs(gwut_lev) < 1.e-15_kind_phys) then
          gwut_lev = 0._kind_phys
        end if
        gwut(icol, k, l) = gwut_lev

        if (k <= tend_level(icol)) then
          tau(icol, l, k + 1) = tau(icol, l, k) + &
                                abs(gwut_lev)*del(icol, k)/gravit
        end if
      end do

      ! Project the mean wind tendency onto the components.
      utgw_lev = utgw(icol, k)
      vtgw_lev = vtgw(icol, k)
      if (k <= tend_level(icol)) then
        utgw_lev = ubt*xv(icol)
        vtgw_lev = ubt*yv(icol)
      end if

      utgw(icol, k) = utgw_lev*vramp(k)
      vtgw(icol, k) = vtgw_lev*vramp(k)

      ! End of level loop.
    end do

    ! Block to undo Sean Santos mods to effgw and limiters.
    ! Here because non-oro GW in WACCM need extensive re-tuning
    ! before Sean's mods can be adopted. --jtb 03/02/16
    !==========================================
    if (.not. (lapply_effgw)) then
      do k = ktop, kbot_tend + 1
        do l = -ngwv, ngwv
          if (k - 1 <= tend_level(icol)) then
            tau(icol, l, k) = tau(icol, l, k)*effgw(icol)
          end if
        end do
      end do
      do k = ktop, kbot_tend
        do l = -ngwv, ngwv
          gwut(icol, k, l) = gwut(icol, k, l)*effgw(icol)
        end do
        utgw(icol, k) = utgw(icol, k)*effgw(icol)
        vtgw(icol, k) = vtgw(icol, k)*effgw(icol)
      end do
    end if
    !===========================================

    if (do_vertical_diffusion) then

      ! Calculate effective diffusivity and LU decomposition for the
      ! vertical diffusion solver.
      call gw_ediff_column(icol, ngwv, kbot_tend, ktop, tend_level, &
                           gwut, ubm, nm, rhoi, dt, prndl, gravit, &
                           del, rdel, rdst, c, vramp, &
                           egwdffi, egwdffm, dpidz_sq, ca, ze, dnom, ro_adjust=ro_adjust)

      ! Calculate tendency on each constituent.
      do m = 1, size(q, 3)

        call gw_diff_tend_column_3d(icol, kbot_tend, ktop, m, q, dt, ca, ze, dnom, &
                                    qtgw, qwork)

      end do

      ! Calculate tendency from diffusing dry static energy (dttdf).
      call gw_diff_tend_column(icol, kbot_tend, ktop, dse, dt, ca, ze, dnom, dttdf, qwork)

    end if

    ! Evaluate second temperature tendency term: Conversion of kinetic
    ! energy into thermal (summed over waves in wave order, as before).
    do k = ktop, kbot_tend
      dttke_lev = dttke(icol, k)
      do l = -ngwv, ngwv
        dttke_lev = dttke_lev - (ubm(icol, k) - c(icol, l))*gwut(icol, k, l)
      end do
      dttke(icol, k) = dttke_lev
    end do

    do k = 1, size(ttgw, 2)
      if (k >= ktop .and. k <= kbot_tend) then
        ttgw(icol, k) = (dttke(icol, k) + dttdf(icol, k))*vramp(k)
      else
        ttgw(icol, k) = dttke(icol, k) + dttdf(icol, k)
      end if
    end do

  end subroutine gw_drag_prof_column

!==========================================================================

! Calculate Reynolds stress for waves propagating in each cardinal
! direction.

  function calc_taucd(ncol, ngwv, tend_level, tau, c, xv, yv, ubi) result(taucd)

    ! Column, level, and gravity wave wavenumber dimensions.
    integer, intent(in) :: ncol, ngwv
    ! Lowest level where wind tendencies are calculated.
    integer, intent(in) :: tend_level(ncol)
    ! Wave Reynolds stress.
    real(kind_phys), intent(in) :: tau(ncol, -ngwv:ngwv, pver+1)
    ! Wave phase speeds for each column.
    real(kind_phys), intent(in) :: c(ncol, -ngwv:ngwv)
    ! Unit vectors of source wind (zonal and meridional components).
    real(kind_phys), intent(in) :: xv(ncol), yv(ncol)
    ! Projection of wind at interfaces.
    real(kind_phys), intent(in) :: ubi(ncol, pver+1)

    real(kind_phys) :: taucd(ncol, pver + 1, 4)

    ! Indices.
    integer :: i, k, l

    ! ubi at tend_level.
    real(kind_phys) :: ubi_tend(ncol)

    ! Signed wave Reynolds stress.
    real(kind_phys) :: tausg(ncol)

    ! Reynolds stress for waves propagating behind and forward of the wind.
    real(kind_phys) :: taub(ncol)
    real(kind_phys) :: tauf(ncol)

    taucd = 0._kind_phys
    tausg = 0._kind_phys

    ubi_tend = [(ubi(i, tend_level(i) + 1), i=1, ncol)]

    do k = ktop, maxval(tend_level) + 1

      taub = 0._kind_phys
      tauf = 0._kind_phys

      do l = -ngwv, ngwv
        where (k - 1 <= tend_level)

          tausg = sign(tau(:, l, k), c(:, l) - ubi(:, k))

          where (c(:, l) < ubi_tend)
            taub = taub + tausg
          elsewhere
            tauf = tauf + tausg
          end where

        end where
      end do

      where (k - 1 <= tend_level)
        where (xv > 0._kind_phys)
          taucd(:, k, east) = tauf*xv
          taucd(:, k, west) = taub*xv
        elsewhere
          taucd(:, k, east) = taub*xv
          taucd(:, k, west) = tauf*xv
        end where

        where (yv > 0._kind_phys)
          taucd(:, k, north) = tauf*yv
          taucd(:, k, south) = taub*yv
        elsewhere
          taucd(:, k, north) = taub*yv
          taucd(:, k, south) = tauf*yv
        end where
      end where

    end do

  end function calc_taucd

!==========================================================================

! Calculate the amount of momentum conveyed from below the gravity wave
! region, to the region where gravity waves are calculated.
  subroutine momentum_flux(tend_level, taucd, um_flux, vm_flux)

    ! Bottom stress level.
    integer, intent(in) :: tend_level(:)
    ! Projected stresses.
    real(kind_phys), intent(in) :: taucd(:, :, :)
    ! Components of momentum change sourced from the bottom.
    real(kind_phys), intent(out) :: um_flux(:), vm_flux(:)

    integer :: i

    ! Tendency for U & V below source level.
    do i = 1, size(tend_level)
      um_flux(i) = taucd(i, tend_level(i) + 1, east) + &
                   taucd(i, tend_level(i) + 1, west)
      vm_flux(i) = taucd(i, tend_level(i) + 1, north) + &
                   taucd(i, tend_level(i) + 1, south)
    end do

  end subroutine momentum_flux

  ! Subtracts a change in momentum in the gravity wave levels from wind
  ! tendencies in lower levels, ensuring momentum conservation.
  subroutine momentum_fixer(tend_level, p, um_flux, vm_flux, utgw, vtgw)

    use coords_1d, only: coords1d

    ! Bottom stress level.
    integer, intent(in) :: tend_level(:)
    ! Pressure coordinates.
    type(coords1d), intent(in) :: p
    ! Components of momentum change sourced from the bottom.
    real(kind_phys), intent(in) :: um_flux(:), vm_flux(:)
    ! Wind tendencies.
    real(kind_phys), intent(inout) :: utgw(:, :), vtgw(:, :)

    ! Indices.
    integer :: i, k
    ! Reciprocal of total mass.
    real(kind_phys) :: rdm(size(tend_level))
    ! Average changes in velocity from momentum change being spread over
    ! total mass.
    real(kind_phys) :: du(size(tend_level)), dv(size(tend_level))

    ! Total mass from ground to source level: rho*dz = dp/gravit
    do i = 1, size(tend_level)
      rdm(i) = gravit/(p%ifc(i, pver + 1) - p%ifc(i, tend_level(i) + 1))
    end do

    ! Average velocity changes.
    du = -um_flux*rdm
    dv = -vm_flux*rdm

    do k = minval(tend_level) + 1, pver
      where (k > tend_level)
        utgw(:, k) = utgw(:, k) + du
        vtgw(:, k) = vtgw(:, k) + dv
      end where
    end do

  end subroutine momentum_fixer

  ! Calculate the change in total energy from tendencies up to this point,
  ! for all columns (host loop over energy_change_column).
  subroutine energy_change(dt, p, u, v, dudt, dvdt, dsdt, de)

    use coords_1d, only: coords1d

    ! Time step.
    real(kind_phys), intent(in) :: dt
    ! Pressure coordinates.
    type(coords1d), intent(in) :: p
    ! Winds at start of time step.
    real(kind_phys), intent(in) :: u(:, :), v(:, :)
    ! Wind tendencies.
    real(kind_phys), intent(in) :: dudt(:, :), dvdt(:, :)
    ! Heating tendency.
    real(kind_phys), intent(in) :: dsdt(:, :)
    ! Change in energy.
    real(kind_phys), intent(out) :: de(:)

    ! Column index.
    integer :: i

    do i = 1, size(de)
      call energy_change_column(i, pver, dt, p%del, u, v, dudt, dvdt, dsdt, de)
    end do

  end subroutine energy_change

  ! Calculate the change in total energy from tendencies up to this point,
  ! for column icol (levels summed in order).
  subroutine energy_change_column(icol, nlev, dt, del, u, v, dudt, dvdt, dsdt, de)
    !$acc routine seq

    ! Column index and number of levels.
    integer, intent(in) :: icol, nlev
    ! Time step.
    real(kind_phys), intent(in) :: dt
    ! Layer thickness (Coords1D's del) [Pa].
    real(kind_phys), intent(in) :: del(:, :)
    ! Winds at start of time step.
    real(kind_phys), intent(in) :: u(:, :), v(:, :)
    ! Wind tendencies.
    real(kind_phys), intent(in) :: dudt(:, :), dvdt(:, :)
    ! Heating tendency.
    real(kind_phys), intent(in) :: dsdt(:, :)
    ! Change in energy.
    real(kind_phys), intent(inout) :: de(:)

    ! Level index.
    integer :: k
    ! Net gain/loss of total energy in the column.
    real(kind_phys) :: de_col

    de_col = 0.0_kind_phys
    do k = 1, nlev
      de_col = de_col + del(icol, k)/gravit*(dsdt(icol, k) + &
                                             dudt(icol, k)*(u(icol, k) + dudt(icol, k)*0.5_kind_phys*dt) + &
                                             dvdt(icol, k)*(v(icol, k) + dvdt(icol, k)*0.5_kind_phys*dt))
    end do
    de(icol) = de_col

  end subroutine energy_change_column

  ! Subtract change in energy from the heating tendency in the levels below
  ! the gravity wave region.
  subroutine energy_fixer(tend_level, p, de, ttgw)

    use coords_1d, only: coords1d

    ! Bottom stress level.
    integer, intent(in) :: tend_level(:)
    ! Pressure coordinates.
    type(coords1d), intent(in) :: p
    ! Change in energy.
    real(kind_phys), intent(in) :: de(:)
    ! Heating tendency.
    real(kind_phys), intent(inout) :: ttgw(:, :)

    ! Column/level indices.
    integer :: i, k
    ! Energy change to apply divided by all the mass it is spread across.
    real(kind_phys) :: de_dm(size(tend_level))

    do i = 1, size(tend_level)
      de_dm(i) = -de(i)*gravit/(p%ifc(i, pver + 1) - p%ifc(i, tend_level(i) + 1))
    end do

    ! Subtract net gain/loss of total energy below tend_level.
    do k = minval(tend_level) + 1, pver
      where (k > tend_level)
        ttgw(:, k) = ttgw(:, k) + de_dm
      end where
    end do

  end subroutine energy_fixer

!==========================================================================

  ! Calculates absolute value of the local Coriolis frequency divided by the
  ! spatial frequency kwv, which gives a characteristic speed in m/s.
  pure function coriolis_speed(band, lat)
    ! Inertial gravity wave lengths.
    type(GWBand), intent(in) :: band
    ! Latitude in radians.
    real(kind_phys), intent(in) :: lat(:)

    real(kind_phys) :: coriolis_speed(size(lat))

    coriolis_speed = abs(sin(lat)*2._kind_phys*omega_earth/band%kwv)

  end function coriolis_speed

!==========================================================================

  subroutine adjust_inertial(band, tend_level, &
                             u_coriolis, c, ubi, tau, ro_adjust)
    ! Inertial gravity wave lengths.
    type(GWBand), intent(in) :: band
    ! Levels above which tau is calculated.
    integer, intent(in) :: tend_level(:)
    ! Absolute value of the Coriolis frequency for each column,
    ! divided by kwv [m/s].
    real(kind_phys), intent(in) :: u_coriolis(:)
    ! Wave propagation speed.
    real(kind_phys), intent(in) :: c(:, -band%ngwv:)
    ! Wind speed in the direction of wave propagation.
    real(kind_phys), intent(in) :: ubi(:, :)

    ! Tau will be adjusted by blocking wave propagation through cells where
    ! the Coriolis effect prevents it.
    real(kind_phys), intent(inout) :: tau(:, -band%ngwv:, :)
    ! Dimensionless Coriolis term used to reduce gravity wave strength.
    ! Equal to max(0, 1 - (1/ro)^2), where ro is the Rossby number of the
    ! wind with respect to inertial waves.
    real(kind_phys), intent(out) :: ro_adjust(:, -band%ngwv:, :)

    ! Column/level/wavenumber indices.
    integer :: i, k, l

    ! For each column and wavenumber, are we clear of levels that block
    ! upward propagation?
    logical :: unblocked_mask(size(tend_level), -band%ngwv:band%ngwv)

    unblocked_mask = .true.
    ro_adjust = 0._kind_phys

    ! Iterate from the bottom up, through every interface level where tau is
    ! set.
    do k = maxval(tend_level) + 1, ktop, -1
      do l = -band%ngwv, band%ngwv
        do i = 1, size(tend_level)
          ! Only operate on valid levels for this column.
          if (k <= tend_level(i) + 1) then
            ! Block waves if Coriolis is too strong.
            ! By setting the mask in this way, we avoid division by zero.
            unblocked_mask(i, l) = unblocked_mask(i, l) .and. &
                                   (abs(ubi(i, k) - c(i, l)) > u_coriolis(i))
            if (unblocked_mask(i, l)) then
              ro_adjust(i, l, k) = &
                1._kind_phys - (u_coriolis(i)/(ubi(i, k) - c(i, l)))**2
            else
              tau(i, l, k) = 0._kind_phys
            end if
          end if
        end do
      end do
    end do

  end subroutine adjust_inertial

  ! Given a value, finds which bin marked by "bounds" the value falls
  ! into.
  pure elemental function find_bin(val) result(idx)
    real(kind_phys), intent(in) :: val

    integer :: idx

    ! We just have to count how many bounds are exceeded.
    if (val >= 0._kind_phys) then
      idx = count(val > bounds) + 1
    else
      idx = count(val >= bounds) + 1
    end if

  end function find_bin

end module gw_common
