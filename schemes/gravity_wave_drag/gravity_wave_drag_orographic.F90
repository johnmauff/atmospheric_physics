! Stationary gravity waves from orographic sources
! Extracted: May 2013
! CCPPized:  August 2025
module gravity_wave_drag_orographic

  use ccpp_kinds, only: kind_phys
  use gw_common, only: GWBand

  implicit none
  private

  public :: gravity_wave_drag_orographic_init
  public :: gravity_wave_drag_orographic_run

  type(GWBand) :: band_oro

  ! Limiters (min/max values)
  ! min surface displacement height for orographic waves [m]
  real(kind_phys), parameter :: orohmin = 10._kind_phys
  ! min wind speed for orographic waves [m s-1]
  real(kind_phys), parameter :: orovmin = 2._kind_phys

contains

!> \section arg_table_gravity_wave_drag_orographic_init Argument Table
!! \htmlinclude arg_table_gravity_wave_drag_orographic_init.html
  subroutine gravity_wave_drag_orographic_init(&
             gw_delta_c, &
             effgw_oro, &
             errmsg, errflg)

    use gw_common, only: wavelength_mid, unset_kind_phys

    real(kind_phys),    intent(in)                :: gw_delta_c
    real(kind_phys),    intent(in)                :: effgw_oro
    character(len=512), intent(out)               :: errmsg
    integer, intent(out)                          :: errflg

    errmsg = ''
    errflg = 0

    if(effgw_oro == unset_kind_phys) then
      errflg = 1
      errmsg = "gravity_wave_drag_orographic_init: orographic gravity wave drag is on but effgw_oro is unset."
      return
    endif

    ! Band to emit orographic waves in.
    ! Regardless, we will only ever emit into l = 0.
    band_oro    = GWBand(0,         gw_delta_c,      1.0_kind_phys, wavelength_mid)

  end subroutine gravity_wave_drag_orographic_init

  ! Orographic source for multiple gravity wave drag parameterization.
  !
  ! The stress is returned for a single wave with c=0, over orography.
  ! For points where the orographic variance is small (including ocean),
  ! the returned stress is zero.
!> \section arg_table_gravity_wave_drag_orographic_run Argument Table
!! \htmlinclude arg_table_gravity_wave_drag_orographic_run.html
  subroutine gravity_wave_drag_orographic_run( &
             ncol, pver, pcnst, &
             dt, &
             rair, &
             cpairv, &
             pint, vramp, &
             piln, rhoi, nm, ni, &
             effgw_oro, gw_lndscl_sgh, gw_oro_south_fac, gw_apply_tndmax, &
             landfrac, &
             lat, &
             u, v, t, q, dse, &
             sgh, zm, &
             kvt_gw, &
             tend_q, tend_u, tend_v, tend_s, &
             src_level, tend_level, ubm, ubi, xv, yv, &
             utgw, vtgw, ttgw, qtgw, &
             dttdf, dttke, egwdffi_tot, &
             flx_heat, &
             tau0x, tau0y, taua, &
             errmsg, errflg)

    use gw_common, only: gw_drag_prof_column, energy_change_column
    use atmos_phys_fin_vol_solver, only: fin_vol_coords_column

    integer,            intent(in)                :: ncol
    integer,            intent(in)                :: pver
    integer,            intent(in)                :: pcnst
    real(kind_phys),    intent(in)                :: dt                       ! Physics timestep [s]
    real(kind_phys),    intent(in)                :: rair                     ! Gas constant for dry air [J kg-1 K-1]
    real(kind_phys),    intent(in)                :: cpairv(:,:)              ! Specific heat of dry air at constant pressure [J kg-1 K-1]
    real(kind_phys),    intent(in)                :: pint(:, :)               ! Air pressure at interfaces [Pa]
    real(kind_phys),    intent(in)                :: vramp(:)                 ! Gravity wave drag tapering coefficients [1]
    real(kind_phys),    intent(in)                :: piln(:, :)               ! Natural log of pressure at interfaces [ln(Pa)]
    real(kind_phys),    intent(in)                :: rhoi(:, :)               ! Density at interfaces [kg m-3]
    real(kind_phys),    intent(in)                :: nm(:, :)                 ! Brunt-Vaisalla frequency at midpoints [s-1]
    real(kind_phys),    intent(in)                :: ni(:, :)                 ! Brunt-Vaisalla frequency at interfaces [s-1]
    real(kind_phys),    intent(in)                :: effgw_oro                ! Efficiency factor for orographic gravity waves [1]
    logical,            intent(in)                :: gw_lndscl_sgh            ! Scale SGH by land fraction? [flag]
    real(kind_phys),    intent(in)                :: gw_oro_south_fac         ! Southern hemisphere scaling factor for orographic waves [1]
    logical,            intent(in)                :: gw_apply_tndmax          ! Maximum wind tendency limiter from stress divergence [flag]
    real(kind_phys),    intent(in)                :: landfrac(:)              ! Land fraction [fraction]
    real(kind_phys),    intent(in)                :: lat(:)                   ! Latitude [rad]
    real(kind_phys),    intent(in)                :: u(:,:)                   ! Zonal wind at midpoints [m s-1]
    real(kind_phys),    intent(in)                :: v(:,:)                   ! Meridional wind at midpoints [m s-1]
    real(kind_phys),    intent(in)                :: t(:,:)                   ! Temperature at midpoints [K]
    real(kind_phys),    intent(in)                :: q(:,:,:)                 ! Constituent mixing ratios [kg kg-1]
    real(kind_phys),    intent(in)                :: dse(:,:)                 ! Dry static energy [J kg-1]
    real(kind_phys),    intent(in)                :: sgh(:)                   ! Standard deviation of orography [m]
    real(kind_phys),    intent(in)                :: zm(:,:)                  ! Geopotential height at midpoints [m]
    real(kind_phys),    intent(in)                :: kvt_gw(:,:)              ! Eddy diffusion coefficient for heat [m2 s-1]

    real(kind_phys),    intent(inout)             :: tend_q(:, :, :)          ! Constituent tendencies [kg kg-1 s-1]
    real(kind_phys),    intent(inout)             :: tend_u(:, :)             ! Zonal wind tendency [m s-2]
    real(kind_phys),    intent(inout)             :: tend_v(:, :)             ! Meridional wind tendency [m s-2]
    real(kind_phys),    intent(inout)             :: tend_s(:, :)             ! Dry static energy tendency [J kg-1 s-1]

    integer,            intent(out)               :: src_level(:)             ! Vertical level index of gravity wave source [index]
    integer,            intent(out)               :: tend_level(:)            ! Lowest vertical level index where tendencies are applied [index]
    real(kind_phys),    intent(out)               :: ubm(:, :)                ! Wind projection at midpoints along source wind direction [m s-1]
    real(kind_phys),    intent(out)               :: ubi(:, :)                ! Wind projection at interfaces along source wind direction [m s-1]
    real(kind_phys),    intent(out)               :: xv(:)                    ! Zonal component of source wind unit vector [1]
    real(kind_phys),    intent(out)               :: yv(:)                    ! Meridional component of source wind unit vector [1]
    real(kind_phys),    intent(out)               :: utgw(:, :)               ! Zonal wind tendency from gravity waves [m s-2]
    real(kind_phys),    intent(out)               :: vtgw(:, :)               ! Meridional wind tendency from gravity waves [m s-2]
    real(kind_phys),    intent(out)               :: ttgw(:, :)               ! Temperature tendency from gravity waves [K s-1]
    real(kind_phys),    intent(out)               :: qtgw(:, :, :)            ! Constituent tendencies from gravity waves [kg kg-1 s-1]
    real(kind_phys),    intent(out)               :: dttdf(:, :)              ! Dry air enthalpy tendency from diffusion [J kg-1 s-1]
    real(kind_phys),    intent(out)               :: dttke(:, :)              ! Dry air enthalpy tendency from kinetic energy dissipation [J kg-1 s-1]
    real(kind_phys),    intent(inout)             :: egwdffi_tot(:, :)        ! Effective diffusivity coefficient from gravity waves, interfaces [m2 s-1]
    real(kind_phys),    intent(out)               :: flx_heat(:)              ! Surface heat flux for energy conservation check [W m-2]
    real(kind_phys),    intent(out)               :: tau0x(:)                 ! Zonal gravity wave surface stress [N m-2]
    real(kind_phys),    intent(out)               :: tau0y(:)                 ! Meridional gravity wave surface stress [N m-2]
    real(kind_phys),    intent(out)               :: taua(:,:)                ! Total stress from orographic scheme, interface [N m-2]
    character(len=512), intent(out)               :: errmsg
    integer,            intent(out)               :: errflg

    integer :: i, k, l, m

    ! Wavenumber fields
    real(kind_phys) :: tau(ncol, -band_oro%ngwv:band_oro%ngwv, pver+1)
    real(kind_phys) :: gwut(ncol, pver, -band_oro%ngwv:band_oro%ngwv)
    real(kind_phys) :: phase_speeds(ncol, -band_oro%ngwv:band_oro%ngwv)

    real(kind_phys) :: effgw(ncol)

    real(kind_phys) :: egwdffi(ncol, pver+1)

    ! Surface streamline displacement height (2*sgh).
    real(kind_phys) :: hdsp(ncol)

    ! Energy change used by fixer.
    real(kind_phys) :: de(ncol)

    ! Pressure coordinates per column (Coords1D(pint)'s mid, del, rdel, rdst).
    real(kind_phys) :: mid(ncol, pver), del(ncol, pver), rdel(ncol, pver)
    real(kind_phys) :: rdst(ncol, pver-1)

    ! Whether any column extends its source region at level k.
    logical :: anyc(pver)
    logical :: found

    ! Work arrays for gw_drag_prof_column.
    real(kind_phys) :: egwdffm(ncol, pver), dpidz_sq(ncol, pver+1)
    real(kind_phys) :: ca(ncol, pver), ze(ncol, pver), dnom(ncol, pver)
    real(kind_phys) :: qwork(ncol, pver)

    ! Band scalars (band_oro stays on the host).
    integer         :: ngwv
    real(kind_phys) :: kwv, effkwv, fcrit2

    ! Lowest source and tendency levels over all columns.
    integer :: kbot_src, kbot_tend

    ! Per-column scalars.
    integer         :: srcl                  ! Source level.
    real(kind_phys) :: sgh_local             ! Scaled sgh when gw_lndscl_sgh.
    real(kind_phys) :: sghmax                ! Max orographic standard deviation to use.
    real(kind_phys) :: tauoro                ! c=0 stress from orography.
    real(kind_phys) :: tauzero               ! Zero stress (scaled like tau).
    real(kind_phys) :: nsrc, rsrc, usrc, vsrc  ! Averages over source region: B-V frequency, density, winds.
    real(kind_phys) :: dpsrc                 ! Difference in interface pressure across source region.
    real(kind_phys) :: mag, xvi, yvi         ! Source wind magnitude and unit vector.
    real(kind_phys) :: ubm_k, ubm_km1        ! Wind projection at midpoints k, k-1.
    real(kind_phys) :: tg                    ! Heating tendency at one point.

    errmsg = ''
    errflg = 0

    ngwv   = band_oro%ngwv
    kwv    = band_oro%kwv
    effkwv = band_oro%effkwv
    fcrit2 = band_oro%fcrit2

    ! OpenACC: one data region; per-column kernels (one thread per column,
    ! levels in order), plus one per-level reduction for the source search.
    ! The arithmetic per column is that of the former array code; level
    ! values are kept in scalars and stored once.
    !$acc data copyin(pint, cpairv, vramp, piln, rhoi, nm, ni, landfrac, lat, &
    !$acc             u, v, t, q, dse, sgh, zm, kvt_gw) &
    !$acc      copy(tend_q, tend_u, tend_v, tend_s, egwdffi_tot) &
    !$acc      copyout(src_level, tend_level, ubm, ubi, xv, yv, utgw, vtgw, ttgw, qtgw, &
    !$acc              dttdf, dttke, flx_heat, tau0x, tau0y, taua) &
    !$acc      create(tau, gwut, phase_speeds, effgw, egwdffi, hdsp, de, mid, del, rdel, &
    !$acc             rdst, anyc, egwdffm, dpidz_sq, ca, ze, dnom, qwork)

    ! Scale sgh by land fraction? Surface streamline displacement height,
    ! and the pressure coordinates.
    !$acc parallel loop gang vector private(sgh_local)
    do i = 1, ncol
      if (gw_lndscl_sgh) then
        if (landfrac(i) >= epsilon(1._kind_phys)) then
          effgw(i) = effgw_oro * landfrac(i)
          sgh_local = sgh(i) / sqrt(landfrac(i))
        else
          effgw(i) = 0._kind_phys
          sgh_local = 0._kind_phys
        end if
      else
        effgw(i) = effgw_oro
        sgh_local = sgh(i)
      endif
      hdsp(i) = 2.0_kind_phys*sgh_local
      call fin_vol_coords_column(i, pver, pint, mid, del, rdel, rdst)
    end do

    !--------------------------------------------------------------------------
    ! Average the basic state variables for the wave source over the depth of
    ! the orographic standard deviation. Here we assume that the appropiate
    ! values of wind, stability, etc. for determining the wave source are
    ! averages over the depth of the atmosphere penterated by the typical
    ! mountain.
    ! Reduces to the bottom midpoint values when sgh=0, such as over ocean.
    !
    ! The search walks up from the bottom and stops at the first level where
    ! no column extends its source region (as the former loop over all
    ! columns did): anyc(k) records, per level, whether any column does.
    !--------------------------------------------------------------------------
    !$acc parallel loop gang private(found)
    do k = 1, pver - 1
      found = .false.
      !$acc loop vector reduction(.or.:found)
      do i = 1, ncol
        found = found .or. (hdsp(i) > sqrt(zm(i, k)*zm(i, k + 1)))
      end do
      anyc(k) = found
    end do

    kbot_src = -huge(1)
    !$acc parallel loop gang vector reduction(max:kbot_src) &
    !$acc   private(k, l, srcl, rsrc, usrc, vsrc, nsrc, dpsrc, mag, xvi, yvi, &
    !$acc           ubm_k, ubm_km1, sghmax, tauoro, tauzero)
    do i = 1, ncol
      k = pver
      srcl = k - 1
      rsrc = mid(i, k)/(rair*t(i, k))*del(i, k)
      usrc = u(i, k)*del(i, k)
      vsrc = v(i, k)*del(i, k)
      nsrc = nm(i, k)*del(i, k)

      !$acc loop seq
      do k = pver - 1, 1, -1
        if (hdsp(i) > sqrt(zm(i, k)*zm(i, k + 1))) then
          srcl = k - 1
          rsrc = rsrc + &
                 mid(i, k)/(rair*t(i, k))*del(i, k)
          usrc = usrc + u(i, k)*del(i, k)
          vsrc = vsrc + v(i, k)*del(i, k)
          nsrc = nsrc + nm(i, k)*del(i, k)
        end if
        ! Break the loop when all source levels found.
        if (.not. anyc(k)) exit
      end do

      dpsrc = pint(i, pver + 1) - pint(i, srcl + 1)

      rsrc = rsrc/dpsrc
      usrc = usrc/dpsrc
      vsrc = vsrc/dpsrc
      nsrc = nsrc/dpsrc

      ! Get the unit vector components and magnitude at the surface
      ! (gw_utils' get_unit_vector).
      mag = sqrt(usrc*usrc + vsrc*vsrc)
      if (mag > 0._kind_phys) then
        xvi = usrc/mag
        yvi = vsrc/mag
      else
        xvi = 0._kind_phys
        yvi = 0._kind_phys
      end if
      xv(i) = xvi
      yv(i) = yvi
      ubi(i, pver + 1) = mag

      ! Project the local wind at midpoints onto the source wind (dot_2d).
      ! Compute the interface wind projection by averaging the midpoint winds
      ! (midpoint_interp). Use the top level wind at the top interface.
      ubm_k = u(i, 1)*xvi + v(i, 1)*yvi
      ubm(i, 1) = ubm_k
      ubi(i, 1) = ubm_k
      !$acc loop seq
      do k = 2, pver
        ubm_km1 = ubm_k
        ubm_k = u(i, k)*xvi + v(i, k)*yvi
        ubm(i, k) = ubm_k
        ubi(i, k) = 0.5_kind_phys*(ubm_km1 + ubm_k)
      end do

      ! Determine the orographic c=0 source term following McFarlane (1987).
      ! Set the source top interface index to pver, if the orographic term is
      ! zero.
      if ((mag > orovmin) .and. (hdsp(i) > orohmin)) then
        sghmax = fcrit2*(mag/nsrc)**2
        tauoro = 0.5_kind_phys*kwv*min(hdsp(i)**2, sghmax)* &
                 rsrc*nsrc*mag
      else
        tauoro = 0._kind_phys
        srcl = pver
      end if
      src_level(i) = srcl

      ! Set the phase speeds and wave numbers in the direction of the source
      ! wind. Set the source stress magnitude (positive only, note that the
      ! sign of the stress is the same as (c-u). Scale tau in the southern
      ! hemisphere.
      tauzero = 0._kind_phys
      if (lat(i) < 0._kind_phys) then
        tauzero = tauzero*gw_oro_south_fac
        tauoro = tauoro*gw_oro_south_fac
      end if
      !$acc loop seq
      do k = 1, pver + 1
        !$acc loop seq
        do l = -ngwv, ngwv
          if (l == 0 .and. k >= srcl + 1) then
            tau(i, l, k) = tauoro
          else
            tau(i, l, k) = tauzero
          end if
        end do
      end do

      ! Allow wind tendencies all the way to the model bottom.
      tend_level(i) = pver

      ! No spectrum; phase speed is just 0.
      !$acc loop seq
      do l = -ngwv, ngwv
        phase_speeds(i, l) = 0._kind_phys
      end do

      kbot_src = max(kbot_src, srcl)
    end do

    ! Lowest tendency level over all columns (tend_level = pver everywhere).
    kbot_tend = pver

    !$acc parallel loop gang vector private(k, m, tg)
    do i = 1, ncol
      ! Solve for the drag profile with orographic sources.
      call gw_drag_prof_column(i, ngwv, kwv, effkwv, kbot_src, kbot_tend, &
                               src_level, tend_level, dt, &
                               t, vramp, piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
                               effgw, phase_speeds, kvt_gw, q, dse, del, rdel, rdst, &
                               tau, utgw, vtgw, ttgw, qtgw, egwdffi, gwut, dttdf, dttke, &
                               egwdffm, dpidz_sq, ca, ze, dnom, qwork, &
                               lapply_effgw_in=gw_apply_tndmax)

      ! For orographic waves, don't bother with taucd, since there are no
      ! momentum conservation routines or directional diagnostics.

      ! Add the diffusion coefficients
      !$acc loop seq
      do k = 1, pver + 1
        egwdffi_tot(i, k) = egwdffi_tot(i, k) + egwdffi(i, k)
      end do

      ! Add the orographic tendencies to the spectrum tendencies.
      ! Don't calculate fixers, since we are too close to the ground to
      ! spread momentum/energy differences across low layers.
      !$acc loop seq
      do k = 1, pver
        tend_u(i, k) = tend_u(i, k) + utgw(i, k)
        tend_v(i, k) = tend_v(i, k) + vtgw(i, k)
        tg = ttgw(i, k)
        tend_s(i, k) = tend_s(i, k) + tg
        ! Convert to temperature tendency for output.
        ttgw(i, k) = tg/cpairv(i, k)
      end do

      ! Calculate energy change for output to CAM's energy checker.
      ! This is sort of cheating; we don't have a good a priori idea of the
      ! energy coming from surface stress, so we just integrate what we and
      ! actually have so far and OVERWRITE flx_heat with that.
      call energy_change_column(i, pver, dt, del, u, v, tend_u, tend_v, tend_s, de)
      flx_heat(i) = de(i)            ! this is the only gw scheme that has flx_heat as out and overwritten directly.

      ! Store constituents tendencies
      !$acc loop seq
      do m = 1, pcnst
        !$acc loop seq
        do k = 1, pver
           tend_q(i,k,m) = tend_q(i,k,m) + qtgw(i,k,m)
        end do
      end do

      ! Compute surface stresses for diagnostic
      tau0x(i) = tau(i,0,pver+1) * xv(i)
      tau0y(i) = tau(i,0,pver+1) * yv(i)
      !$acc loop seq
      do k = 1, pver + 1
        taua(i, k) = tau(i,0,k)
      end do
    end do

    !$acc end data

  end subroutine gravity_wave_drag_orographic_run

end module gravity_wave_drag_orographic
