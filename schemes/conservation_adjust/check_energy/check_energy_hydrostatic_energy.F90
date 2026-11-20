! Vertically integrated hydrostatic total energy and water of one column, for
! the energy and water conservation checks.
!
! This is CAM's cam_thermo get_hydrostatic_energy (moist mixing ratios only,
! the form check_energy_chng uses) on whole arrays and a column index, with
! the host-model data it read from module variables (water vapor and the
! thermodynamically active liquid and ice species, the enthalpy reference
! state) passed in as arguments, and errors returned instead of endrun. The
! arithmetic per column is that of the original, operation by operation. It is
! an OpenACC routine seq so it can run in a column kernel.
module check_energy_hydrostatic_energy

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: hydrostatic_energy_column
  public :: enthalpy_reference_code

  ! Total energy formulas: the values of the host's total_energy_formula_for_
  ! physics / _for_dycore (CAM-SIMA cam_thermo_formula ENERGY_FORMULA_DYCORE_*).
  integer, public, parameter :: ENERGY_FORMULA_DYCORE_FV   = 0  ! vc_moist_pressure
  integer, public, parameter :: ENERGY_FORMULA_DYCORE_SE   = 1  ! vc_dry_pressure
  integer, public, parameter :: ENERGY_FORMULA_DYCORE_MPAS = 2  ! vc_height

  ! Enthalpy reference states (enthalpy_reference_code of 'ice', 'liq', 'wv').
  integer, public, parameter :: ENTHALPY_REFERENCE_ICE = 1
  integer, public, parameter :: ENTHALPY_REFERENCE_LIQ = 2
  integer, public, parameter :: ENTHALPY_REFERENCE_WV  = 3

  ! Error codes of hydrostatic_energy_column.
  integer, public, parameter :: HYDROSTATIC_ENERGY_BAD_FORMULA   = 1
  integer, public, parameter :: HYDROSTATIC_ENERGY_BAD_REFERENCE = 2

contains

  ! Integer code of an enthalpy reference state name ('ice', 'liq', 'wv'),
  ! 0 if the name is not supported.
  pure integer function enthalpy_reference_code(state)
    character(len=*), intent(in) :: state

    select case (trim(state))
    case ('ice')
      enthalpy_reference_code = ENTHALPY_REFERENCE_ICE
    case ('liq')
      enthalpy_reference_code = ENTHALPY_REFERENCE_LIQ
    case ('wv')
      enthalpy_reference_code = ENTHALPY_REFERENCE_WV
    case default
      enthalpy_reference_code = 0
    end select
  end function enthalpy_reference_code

  ! Vertical integrals of column icol (moist mixing ratios: pdel is the moist
  ! layer thickness). Outputs are written at index icol of the optional
  ! arrays given. errflg is 0, HYDROSTATIC_ENERGY_BAD_FORMULA (vcoord not FV,
  ! SE or MPAS) or HYDROSTATIC_ENERGY_BAD_REFERENCE (enthalpy_reference not
  ! one of the ENTHALPY_REFERENCE_* codes; checked only when te is wanted).
  subroutine hydrostatic_energy_column(icol, tracer, pdel, cp_or_cv, U, V, T, &
       vcoord, ptop, phis, wv_idx, liq_idx, ice_idx, enthalpy_reference, &
       rga, latvap, latice, errflg, z_mid, te, se, po, ke, wv, H2O, liq, ice)
    !$acc routine seq

    ! Column index.
    integer,         intent(in)  :: icol
    ! Moist constituent mixing ratios (ncol, nlev, number of constituents).
    real(kind_phys), intent(in)  :: tracer(:,:,:)
    ! Moist pressure level thickness (ncol, nlev) [Pa].
    real(kind_phys), intent(in)  :: pdel(:,:)
    ! Air heat capacity at constant pressure or volume (energy formula
    ! dependent) (ncol, nlev) [J kg-1 K-1].
    real(kind_phys), intent(in)  :: cp_or_cv(:,:)
    real(kind_phys), intent(in)  :: U(:,:)       ! zonal wind [m s-1]
    real(kind_phys), intent(in)  :: V(:,:)       ! meridional wind [m s-1]
    real(kind_phys), intent(in)  :: T(:,:)       ! temperature [K]
    ! Energy formula (ENERGY_FORMULA_DYCORE_*).
    integer,         intent(in)  :: vcoord
    real(kind_phys), intent(in)  :: ptop(:,:)    ! interface pressure; ptop(icol, 1) is used [Pa]
    real(kind_phys), intent(in)  :: phis(:)      ! surface geopotential [m2 s-2]
    ! Constituent indices: water vapor, liquid and ice water species.
    integer,         intent(in)  :: wv_idx
    integer,         intent(in)  :: liq_idx(:)
    integer,         intent(in)  :: ice_idx(:)
    ! Enthalpy reference state (ENTHALPY_REFERENCE_*).
    integer,         intent(in)  :: enthalpy_reference
    real(kind_phys), intent(in)  :: rga          ! 1/gravit [s2 m-1]
    real(kind_phys), intent(in)  :: latvap       ! latent heat of vaporization [J kg-1]
    real(kind_phys), intent(in)  :: latice       ! latent heat of fusion [J kg-1]
    integer,         intent(out) :: errflg
    ! Height above ground at midpoints (ncol, nlev), MPAS formula only [m].
    real(kind_phys), intent(in),    optional :: z_mid(:,:)
    ! Outputs at index icol: total energy, enthalpy (pressure coordinate) or
    ! internal energy (z coordinate), PHIS or potential energy term, kinetic
    ! energy [J m-2]; water vapor, total water, liquid, ice [kg m-2].
    real(kind_phys), intent(inout), optional :: te(:), se(:), po(:), ke(:)
    real(kind_phys), intent(inout), optional :: wv(:), H2O(:), liq(:), ice(:)

    real(kind_phys) :: ke_vint, se_vint, po_vint   ! vertical integrals
    real(kind_phys) :: wv_vint, liq_vint, ice_vint
    real(kind_phys) :: te_col                       ! total energy
    real(kind_phys) :: latsub                       ! latent heat of sublimation
    integer         :: kdx, qdx

    errflg = 0

    ke_vint = 0._kind_phys
    se_vint = 0._kind_phys

    select case (vcoord)
    case (ENERGY_FORMULA_DYCORE_FV, ENERGY_FORMULA_DYCORE_SE)
      po_vint = ptop(icol, 1)
      do kdx = 1, size(tracer, 2)
        ke_vint = ke_vint + (pdel(icol, kdx) *                              &
             0.5_kind_phys * (U(icol, kdx)**2 + V(icol, kdx)**2)) * rga
        se_vint = se_vint + (T(icol, kdx) *                                 &
             cp_or_cv(icol, kdx) * pdel(icol, kdx) * rga)
        po_vint = po_vint + pdel(icol, kdx)
      end do
      po_vint = (phis(icol) * po_vint * rga)
    case (ENERGY_FORMULA_DYCORE_MPAS)
      po_vint = 0._kind_phys
      do kdx = 1, size(tracer, 2)
        ke_vint = ke_vint + (pdel(icol, kdx) *                              &
             0.5_kind_phys * (U(icol, kdx)**2 + V(icol, kdx)**2) * rga)
        se_vint = se_vint + (T(icol, kdx) *                                 &
             cp_or_cv(icol, kdx) * pdel(icol, kdx) * rga)
        ! z_mid is height above ground
        po_vint = po_vint + (z_mid(icol, kdx) +                             &
             phis(icol) * rga) * pdel(icol, kdx)
      end do
    case default
      errflg = HYDROSTATIC_ENERGY_BAD_FORMULA
      return
    end select

    te_col = se_vint + po_vint + ke_vint
    if (present(se)) se(icol) = se_vint
    if (present(po)) po(icol) = po_vint
    if (present(ke)) ke(icol) = ke_vint

    ! vertical integral of water vapor
    wv_vint = 0._kind_phys
    do kdx = 1, size(tracer, 2)
      wv_vint = wv_vint + (tracer(icol, kdx, wv_idx) * pdel(icol, kdx) * rga)
    end do
    if (present(wv)) wv(icol) = wv_vint

    ! vertical integral of total liquid water (species in order)
    liq_vint = 0._kind_phys
    do qdx = 1, size(liq_idx)
      do kdx = 1, size(tracer, 2)
        liq_vint = liq_vint + (pdel(icol, kdx) *                            &
             tracer(icol, kdx, liq_idx(qdx)) * rga)
      end do
    end do
    if (present(liq)) liq(icol) = liq_vint

    ! vertical integral of total frozen (ice) water (species in order)
    ice_vint = 0._kind_phys
    do qdx = 1, size(ice_idx)
      do kdx = 1, size(tracer, 2)
        ice_vint = ice_vint + (pdel(icol, kdx) *                            &
             tracer(icol, kdx, ice_idx(qdx)) * rga)
      end do
    end do
    if (present(ice)) ice(icol) = ice_vint

    ! vertical integral of total water
    if (present(H2O)) H2O(icol) = wv_vint + liq_vint + ice_vint

    ! latent heat terms depend on enthalpy reference state
    if (present(te)) then
      latsub = latvap + latice
      select case (enthalpy_reference)
      case (ENTHALPY_REFERENCE_ICE)
        te_col = te_col + (latsub * wv_vint) + (latice * liq_vint)
      case (ENTHALPY_REFERENCE_LIQ)
        te_col = te_col + (latvap * wv_vint) - (latice * ice_vint)
      case (ENTHALPY_REFERENCE_WV)
        te_col = te_col - (latvap * liq_vint) - (latsub * ice_vint)
      case default
        errflg = HYDROSTATIC_ENERGY_BAD_REFERENCE
        return
      end select
      te(icol) = te_col
    end if

  end subroutine hydrostatic_energy_column

end module check_energy_hydrostatic_energy
