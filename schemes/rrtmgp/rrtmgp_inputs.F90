module rrtmgp_inputs
 implicit none
 private

 public :: rrtmgp_inputs_run

 contains

!> \section arg_table_rrtmgp_inputs_run Argument Table
!! \htmlinclude rrtmgp_inputs_run.html
!!
subroutine rrtmgp_inputs_run(dosw, dolw, do_snow, do_graupel, trick_rrtmgp,    &
                  pmid, pint, t, nday, idxday, coszrs, kdist_sw, kdist_lw,     &
                  lwup, stebol, ncol, ktopcam, ktoprad, nswbands, asdir, asdif,&
                  sw_low_bounds, sw_high_bounds, aldir, aldif, nlay, pverp,    &
                  pver, cld, cldfsnow, cldfgrau, graupel_in_rad, gaslist_lc,   &
                  cldfprime, t_sfc, emis_sfc, t_rad, pmid_rad, pint_rad, t_day,&
                  pmid_day, pint_day, coszrs_day, alb_dir, alb_dif, gas_concs_lw, &
                  aer_lw, atm_optics_lw, sources_lw, aer_sw, atm_optics_sw,       &
                  gas_concs_sw, errmsg, errflg)
     use ccpp_kinds,              only: kind_phys
     use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
     use ccpp_optical_props,      only: ty_optical_props_1scl_ccpp, ty_optical_props_2str_ccpp
     use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
     use ccpp_source_functions,   only: ty_source_func_lw_ccpp
     use atmos_phys_rad_utils,    only: is_visible
     use mo_rte_util_array,       only: zero_array

     ! Inputs
     logical,                              intent(in) :: graupel_in_rad        ! Flag to include graupel in radiation calculation
     integer,                              intent(in) :: ncol                  ! Number of columns
     integer,                              intent(in) :: pver                  ! Number of vertical layers
     integer,                              intent(in) :: pverp                 ! Number of vertical interfaces
     integer,                              intent(in) :: nlay                  ! Number of vertical layers used in radiation calculation
     integer,                              intent(in) :: nswbands              ! Number of shortwave bands
     integer,                              intent(in) :: ktopcam               ! Index in host model arrays of top level (layer or interface) at which RRTMGP is active
     integer,                              intent(in) :: ktoprad               ! Index in RRTMGP array corresponding to top layer or interface of host model arrays
     integer,                              intent(in) :: nday                  ! Number of daylight columns
     logical,                              intent(in) :: dosw                  ! Flag for performing the shortwave calculation
     logical,                              intent(in) :: dolw                  ! Flag for performing the longwave calculation
     logical,                              intent(in) :: do_snow               ! Flag for whether the cloud snow fraction argument should be used
     logical,                              intent(in) :: do_graupel            ! Flag for whether the cloud graupel fraction argument should be used
     logical,                              intent(in) :: trick_rrtmgp          ! Flag for whether to trick RRTMGP levels
     integer,         dimension(:),        intent(in) :: idxday                ! Indices of daylight columns
     real(kind_phys), dimension(:,:),      intent(in) :: pmid                  ! Air pressure at midpoint (Pa)
     real(kind_phys), dimension(:,:),      intent(in) :: pint                  ! Air pressure at interface (Pa)
     real(kind_phys), dimension(:,:),      intent(in) :: t                     ! Air temperature (K)
     real(kind_phys), dimension(:,:),      intent(in) :: cldfsnow              ! Cloud fraction of just "snow clouds"
     real(kind_phys), dimension(:,:),      intent(in) :: cldfgrau              ! Cloud fraction of just "graupel clouds"
     real(kind_phys), dimension(:,:),      intent(in) :: cld                   ! Cloud fraction (liq+ice)
     real(kind_phys), dimension(:),        intent(in) :: sw_low_bounds         ! Lower bounds for shortwave bands
     real(kind_phys), dimension(:),        intent(in) :: sw_high_bounds        ! Upper bounds for shortwave bands
     real(kind_phys), dimension(:),        intent(in) :: coszrs                ! Cosine of solar senith angle (radians)
     real(kind_phys), dimension(:),        intent(in) :: lwup                  ! Longwave up flux (W m-2)
     real(kind_phys), dimension(:),        intent(in) :: asdir                 ! Shortwave direct albedo (fraction)
     real(kind_phys), dimension(:),        intent(in) :: asdif                 ! Shortwave diffuse albedo (fraction)
     real(kind_phys), dimension(:),        intent(in) :: aldir                 ! Longwave direct albedo (fraction)
     real(kind_phys), dimension(:),        intent(in) :: aldif                 ! Longwave diffuse albedo (fraction)
     real(kind_phys),                      intent(in) :: stebol                ! Stefan-Boltzmann constant (W m-2 K-4)
     type(ty_gas_optics_rrtmgp_ccpp),      intent(in) :: kdist_sw              ! Shortwave gas optics object
     type(ty_gas_optics_rrtmgp_ccpp),      intent(in) :: kdist_lw              ! Longwave gas optics object
     character(len=*), dimension(:),       intent(in) :: gaslist_lc            ! Radiatively active gases
     ! Outputs
     real(kind_phys), dimension(:,:),      intent(out) :: t_rad                ! Air temperature with radiation indexing (K)
     real(kind_phys), dimension(:,:),      intent(out) :: pmid_rad             ! Midpoint pressure with radiation indexing (Pa)
     real(kind_phys), dimension(:,:),      intent(out) :: pint_rad             ! Interface pressure with radiation indexing (Pa)
     real(kind_phys), dimension(:,:),      intent(out) :: t_day                ! Air temperature of daylight columns (K)
     real(kind_phys), dimension(:,:),      intent(out) :: pint_day             ! Interface pressure of daylight columns (Pa)
     real(kind_phys), dimension(:,:),      intent(out) :: pmid_day             ! Midpoint pressure of daylight columns (Pa)
     real(kind_phys), dimension(:,:),      intent(out) :: emis_sfc             ! Surface emissivity (fraction)
     real(kind_phys), dimension(:,:),      intent(out) :: alb_dir              ! Surface albedo due to UV and VIS direct (fraction)
     real(kind_phys), dimension(:,:),      intent(out) :: alb_dif              ! Surface albedo due to IR diffused (fraction)
     real(kind_phys), dimension(:,:),      intent(out) :: cldfprime            ! modified cloud fraciton

     real(kind_phys), dimension(:),        intent(out) :: t_sfc                ! Surface temperature (K)
     real(kind_phys), dimension(:),        intent(out) :: coszrs_day           ! Cosine of solar zenith angle for daylight columns (radians)
     type(ty_gas_concs_ccpp),              intent(out) :: gas_concs_lw         ! Gas concentrations object for longwave radiation
     type(ty_optical_props_1scl_ccpp),     intent(out) :: atm_optics_lw        ! Atmosphere optical properties object for longwave radiation
     type(ty_optical_props_1scl_ccpp),     intent(out) :: aer_lw               ! Aerosol optical properties object for longwave radiation
     type(ty_source_func_lw_ccpp),         intent(out) :: sources_lw           ! Longwave sources object for longwave radiation
     type(ty_gas_concs_ccpp),              intent(out) :: gas_concs_sw         ! Gas concentrations object for shortwave radiation
     type(ty_optical_props_2str_ccpp),     intent(out) :: atm_optics_sw        ! Atmosphere optical properties object for shortwave radiation
     type(ty_optical_props_2str_ccpp),     intent(out) :: aer_sw               ! Aerosol optical properties object for shortwave radiation
     character(len=*),                     intent(out) :: errmsg
     integer,                              intent(out) :: errflg

     ! Local variables
     real(kind_phys) :: tref_min
     real(kind_phys) :: tref_max
     real(kind_phys) :: pint_min1, pint_min2
     integer :: idx, kdx, iband, ngpt
     integer :: band_kind(nswbands)  ! 1: visible, 2: near-infrared, 3: straddles both
     logical :: ltrick_rrtmgp, all_layers

     ! Set error variables
     errmsg = ''
     errflg = 0

     if (.not. dosw .and. .not. dolw) then
        return
     end if

     !------------------------------------------------------------------------------
     ! Compile logic to determine whether it is necessary AND possible to 'trick'
     ! rrtmgp to violate its own validity limits by hacking its vertical grid.
     ! Conditions:
     !   1) top CAM interface (k=1) is above 1 Pa
     !   2) next interface (k=2) is below 1 Pa
     !   3) RRTMGP is asked to active over ALL CAM layers
     !      (nlay=pverp, e.g., set by spec p_top_for_equil_rad=0.)
     !   4) dycore is NOT MPAS
     !
     ! These conditions are generally only satisfied in a non-MPAS MT configuration
     !------------------------------------------------------------------------------
     ltrick_rrtmgp = .false.
     if (trick_rrtmgp .and. (nlay==pverp)) then
        ! minval(pint(:,1)) and minval(pint(:,2)), on the device (min is exact
        ! in any order)
        pint_min1 = huge(1._kind_phys)
        pint_min2 = huge(1._kind_phys)
        !$acc parallel loop default(present) reduction(min:pint_min1, pint_min2)
        do idx = 1, size(pint, 1)
           pint_min1 = min(pint_min1, pint(idx,1))
           pint_min2 = min(pint_min2, pint(idx,2))
        end do
        ! we can and need to trick rrtmgp
        ltrick_rrtmgp = (pint_min1 < 1._kind_phys) .and. (pint_min2 > 1._kind_phys)
     end if

     ! RRTMGP set state
     ! Surface temp set based on longwave up flux.
     !$acc parallel loop default(present)
     do idx = 1, size(t_sfc)
        t_sfc(idx) = sqrt(sqrt(lwup(idx)/stebol))
     end do

     ! Set surface emissivity to 1.0.
     ! The land model *does* have its own surface emissivity, but is not spectrally resolved.
     ! The LW upward flux is calculated with that land emissivity, and the "radiative temperature"
     ! t_sfc is derived from that flux. We assume, therefore, that the emissivity is unity
     ! to be consistent with t_sfc.
     !$acc parallel loop collapse(2) default(present)
     do idx = 1, size(emis_sfc, 2)
        do iband = 1, size(emis_sfc, 1)
           emis_sfc(iband,idx) = 1._kind_phys
        end do
     end do

     !-------------------------------------------------------------------------
     ! RRTMGP enforces  P > 1 Pa for validity.
     ! Actual range of RRTMGP in CAM is set with namelist variable p_top_for_equil_rad.
     ! In rrtmg_inputs_setup.F90, active layers for RRTMGP are counted based on
     ! the logical P_ref > p_top_for_equil_rad. Returned as nlay.
     !
     ! If 
     ! 1) entire vertical domain has P_ref> p_top_for_equil_rad then
     !    nlay = pverp
     !    ktoprad = 2
     !    ktopcam = 1
     ! 2) min(P_ref) < p_top_for_equil_rad (e.g. MPAS MT) then
     !    nlay < pverp
     !    ktoprad = 1
     !    ktopcam = pver - nlay + 1
     !
     ! Level ordering is the same for both CAM and RRTMGP (top to bottom)
     ! Note in Case 2, tops of {t,pint,pmid}_rad start at a 'valid' level,
     ! i.e. pref > p_top_for_rrtmgp
     !-------------------------------------------------------
     !$acc parallel loop collapse(2) default(present)
     do kdx = ktoprad, size(t_rad, 2)
        do idx = 1, size(t_rad, 1)
           t_rad(idx,kdx)    = t(idx,kdx-ktoprad+ktopcam)
           pmid_rad(idx,kdx) = pmid(idx,kdx-ktoprad+ktopcam)
        end do
     end do
     !$acc parallel loop collapse(2) default(present)
     do kdx = ktoprad, size(pint_rad, 2)
        do idx = 1, size(pint_rad, 1)
           pint_rad(idx,kdx) = pint(idx,kdx-ktoprad+ktopcam)
        end do
     end do

     ! Deal with vertical grid for RRTMGP 
     ! Case 1 (nlay == pverp): all model layers are within RRTMGP's range of
     ! valid pressures as specified by p_top_for_equil_rad.
     ! Case 2 (nlay < pverp): model min(pref) < p_top_for_rrtmgp.
     all_layers = (nlay == pverp)
     !$acc parallel loop default(present)
     do idx = 1, size(t_rad, 1)
        if (all_layers) then
           t_rad(idx,1)      = t(idx,1)
           if (ltrick_rrtmgp) then
              ! The top reference pressure from the RRTMGP coefficients datasets is 1.005183574463 Pa
              ! Set the top of the extra layer just below that.
              pint_rad(idx,1) = 1.01_kind_phys
              pint_rad(idx,2) = 1.02_kind_phys
              ! set the highest pmid (in the "extra layer") to the midpoint (guarantees > 1Pa)
              pmid_rad(idx,1)   = 1.015_kind_phys ! pint_rad(:,1) + 0.5_kind_phys * (pint_rad(:,2) - pint_rad(:,1))
              pmid_rad(idx,2)   = 0.5_kind_phys*(pint_rad(idx,2) + pint_rad(idx,3))
           else
              ! The top reference pressure from the RRTMGP coefficients datasets is 1.005183574463 Pa
              ! Set the top of the extra layer just below that.
              pint_rad(idx,1) = 1.01_kind_phys
              ! set the highest pmid (in the "extra layer") to the midpoint (guarantees > 1Pa)
              pmid_rad(idx,1)   = pint_rad(idx,1) + 0.5_kind_phys * (pint_rad(idx,2) - pint_rad(idx,1))
           end if
        else
           ! Not sure why is this needed since pint_rad should have been specified
           ! at RRTMGP valid values above. But this is the way it was done in
           ! original RRTMG implementation.
           pint_rad(idx,1) = 1.01_kind_phys
           ! The following *should* work since pint_rad is all in valid range.
           pmid_rad(idx,1) = 0.5_kind_phys * (pint_rad(idx,1) + pint_rad(idx,2))
        end if
     end do

     ! Limit temperatures to be within the limits of RRTMGP validity.
     tref_min = kdist_sw%gas_props%get_temp_min()
     tref_max = kdist_sw%gas_props%get_temp_max()
     !$acc parallel loop collapse(2) default(present)
     do kdx = 1, size(t_rad, 2)
        do idx = 1, size(t_rad, 1)
           t_rad(idx,kdx) = merge(t_rad(idx,kdx), tref_min, t_rad(idx,kdx) > tref_min)
           t_rad(idx,kdx) = merge(t_rad(idx,kdx), tref_max, t_rad(idx,kdx) < tref_max)
        end do
     end do

     ! Construct arrays containing only daylight columns
     !$acc parallel loop collapse(2) default(present)
     do kdx = 1, size(t_day, 2)
        do idx = 1, nday
           t_day(idx,kdx)    = t_rad(idxday(idx),kdx)
           pmid_day(idx,kdx) = pmid_rad(idxday(idx),kdx)
        end do
     end do
     !$acc parallel loop collapse(2) default(present)
     do kdx = 1, size(pint_day, 2)
        do idx = 1, nday
           pint_day(idx,kdx) = pint_rad(idxday(idx),kdx)
        end do
     end do
     !$acc parallel loop default(present)
     do idx = 1, nday
        coszrs_day(idx) = coszrs(idxday(idx))
     end do

     ! Assign albedos to the daylight columns (from E3SM implementation)
     ! Albedos are imported from the surface models as broadband (visible, and near-IR),
     ! and we need to map these to appropriate narrower bands used in RRTMGP. Bands
     ! are categorized broadly as "visible/UV" or "infrared" based on wavenumber.
     ! Loop over bands, and determine for each band whether it is broadly in the
     ! visible or infrared part of the spectrum based on a dividing line of
     ! 0.7 micron, or 14286 cm^-1 (on the host: the band bounds only)
     do iband = 1,nswbands
        if (is_visible(sw_low_bounds(iband)) .and. &
           is_visible(sw_high_bounds(iband))) then
           ! Entire band is in the visible
           band_kind(iband) = 1
        else if (.not.is_visible(sw_low_bounds(iband)) .and. &
                 .not.is_visible(sw_high_bounds(iband))) then
           ! Entire band is in the longwave (near-infrared)
           band_kind(iband) = 2
        else
           ! Band straddles the visible to near-infrared transition, so we take
           ! the albedo to be the average of the visible and near-infrared
           ! broadband albedos
           band_kind(iband) = 3
        end if
     end do
     !$acc parallel loop collapse(2) default(present) copyin(band_kind)
     do idx = 1, nday
        do iband = 1, nswbands
           if (band_kind(iband) == 1) then
              alb_dir(iband,idx) = asdir(idxday(idx))
              alb_dif(iband,idx) = asdif(idxday(idx))
           else if (band_kind(iband) == 2) then
              alb_dir(iband,idx) = aldir(idxday(idx))
              alb_dif(iband,idx) = aldif(idxday(idx))
           else
              alb_dir(iband,idx) = 0.5_kind_phys * (aldir(idxday(idx)) + asdir(idxday(idx)))
              alb_dif(iband,idx) = 0.5_kind_phys * (aldif(idxday(idx)) + asdif(idxday(idx)))
           end if
        end do
     end do
     ! Strictly enforce albedo bounds
     !$acc parallel loop collapse(2) default(present)
     do idx = 1, size(alb_dir, 2)
        do iband = 1, size(alb_dir, 1)
           if (alb_dir(iband,idx) < 0) alb_dir(iband,idx) = 0.0_kind_phys
           if (alb_dir(iband,idx) > 1) alb_dir(iband,idx) = 1.0_kind_phys
        end do
     end do
     !$acc parallel loop collapse(2) default(present)
     do idx = 1, size(alb_dif, 2)
        do iband = 1, size(alb_dif, 1)
           if (alb_dif(iband,idx) < 0) alb_dif(iband,idx) = 0.0_kind_phys
           if (alb_dif(iband,idx) > 1) alb_dif(iband,idx) = 1.0_kind_phys
        end do
     end do

     ! modified cloud fraction
     ! Compute modified cloud fraction, cldfprime.
     ! 1. initialize as cld
     ! 2. modify for snow. use max(cld, cldfsnow)
     ! 3. modify for graupel if graupel_in_rad is true.
     !    use max(cldfprime, cldfgrau)
     if (do_snow) then
        !$acc parallel loop collapse(2) default(present)
        do kdx = 1, pver
           do idx = 1, ncol
              cldfprime(idx,kdx) = max(cld(idx,kdx), cldfsnow(idx,kdx))
           end do
        end do
     else
        !$acc parallel loop collapse(2) default(present)
        do kdx = 1, size(cldfprime, 2)
           do idx = 1, size(cldfprime, 1)
              cldfprime(idx,kdx) = cld(idx,kdx)
           end do
        end do
     end if

     if (do_graupel .and. graupel_in_rad) then
        !$acc parallel loop collapse(2) default(present)
        do kdx = 1, pver
           do idx = 1, ncol
              cldfprime(idx,kdx) = max(cldfprime(idx,kdx), cldfgrau(idx,kdx))
           end do
        end do
     end if

     ! If no daylight columns, can't create empty RRTMGP objects
     if (dosw .and. nday > 0) then
        ! Initialize object for gas concentrations.
        errmsg = gas_concs_sw%init(gaslist_lc)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

        ! Initialize object for combined gas + aerosol + cloud optics.
        ! Allocates arrays for properties represented on g-points.
        errmsg = atm_optics_sw%optical_props%alloc_2str(nday, nlay, kdist_sw%gas_props)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if
        ! Keep the g-point arrays on the device for the whole radiation call:
        ! created here, deleted in rrtmgp_post (rrtmgp-rte-object-residency)
        !$acc enter data copyin(atm_optics_sw%optical_props)
        !$acc enter data create(atm_optics_sw%optical_props%tau, atm_optics_sw%optical_props%ssa, &
        !$acc                   atm_optics_sw%optical_props%g)
        ! Allocation leaves the values undefined; define them (on the device)
        ngpt = atm_optics_sw%optical_props%get_ngpt()
        call zero_array(nday, nlay, ngpt, atm_optics_sw%optical_props%tau)
        call zero_array(nday, nlay, ngpt, atm_optics_sw%optical_props%ssa)
        call zero_array(nday, nlay, ngpt, atm_optics_sw%optical_props%g)

        ! Initialize object for SW aerosol optics.  Allocates arrays
        ! for properties represented by band.
        errmsg = aer_sw%optical_props%alloc_2str(nday, nlay, kdist_sw%gas_props%get_band_lims_wavenumber())
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if
        ! Keep the g-point arrays on the device for the whole radiation call:
        ! created here, deleted in rrtmgp_post (rrtmgp-rte-object-residency)
        !$acc enter data copyin(aer_sw%optical_props)
        !$acc enter data create(aer_sw%optical_props%tau, aer_sw%optical_props%ssa, &
        !$acc                   aer_sw%optical_props%g)
        ! Allocation leaves the values undefined; define them (on the device)
        ngpt = aer_sw%optical_props%get_ngpt()
        call zero_array(nday, nlay, ngpt, aer_sw%optical_props%tau)
        call zero_array(nday, nlay, ngpt, aer_sw%optical_props%ssa)
        call zero_array(nday, nlay, ngpt, aer_sw%optical_props%g)
     end if

     if (dolw) then
        ! Initialize object for gas concentrations
        errmsg = gas_concs_lw%init(gaslist_lc)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if

        ! Initialize object for combined gas + aerosol + cloud optics.
        errmsg = atm_optics_lw%optical_props%alloc_1scl(ncol, nlay, kdist_lw%gas_props)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if
        ! Keep the g-point arrays on the device for the whole radiation call:
        ! created here, deleted in rrtmgp_post (rrtmgp-rte-object-residency)
        !$acc enter data copyin(atm_optics_lw%optical_props)
        !$acc enter data create(atm_optics_lw%optical_props%tau)
        ! Allocation leaves the values undefined; define them (on the device)
        ngpt = atm_optics_lw%optical_props%get_ngpt()
        call zero_array(ncol, nlay, ngpt, atm_optics_lw%optical_props%tau)

        ! Initialize object for LW aerosol optics.
        errmsg = aer_lw%optical_props%alloc_1scl(ncol, nlay, kdist_lw%gas_props%get_band_lims_wavenumber())
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if
        ! Keep the g-point arrays on the device for the whole radiation call:
        ! created here, deleted in rrtmgp_post (rrtmgp-rte-object-residency)
        !$acc enter data copyin(aer_lw%optical_props)
        !$acc enter data create(aer_lw%optical_props%tau)
        ! Allocation leaves the values undefined; define them (on the device)
        ngpt = aer_lw%optical_props%get_ngpt()
        call zero_array(ncol, nlay, ngpt, aer_lw%optical_props%tau)

        ! Initialize object for Planck sources.
        errmsg = sources_lw%sources%alloc(ncol, nlay, kdist_lw%gas_props)
        if (len_trim(errmsg) > 0) then
           errflg = 1
           return
        end if
        ! Keep the source arrays on the device for the whole radiation call:
        ! created here, deleted in rrtmgp_post (rrtmgp-rte-object-residency)
        !$acc enter data copyin(sources_lw%sources)
        !$acc enter data create(sources_lw%sources%lay_source, sources_lw%sources%lev_source, &
        !$acc                   sources_lw%sources%sfc_source, sources_lw%sources%sfc_source_Jac)
        ! Allocation leaves the values undefined; define them (on the device)
        ngpt = sources_lw%sources%get_ngpt()
        call zero_array(ncol, nlay,   ngpt, sources_lw%sources%lay_source)
        call zero_array(ncol, nlay+1, ngpt, sources_lw%sources%lev_source)
        call zero_array(ncol, ngpt, sources_lw%sources%sfc_source)
        call zero_array(ncol, ngpt, sources_lw%sources%sfc_source_Jac)
     end if

  end subroutine rrtmgp_inputs_run

end module rrtmgp_inputs
