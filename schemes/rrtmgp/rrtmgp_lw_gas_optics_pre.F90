module rrtmgp_lw_gas_optics_pre
  implicit none
  private

  public :: rrtmgp_lw_gas_optics_pre_run

contains

!> \section arg_table_rrtmgp_lw_gas_optics_pre_run Argument Table
!! \htmlinclude rrtmgp_lw_gas_optics_pre_run.html
!!
  subroutine rrtmgp_lw_gas_optics_pre_run(rad_const_array, pmid, pint, nlay, ncol, gaslist, &
                  pverp, ktoprad, ktopcam, dolw, nradgas, gas_concs, errmsg, errflg)
    use ccpp_kinds,              only: kind_phys
    use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
    use radiation_utils,         only: get_molar_mass_ratio

    ! Set gas vmr for the gases in the radconstants module's gaslist.

    character(len=5),            intent(in) :: gaslist(:)             ! Radiatively active gases
    integer,                     intent(in) :: nlay                   ! Number of layers in radiation calculation
    integer,                     intent(in) :: ncol                   ! Total number of columns
    integer,                     intent(in) :: pverp                  ! Total number of layer interfaces
    integer,                     intent(in) :: ktoprad                ! Index in RRTMGP array corresponding to top layer or interface of host model arrays
    integer,                     intent(in) :: ktopcam                ! Index in host model arrays of top level (layer or interface) at which RRTMGP is active
    integer,                     intent(in) :: nradgas                ! Number of radiatively active gases
    logical,                     intent(in) :: dolw                   ! Flag for whether to perform longwave calculaion
    real(kind_phys),             intent(in) :: pmid(:,:)              ! Air pressure at midpoints [Pa]
    real(kind_phys),             intent(in) :: pint(:,:)              ! Air pressure at interfaces [Pa]
    real(kind_phys),             intent(in) :: rad_const_array(:,:,:) ! array of radiatively-active constituent vmrs
                                                                      !  last index corresponds to index in gaslist

    type(ty_gas_concs_ccpp),     intent(inout) :: gas_concs           ! the result is VMR inside gas_concs
    character(len=*),            intent(out)   :: errmsg
    integer,                     intent(out)   :: errflg

    ! Local variables
    integer :: i, k, gas_idx
    integer :: istat
    logical :: is_h2o
    real(kind_phys)              :: gas_vmr(ncol, nlay)
    real(kind_phys)              :: mmr(ncol, nlay)
    real(kind_phys) :: massratio
    character(len=256) :: alloc_errmsg

    ! For ozone profile above model
    real(kind_phys) :: P_top, P_int, P_mid, alpha, beta, a, b, chi_mid, chi_0, chi_eff

    character(len=*), parameter :: sub = 'rrtmgp_lw_gas_optics_pre_run'
    !----------------------------------------------------------------------------

    ! Set error variables
    errmsg = ''
    errflg = 0

    if (.not. dolw) then
       return
    end if

    ! gas_vmr is computed on the device, where set_vmr copies it into
    ! gas_concs (RRTMGP keeps the concentrations on the device).
    !$acc data create(mmr, gas_vmr)

    do gas_idx = 1, nradgas

       ! grab mass mixing ratio of gas (all columns for longwave)
       !$acc parallel loop collapse(2) default(present)
       do k = ktoprad, nlay
          do i = 1, ncol
             mmr(i,k) = rad_const_array(i,k-ktoprad+ktopcam,gas_idx)
          end do
       end do

       ! If an extra layer is being used, copy mmr from the top layer of CAM to the extra layer.
       if (nlay == pverp) then
          !$acc parallel loop default(present)
          do i = 1, ncol
             mmr(i,1) = mmr(i,2)
          end do
       end if

       ! convert MMR to VMR, multipy by ratio of dry air molar mas to gas molar mass.
       call get_molar_mass_ratio(gaslist(gas_idx), massratio, errmsg, errflg)
       if (errflg /= 0) then
          exit
       end if

       ! special case: H2O is specific humidity, not mixing ratio. Use r = q/(1-q):
       is_h2o = gaslist(gas_idx) == 'H2O'
       !$acc parallel loop collapse(2) default(present)
       do k = 1, nlay
          do i = 1, ncol
             if (is_h2o) then
                mmr(i,k) = mmr(i,k) / (1._kind_phys - mmr(i,k))
             end if
             gas_vmr(i,k) = mmr(i,k) * massratio
          end do
       end do

       ! special case: Setting O3 in the extra layer:
       ! 
       ! For the purpose of attenuating solar fluxes above the CAM model top, we assume that ozone 
       ! mixing decreases linearly in each column from the value in the top layer of CAM to zero at 
       ! the pressure level set by P_top. P_top has been set to 50 Pa (0.5 hPa) based on model tuning 
       ! to produce temperatures at the top of CAM that are most consistent with WACCM at similar pressure levels. 

       if ((gaslist(gas_idx) == 'O3') .and. (nlay == pverp)) then
          P_top = 50.0_kind_phys
          !$acc parallel loop default(present) &
          !$acc private(P_int, P_mid, alpha, beta, a, b, chi_mid, chi_0, chi_eff)
          do i = 1, ncol
             P_int = pint(i,1) ! pressure (Pa) at upper interface of CAM
             P_mid = pmid(i,1) ! pressure (Pa) at midpoint of top layer of CAM
             alpha = log(P_int/P_top)
             beta =  log(P_mid/P_int)/log(P_mid/P_top)

             a =  ( (1._kind_phys + alpha) * exp(-alpha) - 1._kind_phys ) / alpha
             b =  1._kind_phys - exp(-alpha)
   
             if (alpha .gt. 0) then             ! only apply where top level is below 80 km
                chi_mid = gas_vmr(i,1)          ! molar mixing ratio of O3 at midpoint of top layer
                chi_0 = chi_mid /  (1._kind_phys + beta)
                chi_eff = chi_0 * (a + b)
                gas_vmr(i,1) = chi_eff
             end if
          end do
       end if

       errmsg = gas_concs%set_vmr(gaslist(gas_idx), gas_vmr)
       if (len_trim(errmsg) > 0) then
          errflg = 1
          exit
       end if

    end do

    !$acc end data

  end subroutine rrtmgp_lw_gas_optics_pre_run

end module rrtmgp_lw_gas_optics_pre
