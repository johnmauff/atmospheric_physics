module rrtmgp_post

  implicit none
  private

  public :: rrtmgp_post_run

contains
!> \section arg_table_rrtmgp_post_run Argument Table
!! \htmlinclude rrtmgp_post_run.html
!!
subroutine rrtmgp_post_run(nlay, dolw, qrs_prime, qrl_prime, fsns, pdel, atm_optics_sw, cloud_sw, aer_sw,  &
                  fsw, fswc, atm_optics_lw, sources_lw, cloud_lw, aer_lw, flw, flwc, flwds, qrs, qrl, &
                  netsw, gas_concs_sw, gas_concs_lw, errmsg, errflg)
   use ccpp_kinds,             only: kind_phys
   use ccpp_optical_props,     only: ty_optical_props_1scl_ccpp, ty_optical_props_2str_ccpp
   use ccpp_source_functions,  only: ty_source_func_lw_ccpp
   use ccpp_fluxes,            only: ty_fluxes_broadband_ccpp
   use ccpp_fluxes_byband,     only: ty_fluxes_byband_ccpp
   use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
   use mo_rte_kind,            only: wl
   use mo_rte_config,          only: rte_config_checks
   integer,                          intent(in)    :: nlay           ! Number of layers in radiation calculation
   logical,                          intent(in)    :: dolw           ! Flag for whether to perform longwave calculation
   real(kind_phys), dimension(:,:),  intent(in)    :: pdel           ! Layer thickness [Pa]
   real(kind_phys), dimension(:),    intent(in)    :: fsns           ! Surface net shortwave flux [W m-2]
   real(kind_phys), dimension(:,:),  intent(in)    :: qrs_prime      ! Shortwave heating rate [J kg-1 s-1]
   real(kind_phys), dimension(:,:),  intent(in)    :: qrl_prime      ! Longwave heating rate [J kg-1 s-1]
   type(ty_optical_props_1scl_ccpp), intent(inout) :: atm_optics_lw  ! Atmosphere optical properties object (longwave)
   type(ty_optical_props_2str_ccpp), intent(inout) :: atm_optics_sw  ! Atmosphere optical properties object (shortwave)
   type(ty_optical_props_1scl_ccpp), intent(inout) :: aer_lw         ! Aerosol optical properties object (longwave)
   type(ty_optical_props_2str_ccpp), intent(inout) :: aer_sw         ! Aerosol optical properties object (shortwave)
   type(ty_optical_props_1scl_ccpp), intent(inout) :: cloud_lw       ! Cloud optical properties object (longwave)
   type(ty_optical_props_2str_ccpp), intent(inout) :: cloud_sw       ! Cloud optical properties object (shortwave)
   type(ty_fluxes_broadband_ccpp),   intent(inout) :: fswc           ! Shortwave clear-sky flux object
   type(ty_fluxes_broadband_ccpp),   intent(inout) :: flwc           ! Longwave clear-sky flux object
   type(ty_fluxes_byband_ccpp),      intent(inout) :: fsw            ! Shortwave all-sky flux object
   type(ty_fluxes_byband_ccpp),      intent(inout) :: flw            ! Longwave all-sky flux object
   type(ty_source_func_lw_ccpp),     intent(inout) :: sources_lw     ! Longwave sources object
   real(kind_phys), dimension(:,:),  intent(out)   :: qrs            ! Shortwave heating rate adjusted by air pressure thickness [J Pa kg-1 s-1]
   real(kind_phys), dimension(:,:),  intent(out)   :: qrl            ! Longwave heating rate adjusted by air pressure thickness [J Pa kg-1 s-1]
   real(kind_phys), dimension(:),    intent(out)   :: netsw          ! Net shortwave flux to be sent to coupler [W m-2]
   real(kind_phys), dimension(:),    intent(inout) :: flwds          ! Down longwave flux at surface [W m-2]
   type(ty_gas_concs_ccpp),          intent(inout) :: gas_concs_sw   ! Gas concentrations object (shortwave)
   type(ty_gas_concs_ccpp),          intent(inout) :: gas_concs_lw   ! Gas concentrations object (longwave)
   character(len=*),                 intent(out)   :: errmsg
   integer,                          intent(out)   :: errflg

   integer :: i, k

   ! Set error varaibles
   errflg = 0
   errmsg = ''
   ! The radiative heating rates are maintained across multiple physics timesteps
   ! as Q*dp (for energy conservation).
   !$acc data copyin(qrs_prime, qrl_prime, pdel, fsns) copyout(qrs, qrl, netsw)
   !$acc parallel loop collapse(2)
   do k = 1, size(qrs, 2)
      do i = 1, size(qrs, 1)
         qrs(i, k) = qrs_prime(i, k) * pdel(i, k)
         qrl(i, k) = qrl_prime(i, k) * pdel(i, k)
      end do
   end do

   ! Set the netsw to be sent to the coupler
   !$acc parallel loop
   do i = 1, size(netsw)
      netsw(i) = fsns(i)
   end do
   !$acc end data

   ! On the host, with the RTE objects this scheme frees below: their device
   ! copies are a later step (rrtmgp-rte-object-residency).
   if (dolw) then
      flwds(:) = flw%fluxes%flux_dn(:, nlay+1)
   end if

   ! RTE-RRTMGP checks its input values (check_values: range checks of
   ! every optical property, source and boundary condition, each a
   ! reduction returned to the host) on the first radiation call only, as
   ! RTE's own all-sky example does; the array extent checks (host-only)
   ! stay on. The value checks change no answers, only report errors.
   if (dolw) then
      call rte_config_checks(logical(.true., wl), logical(.false., wl))
   end if

   call free_optics_sw(atm_optics_sw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_optics_sw(cloud_sw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_optics_sw(aer_sw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_fluxes_byband(fsw)
   call free_fluxes_broadband(fswc)
   ! Free the gas concentrations here, like the other objects: rrtmgp_inputs
   ! takes them intent(out), and the compiler's deallocation on entry need not
   ! run RRTMGP's finalizer, leaving their OpenACC device copies behind.
   call gas_concs_sw%reset()

   call free_sources_lw(sources_lw)
   call free_optics_lw(atm_optics_lw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_optics_lw(cloud_lw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_optics_lw(aer_lw, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   call free_fluxes_byband(flw)
   call free_fluxes_broadband(flwc)
   call gas_concs_lw%reset()

end subroutine rrtmgp_post_run

  !=========================================================================================

subroutine free_optics_sw(optics, errmsg, errflg)
   use ccpp_optical_props,      only: ty_optical_props_2str_ccpp

   type(ty_optical_props_2str_ccpp), intent(inout) :: optics
   character(len=*),                 intent(out)   :: errmsg
   integer,                          intent(out)   :: errflg

   errmsg = ''
   errflg = 0

   ! Delete the device copies first (created where the object is allocated;
   ! rrtmgp-rte-object-residency), so no device mapping outlives the host
   ! memory deallocated below
   if (allocated(optics%optical_props%tau)) then
      !$acc exit data delete(optics%optical_props%tau, optics%optical_props%ssa, &
      !$acc                  optics%optical_props%g)
      !$acc exit data delete(optics%optical_props)
   end if

   errmsg = optics%optical_props%finalize_2str()
   if (len_trim(errmsg) /= 0) then
      errflg = 1
      return
   end if

   call optics%optical_props%finalize()

end subroutine free_optics_sw

!=========================================================================================

subroutine free_optics_lw(optics, errmsg, errflg)
   use ccpp_optical_props,      only: ty_optical_props_1scl_ccpp

   type(ty_optical_props_1scl_ccpp), intent(inout) :: optics
   character(len=*),                 intent(out)   :: errmsg
   integer,                          intent(out)   :: errflg

   ! Delete the device copies first (created where the object is allocated;
   ! rrtmgp-rte-object-residency), so no device mapping outlives the host
   ! memory deallocated below
   if (allocated(optics%optical_props%tau)) then
      !$acc exit data delete(optics%optical_props%tau)
      !$acc exit data delete(optics%optical_props)
   end if

   errmsg = optics%optical_props%finalize_1scl()
   if (len_trim(errmsg) /= 0) then
      errflg = 1
      return
   end if

   call optics%optical_props%finalize()

end subroutine free_optics_lw

!=========================================================================================

subroutine free_sources_lw(sources)
   use ccpp_source_functions,   only: ty_source_func_lw_ccpp

   type(ty_source_func_lw_ccpp), intent(inout) :: sources

   ! Delete the device copies first (created in rrtmgp_inputs;
   ! rrtmgp-rte-object-residency), so no device mapping outlives the host
   ! memory deallocated below
   if (allocated(sources%sources%lay_source)) then
      !$acc exit data delete(sources%sources%lay_source, sources%sources%lev_source, &
      !$acc                  sources%sources%sfc_source, sources%sources%sfc_source_Jac)
      !$acc exit data delete(sources%sources)
   end if

   call sources%sources%finalize()

end subroutine free_sources_lw

!=========================================================================================

subroutine free_fluxes_broadband(fluxes)
   use ccpp_fluxes,            only: ty_fluxes_broadband_ccpp

   class(ty_fluxes_broadband_ccpp), intent(inout) :: fluxes

   if (associated(fluxes%fluxes%flux_up)) deallocate(fluxes%fluxes%flux_up)
   if (associated(fluxes%fluxes%flux_dn)) deallocate(fluxes%fluxes%flux_dn)
   if (associated(fluxes%fluxes%flux_net)) deallocate(fluxes%fluxes%flux_net)
   if (associated(fluxes%fluxes%flux_dn_dir)) deallocate(fluxes%fluxes%flux_dn_dir)

end subroutine free_fluxes_broadband

!=========================================================================================

subroutine free_fluxes_byband(fluxes)
   use ccpp_fluxes_byband,     only: ty_fluxes_byband_ccpp

   class(ty_fluxes_byband_ccpp), intent(inout) :: fluxes

   if (associated(fluxes%fluxes%flux_up)) deallocate(fluxes%fluxes%flux_up)
   if (associated(fluxes%fluxes%flux_dn)) deallocate(fluxes%fluxes%flux_dn)
   if (associated(fluxes%fluxes%flux_net)) deallocate(fluxes%fluxes%flux_net)
   if (associated(fluxes%fluxes%flux_dn_dir)) deallocate(fluxes%fluxes%flux_dn_dir)

   if (associated(fluxes%fluxes%bnd_flux_up)) deallocate(fluxes%fluxes%bnd_flux_up)
   if (associated(fluxes%fluxes%bnd_flux_dn)) deallocate(fluxes%fluxes%bnd_flux_dn)
   if (associated(fluxes%fluxes%bnd_flux_net)) deallocate(fluxes%fluxes%bnd_flux_net)
   if (associated(fluxes%fluxes%bnd_flux_dn_dir)) deallocate(fluxes%fluxes%bnd_flux_dn_dir)

end subroutine free_fluxes_byband

end module rrtmgp_post
