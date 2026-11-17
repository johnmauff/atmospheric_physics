!> \file rrtmgp_sw_cloud_optics.F90
!!
module rrtmgp_sw_cloud_optics

!--------------------------------------------------------------------------------
! Transform data for inputs from CAM's data structures to those used by
! RRTMGP.  Subset the number of model levels if CAM's top exceeds RRTMGP's
! valid domain.  Add an extra layer if CAM's top is below 1 Pa.
! The vertical indexing increases from top to bottom of atmosphere in both
! CAM and RRTMGP arrays.   
!--------------------------------------------------------------------------------

implicit none
private
save

public :: rrtmgp_sw_cloud_optics_run

! Mapping from RRTMG shortwave bands to RRTMGP.  Currently needed to continue using
! the SW optics datasets from RRTMG (even thought there is a slight mismatch in the
! band boundaries of the 2 bands that overlap with the LW bands).
integer, parameter, dimension(14) :: rrtmg_to_rrtmgp_swbands = &
   [ 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ]

!==================================================================================================
contains
!==================================================================================================

!> \section arg_table_rrtmgp_sw_cloud_optics_run Argument Table
!! \htmlinclude rrtmgp_sw_cloud_optics_run.html
!!
subroutine rrtmgp_sw_cloud_optics_run(dosw, ncol, pver, ktopcam, ktoprad,  nswgpts, nday, idxday, fillvalue, &
   nswbands, iulog, pgam, lamc, rel, rei, nnite, idxnite, cld, cldfsnow, cldfgrau, cldfprime, &
   degrau, dei, des, iclwpth, iciwpth, icswpth, icgrauwpth, tiny_in, idx_sw_diag, do_graupel, &
   do_snow, kdist_sw, cld_tau, grau_tau, snow_tau, c_cld_tau, c_cld_tau_w, c_cld_tau_w_g, tot_cld_vistau,    &
   tot_icld_vistau, liq_icld_vistau, ice_icld_vistau, snow_icld_vistau, grau_icld_vistau, errmsg, errflg)
   use ccpp_gas_optics_rrtmgp,    only: ty_gas_optics_rrtmgp_ccpp
   use ccpp_optical_props,        only: ty_optical_props_2str_ccpp
   use rrtmgp_cloud_optics_setup, only: g_mu, g_lambda, nmu, nlambda, g_d_eff, n_g_d
   use rrtmgp_cloud_optics_setup, only: ext_sw_liq, asm_sw_liq, ssa_sw_liq
   use rrtmgp_cloud_optics_setup, only: ext_sw_ice, asm_sw_ice, ssa_sw_ice
   use rrtmgp_cloud_optics_setup, only: liq_cld_optics, ice_cld_optics
   use ccpp_kinds,                only: kind_phys

   ! Compute combined cloud optical properties.

   ! arguments
   integer,  intent(in) :: ncol                      ! Total number of columns
   integer,  intent(in) :: nday                      ! Number of daylight columns
   integer,  intent(in) :: idxday(:)                 ! Indices of daylight columns
   integer,  intent(in) :: nswgpts                   ! Number of shortwave g-points
   integer,  intent(in) :: pver                      ! Number of vertical layers
   integer,  intent(in) :: ktopcam                   ! Index in host model arrays of top level (layer or interface) at which RRTMGP is active
   integer,  intent(in) :: ktoprad                   ! Index in RRTMGP array corresponding to top layer or interface of host model arrays
   integer,  intent(in) :: nswbands                  ! Number of shortwve bands
   integer,  intent(in) :: nnite                     ! Number of night columns
   integer,  intent(in) :: idxnite(:)                ! Indices of night columns in the chunk
   integer,  intent(in) :: iulog                     ! Logging unit
   integer,  intent(in) :: idx_sw_diag               ! Index for band that contains 500-nm wave

   logical,  intent(in) :: do_snow                   ! Flag to include snow in radiation calculation
   logical,  intent(in) :: do_graupel                ! Flag to include graupel in radiation calculation
   logical,  intent(in) :: dosw                      ! Flag to do shortwave radiation this timestep

   real(kind_phys), intent(in) :: fillvalue          ! Fill value for night columns
   real(kind_phys), intent(in) :: tiny_in            ! Definition of tiny for RRTMGP

   real(kind_phys), intent(in) :: lamc(:,:)          ! Prognosed value of lambda for cloud [1]
   real(kind_phys), intent(in) :: pgam(:,:)          ! Prognosed value of mu for cloud [1]
   real(kind_phys), intent(in) :: rel(:,:)           ! Effective radius of stratiform cloud liquid water droplet [um]
   real(kind_phys), intent(in) :: rei(:,:)           ! Effective radius of stratiform cloud ice crystal [um]
   real(kind_phys), intent(in) :: dei(:,:)           ! Mean effective radius for ice cloud [um]
   real(kind_phys), intent(in) :: des(:,:)           ! Mean effective radius for snow [um]
   real(kind_phys), intent(in) :: degrau(:,:)        ! Mean effective radius for graupel [um]
   real(kind_phys), intent(in) :: iclwpth(:,:)       ! In-cloud liquid water path [kg m-2]
   real(kind_phys), intent(in) :: iciwpth(:,:)       ! In-cloud ice water path [kg m-2]
   real(kind_phys), intent(in) :: icswpth(:,:)       ! In-cloud snow water path [kg m-2]
   real(kind_phys), intent(in) :: icgrauwpth(:,:)    ! In-cloud graupel water path [kg m-2]
   real(kind_phys), intent(in) :: cld(:,:)           ! Cloud fraction (liq+ice) [fraction]
   real(kind_phys), intent(in) :: cldfsnow(:,:)      ! Cloud fraction of just "snow clouds" [fraction]
   real(kind_phys), intent(in) :: cldfgrau(:,:)      ! Cloud fraction of just "graupel clouds" [fraction]
   real(kind_phys), intent(in) :: cldfprime(:,:)     ! Combined cloud fraction [fraction]

   class(ty_gas_optics_rrtmgp_ccpp), intent(in)  :: kdist_sw             ! shortwave gas optics object
   real(kind_phys),                  intent(out) :: cld_tau(:,:,:)       ! liquid + ice optical depth
   real(kind_phys),                  intent(out) :: snow_tau(:,:,:)      ! snow optical depth
   real(kind_phys),                  intent(out) :: grau_tau(:,:,:)      ! graupel optical depth
   real(kind_phys),                  intent(out) :: c_cld_tau(:,:,:)     ! combined cloud extinction optical depth
   real(kind_phys),                  intent(out) :: c_cld_tau_w  (:,:,:) ! combined cloud single scattering albedo * tau
   real(kind_phys),                  intent(out) :: c_cld_tau_w_g(:,:,:) ! combined cloud asymmetry parameter * w * tau

   ! Diagnostic outputs
   real(kind_phys), intent(out) :: tot_cld_vistau(:,:)   ! gbx total cloud optical depth
   real(kind_phys), intent(out) :: tot_icld_vistau(:,:)  ! in-cld total cloud optical depth
   real(kind_phys), intent(out) :: liq_icld_vistau(:,:)  ! in-cld liq cloud optical depth
   real(kind_phys), intent(out) :: ice_icld_vistau(:,:)  ! in-cld ice cloud optical depth
   real(kind_phys), intent(out) :: snow_icld_vistau(:,:) ! snow in-cloud visible sw optical depth
   real(kind_phys), intent(out) :: grau_icld_vistau(:,:) ! Graupel in-cloud visible sw optical depth

   ! Error variables
   character(len=*),   intent(out) :: errmsg
   integer,            intent(out) :: errflg


   ! Local variables

   integer :: i, k, ns, pb
   integer :: ngrau_warn
   real(kind_phys) :: ctau, ctau_w, ctau_w_g

   ! cloud radiative parameters are "in cloud" not "in cell"
   ! All of these are in RRTMG band order (that of the optics datasets); the
   ! outputs are written in RRTMGP band order through rrtmg_to_rrtmgp_swbands,
   ! except cld_tau, which stays in RRTMG order.
   real(kind_phys) :: liq_tau    (nswbands,ncol,pver)  ! liquid extinction optical depth
   real(kind_phys) :: liq_tau_w  (nswbands,ncol,pver)  ! liquid single scattering albedo * tau
   real(kind_phys) :: liq_tau_w_g(nswbands,ncol,pver)  ! liquid asymmetry parameter * tau * w
   real(kind_phys) :: ice_tau    (nswbands,ncol,pver)  ! ice extinction optical depth
   real(kind_phys) :: ice_tau_w  (nswbands,ncol,pver)  ! ice single scattering albedo * tau
   real(kind_phys) :: ice_tau_w_g(nswbands,ncol,pver)  ! ice asymmetry parameter * tau * w
   real(kind_phys) :: snow_tau_r (nswbands,ncol,pver)  ! snow extinction optical depth
   real(kind_phys) :: snow_tau_w (nswbands,ncol,pver)  ! snow single scattering albedo * tau
   real(kind_phys) :: snow_tau_w_g(nswbands,ncol,pver) ! snow asymmetry parameter * tau * w
   real(kind_phys) :: cld_tau_w  (nswbands,ncol,pver)  ! cloud single scattering albedo * tau
   real(kind_phys) :: cld_tau_w_g(nswbands,ncol,pver)  ! cloud asymmetry parameter * w * tau
   real(kind_phys) :: grau_tau_r  (nswbands,ncol,pver) ! graupel extinction optical depth
   real(kind_phys) :: grau_tau_w  (nswbands,ncol,pver) ! graupel single scattering albedo * tau
   real(kind_phys) :: grau_tau_w_g(nswbands,ncol,pver) ! graupel asymmetry parameter * tau * w

   ! RRTMGP does not use this property in its 2-stream calculations.
   real(kind_phys) :: sw_tau_w_f(nswbands,ncol,pver) ! Forward scattered fraction * tau * w.

   character(len=*), parameter :: sub = 'rrtmgp_set_cloud_sw'
   !--------------------------------------------------------------------------------

   if (.not. dosw) then
      return
   end if

   ! Combine the cloud optical properties.

   select case (trim(liq_cld_optics))
   case ('slingo')
      ! Slingo (1989) liquid optics
      call slingo_liq_optics_sw(ncol, pver, nswbands, cld, rel, iclwpth, liq_tau, liq_tau_w, liq_tau_w_g, sw_tau_w_f, errmsg, errflg)
   case ('gammadist')
      ! gammadist liquid optics
      call get_liquid_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_liq, asm_sw_liq, ssa_sw_liq, lamc, pgam, g_lambda, g_mu, iclwpth, liq_tau, liq_tau_w, liq_tau_w_g, sw_tau_w_f, errmsg, errflg)
   case default
      write(errmsg,'(a,a)') sub, ': liq_cld_optics must be either slingo or gammadist'
      errflg = 1
   end select
   if (errflg /= 0) then
      return
   end if

   select case (trim(ice_cld_optics))
   case ('ebertcurry')
      ! Ebert and Curry (1992) ice optics
      call ec_ice_optics_sw(ncol, pver, nswbands, cld, rei, iciwpth, ice_tau, ice_tau_w, ice_tau_w_g, sw_tau_w_f, errmsg, errflg)
   case ('mitchell')
      ! Mitchell ice optics
      call interpolate_ice_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_ice, asm_sw_ice, ssa_sw_ice, iciwpth, dei, g_d_eff, &
           ice_tau, ice_tau_w, ice_tau_w_g, sw_tau_w_f, errmsg, errflg)
   case default
      write(errmsg,'(a,a)') sub, ': ice_cld_optics must be either ebertcurry or mitchell'
      errflg = 1
   end select
   if (errflg /= 0) then
      return
   end if

   do k = 1, pver
      do i = 1, ncol
         do ns = 1, nswbands
            cld_tau(ns,i,k)     =  liq_tau(ns,i,k)     + ice_tau(ns,i,k)
            cld_tau_w(ns,i,k)   =  liq_tau_w(ns,i,k)   + ice_tau_w(ns,i,k)
            cld_tau_w_g(ns,i,k) =  liq_tau_w_g(ns,i,k) + ice_tau_w_g(ns,i,k)
         end do
      end do
   end do

   ! snow
   if (do_snow) then
      call interpolate_ice_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_ice, asm_sw_ice, ssa_sw_ice, icswpth, des, g_d_eff, &
           snow_tau_r, snow_tau_w, snow_tau_w_g, sw_tau_w_f, errmsg, errflg)
      if (errflg /= 0) then
         return
      end if
   end if

   ! graupel
   if (do_graupel) then
      call get_grau_optics_sw(ncol, pver, nswbands, tiny_in, g_d_eff, ext_sw_ice, asm_sw_ice, ssa_sw_ice, iulog, icgrauwpth, &
           degrau, idx_sw_diag, grau_tau_r, grau_tau_w, grau_tau_w_g, sw_tau_w_f, ngrau_warn, errmsg, errflg)
      if (errflg /= 0) then
         return
      end if
   end if

   ! Combine with snow and graupel, and re-order the cloud optical properties
   ! from the RRTMG spectral bands (assumed in the optics datasets) to RRTMGP's:
   ! band ns of the outputs is band rrtmg_to_rrtmgp_swbands(ns) of the above.
   do k = 1, pver
      do i = 1, ncol
         do ns = 1, nswbands
            pb = rrtmg_to_rrtmgp_swbands(ns)
            ctau     = cld_tau(pb,i,k)
            ctau_w   = cld_tau_w(pb,i,k)
            ctau_w_g = cld_tau_w_g(pb,i,k)
            if (do_snow) then
               if (cldfprime(i,k) > 0._kind_phys) then
                  ctau     = ( cldfsnow(i,k)*snow_tau_r(pb,i,k) &
                             + cld(i,k)*ctau )/cldfprime(i,k)
                  ctau_w   = ( cldfsnow(i,k)*snow_tau_w(pb,i,k)  &
                             + cld(i,k)*ctau_w )/cldfprime(i,k)
                  ctau_w_g = ( cldfsnow(i,k)*snow_tau_w_g(pb,i,k) &
                             + cld(i,k)*ctau_w_g )/cldfprime(i,k)
               else
                  ctau     = 0._kind_phys
                  ctau_w   = 0._kind_phys
                  ctau_w_g = 0._kind_phys
               end if
               snow_tau(ns,i,k) = snow_tau_r(pb,i,k)
            else
               snow_tau(ns,i,k) = 0._kind_phys
            end if
            if (do_graupel) then
               if (cldfprime(i,k) > 0._kind_phys) then
                  ctau     = ( cldfgrau(i,k)*grau_tau_r(pb,i,k) &
                             + cld(i,k)*ctau )/cldfprime(i,k)
                  ctau_w   = ( cldfgrau(i,k)*grau_tau_w(pb,i,k)  &
                             + cld(i,k)*ctau_w )/cldfprime(i,k)
                  ctau_w_g = ( cldfgrau(i,k)*grau_tau_w_g(pb,i,k) &
                             + cld(i,k)*ctau_w_g )/cldfprime(i,k)
               else
                  ctau     = 0._kind_phys
                  ctau_w   = 0._kind_phys
                  ctau_w_g = 0._kind_phys
               end if
               grau_tau(ns,i,k) = grau_tau_r(pb,i,k)
            else
               grau_tau(ns,i,k) = 0._kind_phys
            end if
            c_cld_tau(ns,i,k)     = ctau
            c_cld_tau_w(ns,i,k)   = ctau_w
            c_cld_tau_w_g(ns,i,k) = ctau_w_g
         end do
      end do
   end do

   ! Set arrays for diagnostic output.
   ! cloud optical depth fields for the visible band (idx_sw_diag, RRTMGP order)
   pb = rrtmg_to_rrtmgp_swbands(idx_sw_diag)
   do k = 1, pver
      do i = 1, ncol
         tot_icld_vistau(i,k) = c_cld_tau(idx_sw_diag,i,k)
         liq_icld_vistau(i,k) = liq_tau(pb,i,k)
         ice_icld_vistau(i,k) = ice_tau(pb,i,k)
         snow_icld_vistau(i,k) = snow_tau(idx_sw_diag,i,k)
         grau_icld_vistau(i,k) = grau_tau(idx_sw_diag,i,k)
         ! multiply by total cloud fraction to get gridbox value
         tot_cld_vistau(i,k) = c_cld_tau(idx_sw_diag,i,k)*cldfprime(i,k)
      end do
   end do

   ! overwrite night columns with fillvalue
   do k = 1, pver
      do i = 1, nnite
         tot_cld_vistau(idxnite(i),k)   = fillvalue
         tot_icld_vistau(idxnite(i),k)  = fillvalue
         liq_icld_vistau(idxnite(i),k)  = fillvalue
         ice_icld_vistau(idxnite(i),k)  = fillvalue
         snow_icld_vistau(idxnite(i),k) = fillvalue
         grau_icld_vistau(idxnite(i),k) = fillvalue
      end do
   end do

end subroutine rrtmgp_sw_cloud_optics_run

!==============================================================================

subroutine get_grau_optics_sw(ncol, pver, nswbands, tiny_in, g_d_eff, ext_sw_ice, asm_sw_ice, ssa_sw_ice, &
                iulog, icgrauwpth, degrau, idx_sw_diag, tau, tau_w, tau_w_g, tau_w_f, nwarn, errmsg, errflg)
   use ccpp_kinds,                only: kind_phys

   integer, intent(in)  :: ncol
   integer, intent(in)  :: pver
   integer, intent(in)  :: nswbands
   integer, intent(in)  :: iulog
   integer, intent(in)  :: idx_sw_diag
   real(kind_phys), intent(in) :: tiny_in
   real(kind_phys), intent(in) :: ext_sw_ice(:,:)
   real(kind_phys), intent(in) :: asm_sw_ice(:,:)
   real(kind_phys), intent(in) :: ssa_sw_ice(:,:)
   real(kind_phys), intent(in) :: degrau(:,:)
   real(kind_phys), intent(in) :: g_d_eff(:)
   real(kind_phys), intent(in) :: icgrauwpth(:,:)

   real(kind_phys),intent(out) :: tau    (:,:,:) ! extinction optical depth
   real(kind_phys),intent(out) :: tau_w  (:,:,:) ! single scattering albedo * tau
   real(kind_phys),intent(out) :: tau_w_g(:,:,:) ! asymmetry parameter * tau * w
   real(kind_phys),intent(out) :: tau_w_f(:,:,:) ! forward scattered fraction * tau * w
   integer,         intent(out) :: nwarn         ! number of cells with tau > 100
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   integer :: i,k

   ! This does the same thing as get_ice_optics_sw, except with a different
   ! water path and effective diameter.
   call interpolate_ice_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_ice, asm_sw_ice, ssa_sw_ice, icgrauwpth, degrau, g_d_eff, tau, tau_w, &
        tau_w_g, tau_w_f, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if

   ! Count the cells to warn about; write the warnings only if there are any.
   nwarn = 0
   do k = 1, pver
      do i = 1, ncol
         if (tau(idx_sw_diag,i,k).gt.100._kind_phys) then
            nwarn = nwarn + 1
         end if
      enddo
   enddo
   if (nwarn > 0) then
      do i = 1, ncol
         do k = 1, pver
            if (tau(idx_sw_diag,i,k).gt.100._kind_phys) then
               write(iulog,*) 'WARNING: SW Graupel Tau > 100  (i,k,icgrauwpth,degrau,tau):'
               write(iulog,*) i,k,icgrauwpth(i,k), degrau(i,k), tau(idx_sw_diag,i,k)
            end if
         enddo
      enddo
   end if

end subroutine get_grau_optics_sw

!==============================================================================

subroutine get_liquid_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_liq, asm_sw_liq, ssa_sw_liq, lamc, pgam, g_lambda, &
                 g_mu, iclwpth, tau, tau_w, tau_w_g, tau_w_f, errmsg, errflg)
   use ccpp_kinds,                only: kind_phys
   use rrtmgp_cloud_optics_setup, only: nmu, nlambda
   use atmos_phys_interp_point,   only: interp_point_wgts, interp_point_wgts_blend, &
                                        interp_point_2d, interp_point_errmsg

   integer, intent(in)  :: ncol
   integer, intent(in)  :: pver
   integer, intent(in)  :: nswbands
   real(kind_phys), intent(in)  :: tiny_in
   real(kind_phys), intent(in)  :: g_lambda(:,:)
   real(kind_phys), intent(in)  :: g_mu(:)
   real(kind_phys), intent(in)  :: ext_sw_liq(:,:,:)
   real(kind_phys), intent(in)  :: asm_sw_liq(:,:,:)
   real(kind_phys), intent(in)  :: ssa_sw_liq(:,:,:)
   real(kind_phys), intent(in)  :: iclwpth(:,:)
   real(kind_phys), intent(in)  :: lamc(:,:)
   real(kind_phys), intent(in)  :: pgam(:,:)

   real(kind_phys), intent(out) :: tau    (:,:,:) ! extinction optical depth
   real(kind_phys), intent(out) :: tau_w  (:,:,:) ! single scattering albedo * tau
   real(kind_phys), intent(out) :: tau_w_g(:,:,:) ! asymmetry parameter * tau * w
   real(kind_phys), intent(out) :: tau_w_f(:,:,:) ! forward scattered fraction * tau * w
   character(len=*),   intent(out) :: errmsg
   integer,            intent(out) :: errflg

   integer :: i, k, swband
   integer :: mu_m, mu_p, lam_m, lam_p, ierr, ierr_max
   real(kind_phys) :: mu_ws, mu_wn, lam_ws, lam_wn
   real(kind_phys) :: ext, ssa, asm

   errmsg = ''
   errflg = 0

   ! For each cell: weights in mu (pgam), then in lambda (lamc) on the
   ! lambda row interpolated to that mu, then the bilinear table lookups
   ! (as get_mu_lambda_weights_ccpp and lininterp do).
   ierr_max = 0
   do k = 1,pver
      do i = 1,ncol
         ! lamc <= 0 seems to be clue from microphysics of no cloud
         if (lamc(i,k) > 0._kind_phys .and. .not. (iclwpth(i,k) < tiny_in)) then
            call interp_point_wgts(g_mu, nmu, pgam(i,k), mu_m, mu_p, mu_ws, mu_wn, ierr)
            if (ierr == 0) then
               call interp_point_wgts_blend(g_lambda, nlambda, mu_m, mu_p, mu_ws, mu_wn, &
                       lamc(i,k), lam_m, lam_p, lam_ws, lam_wn, ierr)
            end if
            ierr_max = max(ierr_max, ierr)
         else
            ierr = -1
         end if
         if (ierr == 0) then
            do swband = 1, nswbands
               ext = interp_point_2d(ext_sw_liq(:,:,swband), &
                       mu_m, mu_p, mu_ws, mu_wn, lam_m, lam_p, lam_ws, lam_wn)
               ssa = interp_point_2d(ssa_sw_liq(:,:,swband), &
                       mu_m, mu_p, mu_ws, mu_wn, lam_m, lam_p, lam_ws, lam_wn)
               asm = interp_point_2d(asm_sw_liq(:,:,swband), &
                       mu_m, mu_p, mu_ws, mu_wn, lam_m, lam_p, lam_ws, lam_wn)
               ! compute radiative properties
               tau(swband,i,k) = iclwpth(i,k) * ext
               tau_w(swband,i,k) = tau(swband,i,k) * ssa
               tau_w_g(swband,i,k) = tau_w(swband,i,k) * asm
               tau_w_f(swband,i,k) = tau_w_g(swband,i,k) * asm
            end do
         else
            do swband = 1, nswbands
               tau(swband,i,k) = 0._kind_phys
               tau_w(swband,i,k) = 0._kind_phys
               tau_w_g(swband,i,k) = 0._kind_phys
               tau_w_f(swband,i,k) = 0._kind_phys
            end do
         endif
      enddo
   enddo
   if (ierr_max /= 0) then
      errmsg = interp_point_errmsg(ierr_max)
      errflg = 1
   end if

end subroutine get_liquid_optics_sw

!==============================================================================

subroutine interpolate_ice_optics_sw(ncol, pver, nswbands, tiny_in, ext_sw_ice, asm_sw_ice, ssa_sw_ice, &
     iciwpth, dei, g_d_eff, tau, tau_w, tau_w_g, tau_w_f, errmsg, errflg)
  use ccpp_kinds,              only: kind_phys
  use atmos_phys_interp_point, only: interp_point_wgts, interp_point_1d, interp_point_errmsg

  integer, intent(in) :: ncol
  integer, intent(in) :: pver
  integer, intent(in) :: nswbands
  real(kind_phys), intent(in) :: tiny_in
  real(kind_phys), intent(in) :: iciwpth(:,:)
  real(kind_phys), intent(in) :: dei(:,:)
  real(kind_phys), intent(in) :: g_d_eff(:)
  real(kind_phys), intent(in) :: ext_sw_ice(:,:)
  real(kind_phys), intent(in) :: asm_sw_ice(:,:)
  real(kind_phys), intent(in) :: ssa_sw_ice(:,:)

  real(kind_phys),intent(out) :: tau    (:,:,:) ! extinction optical depth
  real(kind_phys),intent(out) :: tau_w  (:,:,:) ! single scattering albedo * tau
  real(kind_phys),intent(out) :: tau_w_g(:,:,:) ! asymmetry parameter * tau * w
  real(kind_phys),intent(out) :: tau_w_f(:,:,:) ! forward scattered fraction * tau * w
  character(len=*), intent(out) :: errmsg
  integer,          intent(out) :: errflg

  integer :: i, k, swband
  integer :: n_g_d
  integer :: jjm, jjp, ierr, ierr_max
  real(kind_phys) :: wgts, wgtn
  real(kind_phys) :: ext, ssa, asm

  errmsg = ''
  errflg = 0
  n_g_d = size(g_d_eff)

  ierr_max = 0
  do k = 1,pver
     do i = 1,ncol
        if( iciwpth(i,k) < tiny_in .or. dei(i,k) == 0._kind_phys) then
           ! if ice water path is too small, OD := 0
           ierr = -1
        else
           ! for each cell interpolate to find weights in g_d_eff grid.
           call interp_point_wgts(g_d_eff, n_g_d, dei(i,k), jjm, jjp, wgts, wgtn, ierr)
           ierr_max = max(ierr_max, ierr)
        endif
        if (ierr == 0) then
           ! interpolate into grid and extract radiative properties
           do swband = 1, nswbands
              ext = interp_point_1d(ext_sw_ice(:,swband), jjm, jjp, wgts, wgtn)
              ssa = interp_point_1d(ssa_sw_ice(:,swband), jjm, jjp, wgts, wgtn)
              asm = interp_point_1d(asm_sw_ice(:,swband), jjm, jjp, wgts, wgtn)
              tau    (swband,i,k) = iciwpth(i,k) * ext
              tau_w  (swband,i,k) = tau(swband,i,k) * ssa
              tau_w_g(swband,i,k) = tau_w(swband,i,k) * asm
              tau_w_f(swband,i,k) = tau_w_g(swband,i,k) * asm
           end do
        else
           do swband = 1, nswbands
              tau    (swband,i,k) = 0._kind_phys
              tau_w  (swband,i,k) = 0._kind_phys
              tau_w_g(swband,i,k) = 0._kind_phys
              tau_w_f(swband,i,k) = 0._kind_phys
           end do
        endif
     enddo
  enddo
  if (ierr_max /= 0) then
     errmsg = interp_point_errmsg(ierr_max)
     errflg = 1
  end if

end subroutine interpolate_ice_optics_sw

!==============================================================================

subroutine sw_band_slingo_index(nswbands, indxsl, errmsg, errflg)
   ! Index of the Slingo / Ebert-Curry coefficient set for each shortwave band,
   ! in RRTMG band order, from the band's upper wavelength boundary
   ! (A. Slingo (1989) equations 1-3):
   ! Use index 1 (0.25 to 0.69 micrometers) for visible
   ! Use index 2 (0.69 - 1.19 micrometers) for near-infrared
   ! Use index 3 (1.19 to 2.38 micrometers) for near-infrared
   ! Use index 4 (2.38 to 4.00 micrometers) for near-infrared
   use radiation_utils, only: get_sw_spectral_boundaries_ccpp
   use ccpp_kinds,      only: kind_phys

   integer,          intent(in)  :: nswbands
   integer,          intent(out) :: indxsl(nswbands)
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   real(kind_phys), dimension(nswbands) :: wavmin_gp  ! boundaries in RRTMGP band order
   real(kind_phys), dimension(nswbands) :: wavmax_gp  ! boundaries in RRTMGP band order
   real(kind_phys), dimension(nswbands) :: wavmin     ! boundaries in RRTMG band order
   real(kind_phys), dimension(nswbands) :: wavmax     ! boundaries in RRTMG band order
   integer :: ns

   ! get_sw_spectral_boundaries_ccpp returns the boundaries in RRTMGP band
   ! order; the optics are computed in RRTMG band order and re-ordered
   ! afterwards (see rrtmg_to_rrtmgp_swbands), so permute the boundaries
   ! into RRTMG band order here.
   call get_sw_spectral_boundaries_ccpp(wavmin_gp, wavmax_gp, 'microns', errmsg, errflg)
   if (errflg /= 0) then
      return
   end if
   wavmin(rrtmg_to_rrtmgp_swbands) = wavmin_gp
   wavmax(rrtmg_to_rrtmgp_swbands) = wavmax_gp

   do ns = 1, nswbands
      if(wavmax(ns) <= 0.7_kind_phys) then
         indxsl(ns) = 1
      else if(wavmax(ns) <= 1.25_kind_phys) then
         indxsl(ns) = 2
      else if(wavmax(ns) <= 2.38_kind_phys) then
         indxsl(ns) = 3
      else if(wavmax(ns) > 2.38_kind_phys) then
         indxsl(ns) = 4
      end if
   end do

end subroutine sw_band_slingo_index

!==============================================================================

subroutine slingo_liq_optics_sw(ncol, pver, nswbands, cldn, rel, iclwpth, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, errmsg, errflg)
   ! Slingo (1989) (https://doi.org/10.1175/1520-0469(1989)046<1419:AGPFTS>2.0.CO;2) shortwave liquid cloud optics,
   ! using the in-cloud liquid water path.
   use ccpp_kinds,      only: kind_phys

   integer, intent(in)  :: ncol
   integer, intent(in)  :: pver
   integer, intent(in)  :: nswbands
   real(kind_phys), intent(in) :: cldn(:,:)          ! cloud fraction [fraction]
   real(kind_phys), intent(in) :: rel(:,:)           ! liquid effective drop radius [um]
   real(kind_phys), intent(in) :: iclwpth(:,:)       ! in-cloud liquid water path [kg m-2]

   real(kind_phys), intent(out) :: liq_tau    (:,:,:) ! extinction optical depth
   real(kind_phys), intent(out) :: liq_tau_w  (:,:,:) ! single scattering albedo * tau
   real(kind_phys), intent(out) :: liq_tau_w_g(:,:,:) ! asymmetry parameter * tau * w
   real(kind_phys), intent(out) :: liq_tau_w_f(:,:,:) ! forward scattered fraction * tau * w
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   ! Minimum cloud amount (as a fraction of the grid-box area) to
   ! distinguish from clear sky
   real(kind_phys), parameter :: cldmin = 1.0e-80_kind_phys

   ! Decimal precision of cloud amount (0 -> preserve full resolution;
   ! 10^-n -> preserve n digits of cloud amount)
   real(kind_phys), parameter :: cldeps = 0.0_kind_phys

   ! A. Slingo's data for cloud particle radiative properties (from 'A GCM
   ! Parameterization for the Shortwave Properties of Water Clouds' JAS
   ! vol. 46 may 1989 pp 1419-1427)
   real(kind_phys), parameter :: abarl(4) = &  ! A coefficient for extinction optical depth
      [ 2.817e-02_kind_phys, 2.682e-02_kind_phys,2.264e-02_kind_phys,1.281e-02_kind_phys]
   real(kind_phys), parameter :: bbarl(4) = &  ! B coefficient for extinction optical depth
      [ 1.305_kind_phys    , 1.346_kind_phys    ,1.454_kind_phys    ,1.641_kind_phys    ]
   real(kind_phys), parameter :: cbarl(4) = &  ! C coefficient for single scat albedo
      [-5.62e-08_kind_phys ,-6.94e-06_kind_phys ,4.64e-04_kind_phys ,0.201_kind_phys    ]
   real(kind_phys), parameter :: dbarl(4) = &  ! D coefficient for single  scat albedo
      [ 1.63e-07_kind_phys , 2.35e-05_kind_phys ,1.24e-03_kind_phys ,7.56e-03_kind_phys ]
   real(kind_phys), parameter :: ebarl(4) = &  ! E coefficient for asymmetry parameter
      [ 0.829_kind_phys    , 0.794_kind_phys    ,0.754_kind_phys    ,0.826_kind_phys    ]
   real(kind_phys), parameter :: fbarl(4) = &  ! F coefficient for asymmetry parameter
      [ 2.482e-03_kind_phys, 4.226e-03_kind_phys,6.560e-03_kind_phys,4.353e-03_kind_phys]

   real(kind_phys) :: abarli        ! A coefficient for current spectral band
   real(kind_phys) :: bbarli        ! B coefficient for current spectral band
   real(kind_phys) :: cbarli        ! C coefficient for current spectral band
   real(kind_phys) :: dbarli        ! D coefficient for current spectral band
   real(kind_phys) :: ebarli        ! E coefficient for current spectral band
   real(kind_phys) :: fbarli        ! F coefficient for current spectral band

   ! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
   ! greater than 20 micro-meters

   integer :: ns, i, k
   integer :: indxsl(nswbands)
   real(kind_phys) :: tmp1l, tmp2l, tmp3l, g

   ! Set error variables
   errmsg = ''
   errflg = 0

   call sw_band_slingo_index(nswbands, indxsl, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if

   do k=1,pver
      do i=1,ncol
         do ns = 1, nswbands
            ! Set cloud extinction optical depth, single scatter albedo,
            ! asymmetry parameter, and forward scattered fraction:
            abarli = abarl(indxsl(ns))
            bbarli = bbarl(indxsl(ns))
            cbarli = cbarl(indxsl(ns))
            dbarli = dbarl(indxsl(ns))
            ebarli = ebarl(indxsl(ns))
            fbarli = fbarl(indxsl(ns))

            ! note that optical properties for liquid valid only
            ! in range of 4.2 > rel > 16 micron (Slingo 89)
            if (cldn(i,k) >= cldmin .and. cldn(i,k) >= cldeps) then
               tmp1l = abarli + bbarli/min(max(4.2_kind_phys,rel(i,k)),16._kind_phys)
               liq_tau(ns,i,k) = 1000._kind_phys*iclwpth(i,k)*tmp1l
            else
               liq_tau(ns,i,k) = 0.0_kind_phys
            end if

            tmp2l = 1._kind_phys - cbarli - dbarli*min(max(4.2_kind_phys,rel(i,k)),16._kind_phys)
            tmp3l = fbarli*min(max(4.2_kind_phys,rel(i,k)),16._kind_phys)
            ! Do not let single scatter albedo be 1.  Delta-eddington solution
            ! for non-conservative case has different analytic form from solution
            ! for conservative case, and raddedmx is written for non-conservative case.
            liq_tau_w(ns,i,k) = liq_tau(ns,i,k) * min(tmp2l,0.999999_kind_phys)
            g = ebarli + tmp3l
            liq_tau_w_g(ns,i,k) = liq_tau_w(ns,i,k) * g
            liq_tau_w_f(ns,i,k) = liq_tau_w(ns,i,k) * g * g

         end do ! nswbands
      end do    ! End do i=1,ncol
   end do       ! End do k=1,pver

end subroutine slingo_liq_optics_sw

!==============================================================================

subroutine ec_ice_optics_sw(ncol, pver, nswbands, cldn, rei, iciwpth, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, errmsg, errflg)
   ! Ebert and Curry (1992) (https://doi.org/10.1029/91JD02472) shortwave ice cloud optics,
   ! using the in-cloud ice water path.
   use ccpp_kinds,      only: kind_phys

   integer, intent(in)  :: ncol
   integer, intent(in)  :: pver
   integer, intent(in)  :: nswbands
   real(kind_phys), intent(in) :: cldn(:,:)          ! cloud fraction [fraction]
   real(kind_phys), intent(in) :: rei(:,:)           ! ice effective drop size [um]
   real(kind_phys), intent(in) :: iciwpth(:,:)       ! in-cloud ice water path [kg m-2]

   real(kind_phys), intent(out) :: ice_tau    (:,:,:) ! extinction optical depth
   real(kind_phys), intent(out) :: ice_tau_w  (:,:,:) ! single scattering albedo * tau
   real(kind_phys), intent(out) :: ice_tau_w_g(:,:,:) ! asymmetry parameter * tau * w
   real(kind_phys), intent(out) :: ice_tau_w_f(:,:,:) ! forward scattered fraction * tau * w
   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   real(kind_phys), parameter :: scalefactor = 1._kind_phys

   ! Minimum cloud amount (as a fraction of the grid-box area) to
   ! distinguish from clear sky
   real(kind_phys), parameter :: cldmin = 1.0e-80_kind_phys

   ! Decimal precision of cloud amount (0 -> preserve full resolution;
   ! 10^-n -> preserve n digits of cloud amount)
   real(kind_phys), parameter :: cldeps = 0.0_kind_phys

   ! ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
   real(kind_phys), parameter :: abari(4) = &     ! a coefficient for extinction optical depth
      [ 3.448e-03_kind_phys, 3.448e-03_kind_phys,3.448e-03_kind_phys,3.448e-03_kind_phys]
   real(kind_phys), parameter :: bbari(4) = &     ! b coefficient for extinction optical depth
      [ 2.431_kind_phys    , 2.431_kind_phys    ,2.431_kind_phys    ,2.431_kind_phys    ]
   real(kind_phys), parameter :: cbari(4) = &     ! c coefficient for single scat albedo
      [ 1.00e-05_kind_phys , 1.10e-04_kind_phys ,1.861e-02_kind_phys,.46658_kind_phys   ]
   real(kind_phys), parameter :: dbari(4) = &     ! d coefficient for single scat albedo
      [ 0.0_kind_phys      , 1.405e-05_kind_phys,8.328e-04_kind_phys,2.05e-05_kind_phys ]
   real(kind_phys), parameter :: ebari(4) = &     ! e coefficient for asymmetry parameter
      [ 0.7661_kind_phys   , 0.7730_kind_phys   ,0.794_kind_phys    ,0.9595_kind_phys   ]
   real(kind_phys), parameter :: fbari(4) = &     ! f coefficient for asymmetry parameter
      [ 5.851e-04_kind_phys, 5.665e-04_kind_phys,7.267e-04_kind_phys,1.076e-04_kind_phys]

   real(kind_phys) :: abarii           ! A coefficient for current spectral band
   real(kind_phys) :: bbarii           ! B coefficient for current spectral band
   real(kind_phys) :: cbarii           ! C coefficient for current spectral band
   real(kind_phys) :: dbarii           ! D coefficient for current spectral band
   real(kind_phys) :: ebarii           ! E coefficient for current spectral band
   real(kind_phys) :: fbarii           ! F coefficient for current spectral band

   integer :: ns, i, k
   integer :: indxsl(nswbands)
   real(kind_phys) :: tmp1i, tmp2i, tmp3i, g

   ! Set error variables
   errmsg = ''
   errflg = 0

   call sw_band_slingo_index(nswbands, indxsl, errmsg, errflg)
   if (errflg /= 0) then
      return
   end if

   do k=1,pver
      do i=1,ncol
         do ns = 1, nswbands
            abarii = abari(indxsl(ns))
            bbarii = bbari(indxsl(ns))
            cbarii = cbari(indxsl(ns))
            dbarii = dbari(indxsl(ns))
            ebarii = ebari(indxsl(ns))
            fbarii = fbari(indxsl(ns))

            ! note that optical properties for ice valid only
            ! in range of 13 > rei > 130 micron (Ebert and Curry 92)
            if (cldn(i,k) >= cldmin .and. cldn(i,k) >= cldeps) then
               tmp1i = abarii + bbarii/max(13._kind_phys,min(scalefactor*rei(i,k),130._kind_phys))
               ice_tau(ns,i,k) = 1000.0_kind_phys*iciwpth(i,k)*tmp1i
            else
               ice_tau(ns,i,k) = 0.0_kind_phys
            end if

            tmp2i = 1._kind_phys - cbarii - dbarii*min(max(13._kind_phys,scalefactor*rei(i,k)),130._kind_phys)
            tmp3i = fbarii*min(max(13._kind_phys,scalefactor*rei(i,k)),130._kind_phys)
            ! Do not let single scatter albedo be 1.  Delta-eddington solution
            ! for non-conservative case has different analytic form from solution
            ! for conservative case, and raddedmx is written for non-conservative case.
            ice_tau_w(ns,i,k) = ice_tau(ns,i,k) * min(tmp2i,0.999999_kind_phys)
            g = ebarii + tmp3i
            ice_tau_w_g(ns,i,k) = ice_tau_w(ns,i,k) * g
            ice_tau_w_f(ns,i,k) = ice_tau_w(ns,i,k) * g * g

         end do ! nswbands
      end do    ! End do i=1,ncol
   end do       ! End do k=1,pver

end subroutine ec_ice_optics_sw

end module rrtmgp_sw_cloud_optics
