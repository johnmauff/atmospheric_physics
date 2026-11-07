module set_general_conv_fluxes_to_shallow

   use ccpp_kinds, only:  kind_phys

   implicit none

contains
!> \section arg_table_set_general_conv_fluxes_to_shallow_run Argument Table
!! \htmlinclude set_general_conv_fluxes_to_shallow_run.html
   subroutine set_general_conv_fluxes_to_shallow_run(tend_s_snwprd_gen, tend_s_snwprd_sh, tend_s_snwevmlt_gen, tend_s_snwevmlt_sh, &
          prec_gen, prec_sh, snow_gen, snow_sh, ntprprd_gen, ntprprd_sh, ntsnprd_gen, ntsnprd_sh, &
          flxprec_gen, flxprec_sh, flxsnow_gen, flxsnow_sh, &
          flx_ice)

   real(kind_phys), intent(in) :: tend_s_snwprd_gen(:,:)
   real(kind_phys), intent(out) :: tend_s_snwprd_sh(:,:)
   real(kind_phys), intent(in) :: tend_s_snwevmlt_gen(:,:)
   real(kind_phys), intent(out) :: tend_s_snwevmlt_sh(:,:)
   real(kind_phys), intent(in) :: prec_gen(:) ! precipitation rate is modified by zm_conv_evap based on computed flux
   real(kind_phys), intent(out) :: prec_sh(:) ! and thus needs to be renamed from generic to shallow.
   real(kind_phys), intent(in) :: snow_gen(:)
   real(kind_phys), intent(out) :: snow_sh(:)
   real(kind_phys), intent(in) :: ntprprd_gen(:,:)
   real(kind_phys), intent(out) :: ntprprd_sh(:,:)
   real(kind_phys), intent(in) :: ntsnprd_gen(:,:)
   real(kind_phys), intent(out) :: ntsnprd_sh(:,:)
   real(kind_phys), intent(in) :: flxprec_gen(:,:)
   real(kind_phys), intent(out) :: flxprec_sh(:,:)
   real(kind_phys), intent(in) :: flxsnow_gen(:,:)
   real(kind_phys), intent(out) :: flxsnow_sh(:,:)

   real(kind_phys), intent(out)    :: flx_ice(:)         ! net_lwe_ice_fluxes_through_top_and_bottom_of_atmosphere_column [m s-1] for check_energy_chng

   integer :: i, k

   ! Copies, with the outputs' extents (the arrays are conformable):
   ! layer fields, interface fields, then column fields.
   !$acc parallel loop collapse(2) &
   !$acc copyin(tend_s_snwprd_gen, tend_s_snwevmlt_gen, ntprprd_gen, ntsnprd_gen) &
   !$acc copy(tend_s_snwprd_sh, tend_s_snwevmlt_sh, ntprprd_sh, ntsnprd_sh)
   do k = 1, size(tend_s_snwprd_sh, 2)
      do i = 1, size(tend_s_snwprd_sh, 1)
         tend_s_snwprd_sh(i,k) = tend_s_snwprd_gen(i,k)
         tend_s_snwevmlt_sh(i,k) = tend_s_snwevmlt_gen(i,k)
         ntprprd_sh(i,k) = ntprprd_gen(i,k)
         ntsnprd_sh(i,k) = ntsnprd_gen(i,k)
      end do
   end do

   !$acc parallel loop collapse(2) copyin(flxprec_gen, flxsnow_gen) copy(flxprec_sh, flxsnow_sh)
   do k = 1, size(flxprec_sh, 2)
      do i = 1, size(flxprec_sh, 1)
         flxprec_sh(i,k) = flxprec_gen(i,k)
         flxsnow_sh(i,k) = flxsnow_gen(i,k)
      end do
   end do

   ! flx_ice: boundary flux for check energy (= snow_sh)
   !$acc parallel loop copyin(prec_gen, snow_gen) copy(prec_sh, snow_sh, flx_ice)
   do i = 1, size(prec_sh, 1)
      prec_sh(i) = prec_gen(i)
      snow_sh(i) = snow_gen(i)
      flx_ice(i) = snow_gen(i)
   end do

   end subroutine set_general_conv_fluxes_to_shallow_run

end module set_general_conv_fluxes_to_shallow
