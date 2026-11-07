module set_general_conv_fluxes_to_deep

   use ccpp_kinds, only:  kind_phys

   implicit none

contains

!===============================================================================
!> \section arg_table_set_general_conv_fluxes_to_deep_run Argument Table
!! \htmlinclude set_general_conv_fluxes_to_deep_run.html
!!

   subroutine set_general_conv_fluxes_to_deep_run(tend_s_snwprd_gen, tend_s_snwprd_dp, tend_s_snwevmlt_gen, tend_s_snwevmlt_dp, &
          prec_gen, prec_dp, snow_gen, snow_dp, ntprprd_gen, ntprprd_dp, ntsnprd_gen, ntsnprd_dp, &
          flxprec_gen, flxprec_dp, flxsnow_gen, flxsnow_dp)

   real(kind_phys), intent(in) :: tend_s_snwprd_gen(:,:)
   real(kind_phys), intent(out) :: tend_s_snwprd_dp(:,:)
   real(kind_phys), intent(in) :: tend_s_snwevmlt_gen(:,:)
   real(kind_phys), intent(out) :: tend_s_snwevmlt_dp(:,:)
   real(kind_phys), intent(in) :: prec_gen(:)
   real(kind_phys), intent(out) :: prec_dp(:)
   real(kind_phys), intent(in) :: snow_gen(:)
   real(kind_phys), intent(out) :: snow_dp(:)
   real(kind_phys), intent(in) :: ntprprd_gen(:,:)
   real(kind_phys), intent(out) :: ntprprd_dp(:,:)
   real(kind_phys), intent(in) :: ntsnprd_gen(:,:)
   real(kind_phys), intent(out) :: ntsnprd_dp(:,:)
   real(kind_phys), intent(in) :: flxprec_gen(:,:)
   real(kind_phys), intent(out) :: flxprec_dp(:,:)
   real(kind_phys), intent(in) :: flxsnow_gen(:,:)
   real(kind_phys), intent(out) :: flxsnow_dp(:,:)

   integer :: i, k

   ! Copies, with the outputs' extents (the arrays are conformable):
   ! layer fields, interface fields, then column fields.
   !$acc parallel loop collapse(2) &
   !$acc copyin(tend_s_snwprd_gen, tend_s_snwevmlt_gen, ntprprd_gen, ntsnprd_gen) &
   !$acc copy(tend_s_snwprd_dp, tend_s_snwevmlt_dp, ntprprd_dp, ntsnprd_dp)
   do k = 1, size(tend_s_snwprd_dp, 2)
      do i = 1, size(tend_s_snwprd_dp, 1)
         tend_s_snwprd_dp(i,k) = tend_s_snwprd_gen(i,k)
         tend_s_snwevmlt_dp(i,k) = tend_s_snwevmlt_gen(i,k)
         ntprprd_dp(i,k) = ntprprd_gen(i,k)
         ntsnprd_dp(i,k) = ntsnprd_gen(i,k)
      end do
   end do

   !$acc parallel loop collapse(2) copyin(flxprec_gen, flxsnow_gen) copy(flxprec_dp, flxsnow_dp)
   do k = 1, size(flxprec_dp, 2)
      do i = 1, size(flxprec_dp, 1)
         flxprec_dp(i,k) = flxprec_gen(i,k)
         flxsnow_dp(i,k) = flxsnow_gen(i,k)
      end do
   end do

   !$acc parallel loop copyin(prec_gen, snow_gen) copy(prec_dp, snow_dp)
   do i = 1, size(prec_dp, 1)
      prec_dp(i) = prec_gen(i)
      snow_dp(i) = snow_gen(i)
   end do

   end subroutine set_general_conv_fluxes_to_deep_run

end module set_general_conv_fluxes_to_deep
