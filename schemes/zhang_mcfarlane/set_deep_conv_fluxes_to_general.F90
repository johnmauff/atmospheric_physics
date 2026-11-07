module set_deep_conv_fluxes_to_general

   use ccpp_kinds, only:  kind_phys

   implicit none

contains

!===============================================================================
!> \section arg_table_set_deep_conv_fluxes_to_general_run Argument Table
!! \htmlinclude set_deep_conv_fluxes_to_general_run.html
!!

   subroutine set_deep_conv_fluxes_to_general_run(prec_gen, prec_dp, prdprec_gen, prdprec_dp)

   real(kind_phys), intent(out) :: prec_gen(:)
   real(kind_phys), intent(in) :: prec_dp(:)
   real(kind_phys), intent(out) :: prdprec_gen(:,:)
   real(kind_phys), intent(in) :: prdprec_dp(:,:)

   integer :: i, k

   ! Copies, with the outputs' extents (the arrays are conformable).
   !$acc parallel loop copyin(prec_dp) copy(prec_gen)
   do i = 1, size(prec_gen, 1)
      prec_gen(i) = prec_dp(i)
   end do
   !$acc parallel loop collapse(2) copyin(prdprec_dp) copy(prdprec_gen)
   do k = 1, size(prdprec_gen, 2)
      do i = 1, size(prdprec_gen, 1)
         prdprec_gen(i,k) = prdprec_dp(i,k)
      end do
   end do

   end subroutine set_deep_conv_fluxes_to_general_run

end module set_deep_conv_fluxes_to_general
