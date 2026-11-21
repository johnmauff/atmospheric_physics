! Set default dust radii for 4 size bins used in contact freezing.
! For modal aerosols, bin 3 is overwritten later with wet diameter.
!
! OpenACC: one fill kernel.
module dust_default_radii

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: dust_default_radii_timestep_init

  real(kind_phys), parameter :: rn_dst1 = 0.258e-6_kind_phys
  real(kind_phys), parameter :: rn_dst2 = 0.717e-6_kind_phys
  real(kind_phys), parameter :: rn_dst3 = 1.576e-6_kind_phys
  real(kind_phys), parameter :: rn_dst4 = 3.026e-6_kind_phys

contains

!> \section arg_table_dust_default_radii_timestep_init Argument Table
!! \htmlinclude dust_default_radii_timestep_init.html
  pure subroutine dust_default_radii_timestep_init(ncol, pver, ndust, rndst, errmsg, errflg)
    integer,            intent(in)  :: ncol
    integer,            intent(in)  :: pver
    integer,            intent(in)  :: ndust

    real(kind_phys),    intent(out) :: rndst(:,:,:)
    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    integer :: i, k

    errmsg = ''
    errflg = 0

    !$acc parallel loop gang vector collapse(2) copyout(rndst)
    do k = 1, pver
      do i = 1, ncol
        rndst(i, k, 1) = rn_dst1
        rndst(i, k, 2) = rn_dst2
        rndst(i, k, 3) = rn_dst3
        rndst(i, k, 4) = rn_dst4
      end do
    end do

  end subroutine dust_default_radii_timestep_init

end module dust_default_radii
