! zeros input fluxes to check_energy
! before running other schemes
module check_energy_zero_fluxes
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: check_energy_zero_fluxes_run

contains

!> \section arg_table_check_energy_zero_fluxes_run Argument Table
!! \htmlinclude arg_table_check_energy_zero_fluxes_run.html
  subroutine check_energy_zero_fluxes_run(ncol, name, flx_vap, flx_cnd, flx_ice, flx_sen, errmsg, errflg)
    ! Input arguments
    integer,            intent(in)     :: ncol           ! number of atmospheric columns

    ! Output arguments
    character(len=64),  intent(out)    :: name           ! parameterization name for fluxes
    real(kind_phys),    intent(out)    :: flx_vap(:)     ! boundary flux of vapor [kg m-2 s-1]
    real(kind_phys),    intent(out)    :: flx_cnd(:)     ! boundary flux of liquid+ice (precip?) [m s-1]
    real(kind_phys),    intent(out)    :: flx_ice(:)     ! boundary flux of ice (snow?) [m s-1]
    real(kind_phys),    intent(out)    :: flx_sen(:)     ! boundary flux of sensible heat [W m-2]
    character(len=512), intent(out)    :: errmsg         ! error message
    integer,            intent(out)    :: errflg         ! error flag

    integer :: i

    errmsg = ''
    errflg = 0

    ! reset values to zero (the fluxes are (horizontal_loop_extent) = ncol)
    name = ''
    !$acc parallel loop copy(flx_vap, flx_cnd, flx_ice, flx_sen)
    do i = 1, ncol
      flx_vap(i) = 0._kind_phys
      flx_cnd(i) = 0._kind_phys
      flx_ice(i) = 0._kind_phys
      flx_sen(i) = 0._kind_phys
    end do
  end subroutine check_energy_zero_fluxes_run

end module check_energy_zero_fluxes
