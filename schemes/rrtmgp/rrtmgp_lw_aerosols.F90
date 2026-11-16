!> \file rrtmgp_lw_aerosols.F90
!!

!> This module sets the RRTMGP aerosol longwave optical depth
module rrtmgp_lw_aerosols
  implicit none
  private

  public rrtmgp_lw_aerosols_run
contains

!> \section arg_table_rrtmgp_lw_aerosols_run Argument Table
!! \htmlinclude rrtmgp_lw_aerosols_run.html
!!
   subroutine rrtmgp_lw_aerosols_run( &
     dolwrad, ncol, ktopcam, ktoprad, &
     aer_lw_abs, &
     aer_lw, errmsg, errflg)

    use ccpp_optical_props, only: ty_optical_props_1scl_ccpp
    use ccpp_kinds,         only: kind_phys

    ! Inputs
    logical,          intent(in)    :: dolwrad                    !< do longwave calculation
    integer,          intent(in)    :: ncol                       !< number of columns
    integer,          intent(in)    :: ktopcam                    !< top CAM level for RRTMGP
    integer,          intent(in)    :: ktoprad                    !< top RRTMGP level

    real(kind_phys),  intent(in)    :: aer_lw_abs(:, :, :)        !< LW absorption OD (ncol,pver,nlwbands)

    ! Outputs
    class(ty_optical_props_1scl_ccpp), intent(inout) :: aer_lw    !< Aerosol optical properties object

    character(len=*), intent(out) :: errmsg                       !< CCPP error message
    integer,          intent(out) :: errflg                       !< CCPP error flag

    integer :: i, k, ib

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. dolwrad) return

    ! The object is only copied in (copying it back could overwrite scalars
    ! set on the host); its tau array is written out.
    !$acc data copyin(aer_lw_abs, aer_lw%optical_props) copyout(aer_lw%optical_props%tau)
    ! Initialize to zero (provides zero optical depth in any extra layer)
    !$acc parallel loop collapse(3)
    do ib = 1, size(aer_lw%optical_props%tau, 3)
      do k = 1, size(aer_lw%optical_props%tau, 2)
        do i = 1, size(aer_lw%optical_props%tau, 1)
          aer_lw%optical_props%tau(i, k, ib) = 0.0_kind_phys
        end do
      end do
    end do

    ! Map CAM vertical levels to RRTMGP vertical levels
    ! No band reorder needed for LW
    !$acc parallel loop collapse(3)
    do ib = 1, size(aer_lw_abs, 3)
      do k = ktoprad, size(aer_lw%optical_props%tau, 2)
        do i = 1, ncol
          aer_lw%optical_props%tau(i, k, ib) = aer_lw_abs(i, k - ktoprad + ktopcam, ib)
        end do
      end do
    end do
    !$acc end data

  end subroutine rrtmgp_lw_aerosols_run
end module rrtmgp_lw_aerosols
