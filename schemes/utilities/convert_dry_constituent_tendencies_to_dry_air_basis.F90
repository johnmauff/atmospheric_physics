! Convert dry constituent tendencies from moist to dry air basis
module convert_dry_constituent_tendencies_to_dry_air_basis
  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  ! public CCPP-compliant subroutines
  public :: convert_dry_constituent_tendencies_to_dry_air_basis_run

contains

!> \section arg_table_convert_dry_constituent_tendencies_to_dry_air_basis_run Argument Table
!! \htmlinclude convert_dry_constituent_tendencies_to_dry_air_basis_run.html
  subroutine convert_dry_constituent_tendencies_to_dry_air_basis_run( &
    ncol, pver, pcnst, &
    pdel, pdeldry, &
    const_props, &
    tend_q, &
    errmsg, errflg)

    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    ! Input arguments
    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: pcnst
    real(kind_phys),  intent(in)    :: pdel(:, :)          ! Layer thickness (moist air) [Pa]
    real(kind_phys),  intent(in)    :: pdeldry(:, :)       ! Layer thickness (dry air) [Pa]
    
    ! Framework dependency for constituent properties
    type(ccpp_constituent_prop_ptr_t), &
                      intent(in)    :: const_props(:)      ! CCPP constituent properties pointer

    ! Input/Output arguments
    real(kind_phys),  intent(inout) :: tend_q(:, :, :)     ! Constituent tendencies [kg kg-1 s-1]

    ! Output arguments
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k, m
    logical :: const_is_dry(pcnst)

    errmsg = ''
    errflg = 0

    ! Which constituents are dry type, on the host: the constituent
    ! properties (a pointer derived type) never go to the device. An error
    ! now returns before any tendency is converted (before, constituents
    ! ahead of the failing one were already converted).
    do m = 1, pcnst
       call const_props(m)%is_dry(const_is_dry(m), errflg, errmsg)
       if (errflg /= 0) return
    end do

    ! Convert the tendencies for the dry constituents to dry air basis
    !$acc parallel loop collapse(3) copyin(const_is_dry, pdel, pdeldry) copy(tend_q)
    do m = 1, pcnst
       do k = 1, pver
          do i = 1, ncol
             if (const_is_dry(m)) then
                tend_q(i, k, m) = tend_q(i, k, m) * pdel(i, k) / pdeldry(i, k)
             end if
          end do
       end do
    end do

  end subroutine convert_dry_constituent_tendencies_to_dry_air_basis_run

end module convert_dry_constituent_tendencies_to_dry_air_basis
