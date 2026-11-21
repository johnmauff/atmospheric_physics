! Apply minimum, maximum, and scaling to the subgrid vertical velocity.
! wsub is the subgrid vertical velocity for droplet (liquid) nucleation.
! wsubi is the subgrid vertical velocity for ice nucleation,
! derived here from the raw wsub before wsub is scaled
! (and capped at 0.2 m s-1 when preexisting ice is not used).
!
! OpenACC: one data region; the fills above top_lev and the scaling are
! pointwise kernels, each value computed in a scalar and stored once.
module scale_subgrid_vertical_velocity

  implicit none
  private

  public :: scale_subgrid_vertical_velocity_run

contains

!> \section arg_table_scale_subgrid_vertical_velocity_run Argument Table
!! \htmlinclude scale_subgrid_vertical_velocity_run.html
  pure subroutine scale_subgrid_vertical_velocity_run( &
    ncol, pver, top_lev,                          &
    wsub_min, wsubi_min,                          &
    wsub_scale, wsubi_scale,                      &
    use_preexisting_ice,                          &
    wsub, wsubi,                                  &
    errmsg, errflg)

    use ccpp_kinds, only: kind_phys

    ! Input arguments
    integer,          intent(in)    :: ncol
    integer,          intent(in)    :: pver
    integer,          intent(in)    :: top_lev             ! top vertical level for cloud physics

    real(kind_phys),  intent(in)    :: wsub_min            ! minimum subgrid vertical velocity for droplet nucleation [m s-1]
    real(kind_phys),  intent(in)    :: wsubi_min           ! minimum subgrid vertical velocity for ice nucleation [m s-1]
    real(kind_phys),  intent(in)    :: wsub_scale          ! scaling factor for subgrid vertical velocity (liquid)
    real(kind_phys),  intent(in)    :: wsubi_scale         ! scaling factor for subgrid vertical velocity (ice)

    logical,          intent(in)    :: use_preexisting_ice ! use preexisting ice? [flag]

    ! Input/Output arguments
    real(kind_phys),  intent(inout) :: wsub(:, :)          ! subgrid vertical velocity [m s-1] (ncol, pver)

    ! Output arguments
    real(kind_phys),  intent(out)   :: wsubi(:, :)         ! subgrid vertical velocity for ice nucleation [m s-1] (ncol, pver)

    character(len=*),   intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i, k
    real(kind_phys) :: wsubi_ik

    errmsg = ''
    errflg = 0

    !$acc data copy(wsub) copyout(wsubi)

    ! Set minimum values above top_lev
    !$acc parallel loop gang vector collapse(2)
    do k = 1, top_lev-1
      do i = 1, size(wsub, 1)
        wsub(i,k) = wsub_min
      end do
    end do
    !$acc parallel loop gang vector collapse(2)
    do k = 1, top_lev-1
      do i = 1, size(wsubi, 1)
        wsubi(i,k) = wsubi_min
      end do
    end do

    !$acc parallel loop gang vector collapse(2) private(wsubi_ik)
    do k = top_lev, pver
      do i = 1, ncol
        ! Derive wsubi from raw wsub BEFORE scaling wsub
        wsubi_ik = max(wsubi_min, wsub(i,k)) * wsubi_scale
        if (.not. use_preexisting_ice) then
          wsubi_ik = min(wsubi_ik, 0.2_kind_phys)
        end if
        wsubi(i,k) = wsubi_ik

        wsub(i,k) = max(wsub_min, wsub(i,k)) * wsub_scale
      end do
    end do

    !$acc end data

  end subroutine scale_subgrid_vertical_velocity_run

end module scale_subgrid_vertical_velocity
