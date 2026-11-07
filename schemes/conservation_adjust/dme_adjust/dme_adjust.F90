module dme_adjust

   use ccpp_kinds, only:  kind_phys

   implicit none

contains
!===============================================================================
!> \section arg_table_dme_adjust_run Argument Table
!! \htmlinclude dme_adjust_run.html
!!
  subroutine dme_adjust_run(ncol, pver, pcnst, ps, pint, pdel, lnpint, rpdel, const_props, const_array, qini, liqini, iceini, &
                            is_dycore_moist, errmsg, errflg)
    !-----------------------------------------------------------------------
    !
    ! Purpose: Adjust the dry mass in each layer back to the value of physics input state
    !
    ! Method: Conserve the integrated mass, momentum and total energy in each layer
    !         by scaling the specific mass of consituents, specific momentum (velocity)
    !         and specific total energy by the relative change in layer mass. Solve for
    !         the new temperature by subtracting the new kinetic energy from total energy
    !         and inverting the hydrostatic equation
    !
    !         The mass in each layer is modified, changing the relationship of the layer
    !         interfaces and midpoints to the surface pressure. The result is no longer in
    !         the original hybrid coordinate.
    !
    !         This procedure cannot be applied to the "eul" or "sld" dycores because they
    !         require the hybrid coordinate.
    !
    ! Author: Byron Boville

    ! !REVISION HISTORY:
    !   03.03.28  Boville    Created, partly from code by Lin in p_d_adjust
    !
    !-----------------------------------------------------------------------

    use ccpp_constituent_prop_mod, only: ccpp_constituent_prop_ptr_t

    implicit none
    !
    ! Arguments
    !
    integer,                           intent(in)    :: ncol
    integer,                           intent(in)    :: pver
    integer,                           intent(in)    :: pcnst
    real(kind_phys),                   intent(inout) :: ps(:)
    real(kind_phys),                   intent(inout) :: pint(:,:)
    real(kind_phys),                   intent(inout) :: pdel(:,:)
    real(kind_phys),                   intent(inout) :: lnpint(:,:)
    real(kind_phys),                   intent(inout) :: rpdel(:,:)
    type(ccpp_constituent_prop_ptr_t), intent(in)    :: const_props(:)
    real(kind_phys),                   intent(inout) :: const_array(:,:,:)

    real(kind_phys),                   intent(in)    :: qini(:,:)    ! initial specific humidity
    real(kind_phys),                   intent(in)    :: liqini(:,:)  ! initial total liquid
    real(kind_phys),                   intent(in)    :: iceini(:,:)  ! initial total ice
    logical,                           intent(in)    :: is_dycore_moist
    character(len=512),                intent(out)   :: errmsg
    integer,                           intent(out)   :: errflg

    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i,k,m         ! Column, level and constituent indices
    real(kind_phys) :: fdq          ! mass adjustment factor
    real(kind_phys) :: tot_water_ini ! total water, initial
    real(kind_phys) :: tot_water_now ! total water, present

    integer :: m_cnst
    logical :: water_flag, active_flag
    logical :: is_thermo_water(pcnst)  ! water species that are thermodynamically active

    errmsg = ' '
    errflg = 0

    !-----------------------------------------------------------------------
    ! if dycore is not moist, return as dme_adjust is only for moist mixing ratios

    if (.not. is_dycore_moist) return

    ! The constituent properties are queried on the host (type-bound
    ! procedures), once, not per level.
    do m_cnst = 1, pcnst
      call const_props(m_cnst)%is_water_species(water_flag)
      call const_props(m_cnst)%is_thermo_active(active_flag)
      is_thermo_water(m_cnst) = water_flag .and. active_flag
    end do

    !
    !-----------------------------------------------------------------------
    ! adjust dry mass in each layer back to input value, while conserving
    ! constituents, momentum, and total energy.
    ! Columns are independent; the levels of a column are not (ps and pint
    ! accumulate downward), so the loop runs over columns, levels in order.
    !$acc parallel loop gang vector copyin(is_thermo_water, qini, liqini, iceini) &
    !$acc   copy(ps, pint, pdel, lnpint, rpdel, const_array)
    do i = 1, ncol
      ps(i) = pint(i,1)
      do k = 1, pver
        tot_water_ini = qini(i,k) + liqini(i,k) + iceini(i,k) ! initial total H2O
        tot_water_now = 0.0_kind_phys
        do m = 1, pcnst
          if (is_thermo_water(m)) then
            tot_water_now = tot_water_now + const_array(i,k,m)
          end if
        end do
        fdq = 1._kind_phys + tot_water_now - tot_water_ini
        ! adjust constituents to conserve mass in each layer
        do m = 1, pcnst
          const_array(i,k,m) = const_array(i,k,m) / fdq
        end do
        ! compute new total pressure variables
        pdel  (i,k  ) = pdel(i,k  ) * fdq
        ps(i)         = ps(i)       + pdel(i,k)
        pint  (i,k+1) = pint(i,k  ) + pdel(i,k)
        lnpint(i,k+1) = log(pint(i,k+1))
        rpdel (i,k  ) = 1._kind_phys / pdel(i,k  )
        !note that mid-level variables (e.g. pmid) are not recomputed
      end do
    end do

  end subroutine dme_adjust_run

end module dme_adjust
