!> \file rrtmgp_sw_aerosols.F90
!!

!> This module sets the RRTMGP aerosol shortwave optical properties
module rrtmgp_sw_aerosols

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public rrtmgp_sw_aerosols_run

  ! Mapping from RRTMG shortwave bands to RRTMGP
  integer, parameter, dimension(14) :: rrtmg_to_rrtmgp_swbands = &
     [ 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 ]

  real(kind_phys), parameter :: tiny = 1.0e-80_kind_phys

contains

!> \section arg_table_rrtmgp_sw_aerosols_run Argument Table
!! \htmlinclude rrtmgp_sw_aerosols_run.html
!!
   subroutine rrtmgp_sw_aerosols_run( &
     doswrad, ncol, pver, nswbands, nday, idxday, ktopcam, ktoprad, &
     aer_tau, aer_tau_w, aer_tau_w_g, &
     aer_sw, errmsg, errflg)

    use ccpp_optical_props, only: ty_optical_props_2str_ccpp

    ! Inputs
    logical,          intent(in)    :: doswrad                    !< do shortwave calculation
    integer,          intent(in)    :: ncol                       !< number of columns
    integer,          intent(in)    :: pver                       !< number of vertical layers
    integer,          intent(in)    :: nswbands                   !< number of SW bands
    integer,          intent(in)    :: nday                       !< daytime columns count
    integer,          intent(in)    :: idxday(:)                  !< daytime column indices
    integer,          intent(in)    :: ktopcam                    !< top CAM level for RRTMGP
    integer,          intent(in)    :: ktoprad                    !< top RRTMGP level

    real(kind_phys),  intent(in)    :: aer_tau(:, :, :)           !< SW extinction OD (ncol,pver,nswbands)
    real(kind_phys),  intent(in)    :: aer_tau_w(:, :, :)         !< SW ssa*tau (ncol,pver,nswbands)
    real(kind_phys),  intent(in)    :: aer_tau_w_g(:, :, :)       !< SW asy*ssa*tau (ncol,pver,nswbands)

    ! Outputs
    class(ty_optical_props_2str_ccpp), intent(inout) :: aer_sw    !< Aerosol optical properties object

    character(len=*), intent(out) :: errmsg                       !< CCPP error message
    integer,          intent(out) :: errflg                       !< CCPP error flag

    ! Local variables
    integer :: i, k, ib, ic, kc, bc
    integer :: band_map(size(rrtmg_to_rrtmgp_swbands))
    real(kind_phys) :: t, tw

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doswrad .or. nday == 0) return

    ! Band map (RRTMG band order -> RRTMGP), as data the kernels read
    band_map = rrtmg_to_rrtmgp_swbands

    !$acc data copyin(aer_tau, aer_tau_w, aer_tau_w_g, idxday, band_map, aer_sw%optical_props) &
    !$acc      copyout(aer_sw%optical_props%tau, aer_sw%optical_props%ssa, aer_sw%optical_props%g)

    ! Initialize defaults: tau=0, ssa=1, g=0
    !$acc parallel loop collapse(3)
    do ib = 1, size(aer_sw%optical_props%tau, 3)
      do k = 1, size(aer_sw%optical_props%tau, 2)
        do i = 1, size(aer_sw%optical_props%tau, 1)
          aer_sw%optical_props%tau(i, k, ib) = 0.0_kind_phys
          aer_sw%optical_props%ssa(i, k, ib) = 1.0_kind_phys
          aer_sw%optical_props%g(i, k, ib)   = 0.0_kind_phys
        end do
      end do
    end do

    ! Subset to daytime columns, reorder from RRTMG band order to RRTMGP band
    ! order and map CAM vertical to RRTMGP vertical
    !$acc parallel loop collapse(3) private(ic, kc, bc, t, tw)
    do ib = 1, nswbands
      do k = ktoprad, size(aer_sw%optical_props%tau, 2)
        do i = 1, nday
          ic = idxday(i)
          kc = k - ktoprad + ktopcam
          bc = band_map(ib)
          t  = aer_tau(ic, kc, bc)
          tw = aer_tau_w(ic, kc, bc)

          ! Set aerosol optical depth, clip to zero
          aer_sw%optical_props%tau(i, k, ib) = max(t, 0._kind_phys)

          ! Set single scattering albedo: ssa = tau_w / tau
          if (t > 0._kind_phys) then
            aer_sw%optical_props%ssa(i, k, ib) = tw / t
          else
            aer_sw%optical_props%ssa(i, k, ib) = 1._kind_phys
          end if

          ! Set asymmetry parameter: g = tau_w_g / tau_w
          if (tw > tiny) then
            aer_sw%optical_props%g(i, k, ib) = aer_tau_w_g(ic, kc, bc) / tw
          else
            aer_sw%optical_props%g(i, k, ib) = 0._kind_phys
          end if
        end do
      end do
    end do

    ! Impose limits on the components
    !$acc parallel loop collapse(3)
    do ib = 1, size(aer_sw%optical_props%tau, 3)
      do k = 1, size(aer_sw%optical_props%tau, 2)
        do i = 1, size(aer_sw%optical_props%tau, 1)
          aer_sw%optical_props%ssa(i, k, ib) = &
            min(max(aer_sw%optical_props%ssa(i, k, ib), 0._kind_phys), 1._kind_phys)
          aer_sw%optical_props%g(i, k, ib) = &
            min(max(aer_sw%optical_props%g(i, k, ib), -1._kind_phys), 1._kind_phys)
        end do
      end do
    end do
    !$acc end data

  end subroutine rrtmgp_sw_aerosols_run
end module rrtmgp_sw_aerosols
