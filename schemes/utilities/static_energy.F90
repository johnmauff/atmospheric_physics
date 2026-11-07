! Updates dry static energy from updated temperature and geopotential height.
! If ran, must be after geopotential_temp scheme.
module static_energy

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: update_dry_static_energy_run

CONTAINS

   !> \section arg_table_update_dry_static_energy_run  Argument Table
   !! \htmlinclude update_dry_static_energy_run.html
   subroutine update_dry_static_energy_run(nz, gravit, temp, zm, phis,        &
        st_energy, cpair, errcode, errmsg)

      ! Dummy arguments
      integer,            intent(in)  :: nz             ! Num vertical layers
      real(kind_phys),    intent(in)  :: gravit         ! gravitational acceleration
      real(kind_phys),    intent(in)  :: temp(:,:)      ! air temperature
      real(kind_phys),    intent(in)  :: zm(:,:)        ! geopotential height
      real(kind_phys),    intent(in)  :: phis(:)        ! surface geopotential
      real(kind_phys),    intent(out) :: st_energy(:,:) ! dry static energy
      real(kind_phys),    intent(in)  :: cpair(:,:)     ! specific heat, dry air
      integer,            intent(out) :: errcode
      character(len=512), intent(out) :: errmsg

      ! Local variables
      integer :: icol, klev

      errcode = 0
      errmsg = ''

      ! st_energy is (horizontal_loop_extent, vertical_layer_dimension).
      !$acc parallel loop collapse(2) copyin(temp, zm, phis, cpair) copy(st_energy)
      do klev = 1, nz
         do icol = 1, size(st_energy, 1)
            st_energy(icol, klev) = (temp(icol, klev) * cpair(icol, klev)) +   &
                 (gravit * zm(icol, klev)) + phis(icol)
         end do
      end do

   end subroutine update_dry_static_energy_run

end module static_energy
