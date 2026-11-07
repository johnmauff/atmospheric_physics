module physics_tendency_updaters

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: apply_tendency_of_eastward_wind_run
  public :: apply_tendency_of_northward_wind_run
  public :: apply_heating_rate_run
  public :: apply_tendency_of_air_temperature_run
  public :: apply_constituent_tendencies_run

CONTAINS

   !> \section arg_table_apply_tendency_of_eastward_wind_run  Argument Table
   !! \htmlinclude apply_tendency_of_eastward_wind_run.html
   subroutine apply_tendency_of_eastward_wind_run(nz, dudt, u, dudt_total, dt,             &
        errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz              ! Num vertical  layers
      real(kind_phys),    intent(inout) :: dudt(:,:)       ! tendency of eastward wind
      real(kind_phys),    intent(inout) :: u(:,:)          ! eastward wind
      real(kind_phys),    intent(inout) :: dudt_total(:,:) ! total tendency of eastward wind
      real(kind_phys),    intent(in)    :: dt              ! physics time step
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variables
      integer :: icol, klev, ncol

      errcode = 0
      errmsg = ''

      ! All three arrays are (horizontal_loop_extent, vertical_layer_dimension),
      ! so one kernel applies the tendency and then zeroes it, point by point.
      ncol = size(u, 1)
      !$acc parallel loop collapse(2) copy(dudt, u, dudt_total)
      do klev = 1, nz
         do icol = 1, ncol
            u(icol, klev) = u(icol, klev) + (dudt(icol, klev) * dt)
            dudt_total(icol, klev) = dudt_total(icol, klev) + dudt(icol, klev)
            dudt(icol, klev) = 0.0_kind_phys
         end do
      end do

   end subroutine apply_tendency_of_eastward_wind_run

   !> \section arg_table_apply_tendency_of_northward_wind_run  Argument Table
   !! \htmlinclude apply_tendency_of_northward_wind_run.html
   subroutine apply_tendency_of_northward_wind_run(nz, dvdt, v, dvdt_total, dt,             &
        errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz              ! Num vertical layers
      real(kind_phys),    intent(inout) :: dvdt(:,:)       ! tendency of northward wind
      real(kind_phys),    intent(inout) :: v(:,:)          ! northward wind
      real(kind_phys),    intent(inout) :: dvdt_total(:,:) ! total tendency of northward wind
      real(kind_phys),    intent(in)    :: dt              ! physics time step
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variables
      integer :: icol, klev, ncol

      errcode = 0
      errmsg = ''

      ! All three arrays are (horizontal_loop_extent, vertical_layer_dimension),
      ! so one kernel applies the tendency and then zeroes it, point by point.
      ncol = size(v, 1)
      !$acc parallel loop collapse(2) copy(dvdt, v, dvdt_total)
      do klev = 1, nz
         do icol = 1, ncol
            v(icol, klev) = v(icol, klev) + (dvdt(icol, klev) * dt)
            dvdt_total(icol, klev) = dvdt_total(icol, klev) + dvdt(icol, klev)
            dvdt(icol, klev) = 0.0_kind_phys
         end do
      end do

   end subroutine apply_tendency_of_northward_wind_run

   !> \section arg_table_apply_heating_rate_run  Argument Table
   !! \htmlinclude apply_heating_rate_run.html
   subroutine apply_heating_rate_run(nz, heating_rate, temp, dTdt_total, dt, cpair,  &
        errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz                ! Num vertical  layers
      real(kind_phys),    intent(inout) :: heating_rate(:,:) ! heating rate
      real(kind_phys),    intent(inout) :: temp(:,:)         ! air temperature
      real(kind_phys),    intent(inout) :: dTdt_total(:,:)   ! total temperature tend.
      real(kind_phys),    intent(in)    :: dt                ! physics time step
      real(kind_phys),    intent(in)    :: cpair(:,:)        ! specific heat, dry air
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variables
      integer :: icol, klev, ncol

      errcode = 0
      errmsg = ''

      ! All four arrays are (horizontal_loop_extent, vertical_layer_dimension),
      ! so one kernel applies the heating rate and then zeroes it, point by point.
      ncol = size(temp, 1)
      !$acc parallel loop collapse(2) copyin(cpair) copy(heating_rate, temp, dTdt_total)
      do klev = 1, nz
         do icol = 1, ncol
            temp(icol, klev) = temp(icol, klev) + (heating_rate(icol, klev) * dt / cpair(icol, klev))
            dTdt_total(icol, klev) = dTdt_total(icol, klev) + (heating_rate(icol, klev) / cpair(icol, klev))
            heating_rate(icol, klev) = 0.0_kind_phys
         end do
      end do

   end subroutine apply_heating_rate_run

   !> \section arg_table_apply_tendency_of_air_temperature_run  Argument Table
   !! \htmlinclude apply_tendency_of_air_temperature_run.html
   subroutine apply_tendency_of_air_temperature_run(nz, t_tend, temp, dtdT_total,   &
        dt, errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz              ! Num vertical  layers
      real(kind_phys),    intent(inout) :: t_tend(:,:)     ! temperature tendency
      real(kind_phys),    intent(inout) :: temp(:,:)       ! air temperature
      real(kind_phys),    intent(inout) :: dTdt_total(:,:) ! total temp. tendency
      real(kind_phys),    intent(in)    :: dt              ! physics time step
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variables
      integer :: icol, klev, ncol

      errcode = 0
      errmsg = ''

      ! All three arrays are (horizontal_loop_extent, vertical_layer_dimension),
      ! so one kernel applies the tendency and then zeroes it, point by point.
      ncol = size(temp, 1)
      !$acc parallel loop collapse(2) copy(t_tend, temp, dTdt_total)
      do klev = 1, nz
         do icol = 1, ncol
            temp(icol, klev) = temp(icol, klev) + (t_tend(icol, klev) * dt)
            dTdt_total(icol, klev) = dTdt_total(icol, klev) + t_tend(icol, klev)
            t_tend(icol, klev) = 0.0_kind_phys
         end do
      end do

   end subroutine apply_tendency_of_air_temperature_run

   !> \section arg_table_apply_constituent_tendencies_run Argument Table
   !!! \htmlinclude apply_constituent_tendencies_run.html
   subroutine apply_constituent_tendencies_run(nz, const_tend, const, dt, errcode, errmsg)
      ! Dummy arguments
      integer,            intent(in)    :: nz                 ! Num vertical layers
      real(kind_phys),    intent(inout) :: const_tend(:,:,:)  ! constituent tendency array
      real(kind_phys),    intent(inout) :: const(:,:,:)       ! constituent state array
      real(kind_phys),    intent(in)    :: dt                 ! physics time step
      integer,            intent(out)   :: errcode
      character(len=512), intent(out)   :: errmsg

      ! Local variables
      integer :: icol, klev, m

      errcode = 0
      errmsg = ''

      ! Both arrays are (horizontal_loop_extent, vertical_layer_dimension,
      ! number_of_ccpp_constituents), so one kernel applies the tendency and
      ! then zeroes it, point by point.
      !$acc parallel loop collapse(3) copy(const_tend, const)
      do m = 1, size(const, 3)
         do klev = 1, nz
            do icol = 1, size(const, 1)
               const(icol, klev, m) = const(icol, klev, m) + (const_tend(icol, klev, m) * dt)
               const_tend(icol, klev, m) = 0._kind_phys
            end do
         end do
      end do

   end subroutine apply_constituent_tendencies_run

end module physics_tendency_updaters
