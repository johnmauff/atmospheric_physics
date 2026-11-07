module rrtmgp_sw_calculate_heating_rate
   public :: rrtmgp_sw_calculate_heating_rate_run

CONTAINS
   !> \section arg_table_rrtmgp_sw_calculate_heating_rate_run  Argument Table
   !! \htmlinclude rrtmgp_sw_calculate_heating_rate_run.html
   subroutine rrtmgp_sw_calculate_heating_rate_run(ktopcam, pver, gravit, rpdel, flux_net, flux_net_clrsky, &
                   hrate, hrate_clrsky, errmsg, errflg)
      use ccpp_kinds, only: kind_phys
      integer,          intent(in) :: ktopcam               ! Vertical index at top level where RRTMGP is active
      integer,          intent(in) :: pver                  ! Number of vertical layers
      real(kind_phys),  intent(in) :: gravit                ! Standard gravitiational acceleration [m s-2]
      real(kind_phys),  intent(in) :: rpdel(:,:)            ! Reciprocal of air pressure thickness [Pa-1]
      real(kind_phys),  intent(in) :: flux_net(:,:)         ! Shortwave net radiative flux [W m-2]
      real(kind_phys),  intent(in) :: flux_net_clrsky(:,:)  ! Shortwave net radiative clear-sky flux [W m-2]
      real(kind_phys), intent(out) :: hrate(:,:)            ! Tendency of dry air enthalpy due to SW radiation [J kg-1 s-1]
      real(kind_phys), intent(out) :: hrate_clrsky(:,:)     ! Tendency of dry air enthalpy due to SW clear-sky radiation [J kg-1 s-1]
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      integer :: i, k

      ! Set error variables
      errmsg = ''
      errflg = 0

      !$acc data copyin(rpdel, flux_net, flux_net_clrsky) copy(hrate, hrate_clrsky)
      !$acc kernels
      hrate = 0.0_kind_phys
      !$acc end kernels

      !$acc parallel loop collapse(2)
      do k = ktopcam, pver
         do i = 1, size(hrate, 1)
            ! top - bottom
            hrate(i,k) = (flux_net(i,k) - flux_net(i,k+1)) * &
                    gravit * rpdel(i,k)
            hrate_clrsky(i,k) = (flux_net_clrsky(i,k) - flux_net_clrsky(i,k+1)) * &
                    gravit * rpdel(i,k)
         end do
      end do
      !$acc end data

   end subroutine rrtmgp_sw_calculate_heating_rate_run

end module rrtmgp_sw_calculate_heating_rate
