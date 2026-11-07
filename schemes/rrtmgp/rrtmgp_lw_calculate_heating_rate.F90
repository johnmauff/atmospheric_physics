module rrtmgp_lw_calculate_heating_rate
   public :: rrtmgp_lw_calculate_heating_rate_run

CONTAINS
   !> \section arg_table_rrtmgp_lw_calculate_heating_rate_run  Argument Table
   !! \htmlinclude rrtmgp_lw_calculate_heating_rate_run.html
   subroutine rrtmgp_lw_calculate_heating_rate_run(ktopcam, pver, gravit, rpdel, flux_net, &
                   flux_net_clrsky, hrate, hrate_clrsky, errmsg, errflg)
      use ccpp_kinds, only: kind_phys
      integer,          intent(in) :: ktopcam              ! Index in host model arrays of top level (layer or interface) at which RRTMGP is active
      integer,          intent(in) :: pver                 ! Number of vertical levels
      real(kind_phys),  intent(in) :: gravit               ! Standard gravitational acceleration [m s-2]
      real(kind_phys),  intent(in) :: rpdel(:,:)           ! Reciprocal of air pressure thickness [Pa-1] 
      real(kind_phys),  intent(in) :: flux_net(:,:)        ! Longwave net radiative flux [W m-2]
      real(kind_phys),  intent(in) :: flux_net_clrsky(:,:) ! Longwave net radiative clear-sky flux [W m-2]
      real(kind_phys), intent(out) :: hrate(:,:)           ! Tendency of dry air enthalpy due to LW radiation [J kg-1 s-1]
      real(kind_phys), intent(out) :: hrate_clrsky(:,:)    ! Tendency of dry air enthalpy due to clear-sky LW radiation [J kg-1 s-1]
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      integer :: i, k

      ! Set error variables
      errmsg = ''
      errflg = 0

      !$acc data copyin(rpdel, flux_net, flux_net_clrsky) copy(hrate, hrate_clrsky)
      !$acc kernels
      hrate = 0.0_kind_phys
      hrate_clrsky = 0.0_kind_phys
      !$acc end kernels

      !$acc parallel loop collapse(2)
      do k = ktopcam, pver
         do i = 1, size(hrate, 1)
            ! (flux divergence as bottom-MINUS-top) * g/dp
            hrate(i,k) = (flux_net(i,k+1) - flux_net(i,k)) * &
                                 gravit * rpdel(i,k)
            hrate_clrsky(i,k) = (flux_net_clrsky(i,k+1) - flux_net_clrsky(i,k)) * &
                                 gravit * rpdel(i,k)
         end do
      end do
      !$acc end data

   end subroutine rrtmgp_lw_calculate_heating_rate_run

end module rrtmgp_lw_calculate_heating_rate
