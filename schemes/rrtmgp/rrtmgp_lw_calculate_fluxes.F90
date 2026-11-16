module rrtmgp_lw_calculate_fluxes

   use ccpp_kinds, only:  kind_phys
   implicit none
   private

   public :: rrtmgp_lw_calculate_fluxes_run  ! main routine


CONTAINS

   !> \section arg_table_rrtmgp_lw_calculate_fluxes_run  Argument Table
   !! \htmlinclude rrtmgp_lw_calculate_fluxes_run.html
   subroutine rrtmgp_lw_calculate_fluxes_run(num_diag_subcycles, icall, ncol, pverp, nlay, ktopcam, ktoprad, &
      active_calls, dolw, flw, flwc, flns, flnt, flwds, fnl, fcnl, errmsg, errflg)

      use ccpp_fluxes,        only: ty_fluxes_broadband_ccpp
      use ccpp_fluxes_byband, only: ty_fluxes_byband_ccpp
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      integer,                        intent(in) :: num_diag_subcycles  ! Number of diagnostics subcycles
      integer,                        intent(in) :: icall               ! Current diagnostic subcycle
      integer,                        intent(in) :: pverp               ! Number of vertical layer interfaces
      integer,                        intent(in) :: ncol                ! Number of horizontal grid points
      integer,                        intent(in) :: nlay                ! Number of vertical layers in RRTMGP
      integer,                        intent(in) :: ktopcam             ! Index in host model arrays of top level (layer or interface) at which RRTMGP is active
      integer,                        intent(in) :: ktoprad             ! Index in RRTMGP array corresponding to top layer or interface of host model arrays
      logical,                        intent(in) :: active_calls(:)     ! Logical array of flags for whether a specified subcycle is active
      logical,                        intent(in) :: dolw                ! Flag for whether to perform longwave calculation
      type(ty_fluxes_byband_ccpp),    intent(in) :: flw                 ! Longwave all-sky flux object
      type(ty_fluxes_broadband_ccpp), intent(in) :: flwc                ! Longwave clear-sky flux object
      ! Output variables
      real(kind_phys),             intent(inout) :: fnl(:,:)            ! Longwave net radiative flux [W m-2]
      real(kind_phys),             intent(inout) :: fcnl(:,:)           ! Longwave net radiative clear-sky flux [W m-2]
      real(kind_phys),             intent(inout) :: flns(:)             ! Longwave net upward flux at surface [W m-2]
      real(kind_phys),             intent(inout) :: flnt(:)             ! Longwave net outgoing flux at model top [W m-2]
      real(kind_phys),             intent(inout) :: flwds(:)            ! Longwave downward radiative flux at surface [W m-2]

      
      ! CCPP error handling variables
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      integer :: diag_index, idx
      integer :: i, k

      errmsg = ''
      errflg = 0

      ! Persist fluxes from the last radiation timestep in non-radiation timesteps:
      if (.not. dolw) then
         return
      end if

      diag_index = num_diag_subcycles - icall + 1

      ! Don't do anything if this subcycle is inactive
      if (.not. active_calls(diag_index)) then
         return
      end if

      !$acc data copyin(flw%fluxes, flw%fluxes%flux_net, flw%fluxes%flux_dn, &
      !$acc             flwc%fluxes, flwc%fluxes%flux_net) &
      !$acc      copy(fnl, fcnl, flns, flnt, flwds)
      !$acc parallel loop collapse(2)
      do k = 1, size(fnl, 2)
         do i = 1, size(fnl, 1)
            fnl(i, k) = 0.0_kind_phys
            fcnl(i, k) = 0.0_kind_phys
         end do
      end do

      ! RTE-RRTMGP convention for net is (down - up) **CAM assumes (up - down) !!
      !$acc parallel loop collapse(2)
      do k = ktopcam, size(fnl, 2)
         do i = 1, size(fnl, 1)
            fnl(i, k) = -1._kind_phys * flw%fluxes%flux_net(i, k - ktopcam + ktoprad)
            fcnl(i, k) = -1._kind_phys * flwc%fluxes%flux_net(i, k - ktopcam + ktoprad)
         end do
      end do

      !$acc parallel loop
      do i = 1, ncol
         flns(i) = fnl(i, pverp)
         flnt(i) = fnl(i, ktopcam)
         flwds(i) = flw%fluxes%flux_dn(i, nlay+1)
      end do
      !$acc end data

   end subroutine rrtmgp_lw_calculate_fluxes_run

end module rrtmgp_lw_calculate_fluxes
