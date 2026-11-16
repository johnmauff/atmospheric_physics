module rrtmgp_sw_calculate_fluxes

   use ccpp_kinds, only:  kind_phys
   implicit none
   private

   public :: rrtmgp_sw_calculate_fluxes_run  ! main routine


CONTAINS

   !> \section arg_table_rrtmgp_sw_calculate_fluxes_run  Argument Table
   !! \htmlinclude rrtmgp_sw_calculate_fluxes_run.html
   subroutine rrtmgp_sw_calculate_fluxes_run(num_diag_subcycles, icall, ncol, pverp, nlay, nday, idxday, ktopcam, ktoprad, &
      active_calls, dosw, fsw, fswc, fns, fcns, fsns, fsnt, soll, sols, solld, solsd, errmsg, errflg)

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
      integer,                        intent(in) :: nday                ! Daytime points dimension
      integer,                        intent(in) :: ktopcam             ! Index in host model arrays of top level (layer or interface) at which RRTMGP is active
      integer,                        intent(in) :: ktoprad             ! Index in RRTMGP array corresponding to top layer or interface of host model arrays
      integer,                        intent(in) :: idxday(:)           ! Daytime points indices
      logical,                        intent(in) :: active_calls(:)     ! Logical array of flags for whether a specified subcycle is active
      logical,                        intent(in) :: dosw                ! Flag for whether to perform shortwave calculation
      type(ty_fluxes_byband_ccpp),    intent(in) :: fsw                 ! Shortwave all-sky flux object
      type(ty_fluxes_broadband_ccpp), intent(in) :: fswc                ! Shortwave clear-sky flux object
      ! Output variables
      real(kind_phys),             intent(inout) :: fns(:,:)            ! Shortwave net radiative flux [W m-2]
      real(kind_phys),             intent(inout) :: fcns(:,:)           ! Shortwave net radiative clear-sky flux [W m-2]
      real(kind_phys),             intent(inout) :: fsns(:)             ! Shortwave net upward flux at surface [W m-2]
      real(kind_phys),             intent(inout) :: fsnt(:)             ! Shortwave net outgoing flux at model top [W m-2]
      real(kind_phys),             intent(inout) :: soll(:)             ! Direct solar radiative flux at surface >= 700nm [W m-2]
      real(kind_phys),             intent(inout) :: sols(:)             ! Direct solar radiative flux at surface < 700nm [W m-2]
      real(kind_phys),             intent(inout) :: solld(:)            ! Diffuse solar radiative flux at surface >= 700nm [W m-2]
      real(kind_phys),             intent(inout) :: solsd(:)            ! Diffuse solar radiative flux at surface < 700nm [W m-2]

      
      ! CCPP error handling variables
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      integer :: diag_index, idx
      integer :: i, k, ib
      real(kind_phys) :: s, s2   ! band sums (direct, diffuse)

      errmsg = ''
      errflg = 0

      ! Persist fluxes from the last radiation timestep in non-radiation timesteps:
      if (.not. dosw) then
         return
      end if

      diag_index = num_diag_subcycles - icall + 1

      ! Don't do anything if this subcycle is inactive
      if (.not. active_calls(diag_index)) then
         return
      end if

      ! Initialize to provide 0.0 values for night columns.
      !$acc data copyin(idxday, fsw%fluxes, fsw%fluxes%flux_net, fsw%fluxes%bnd_flux_dn, &
      !$acc             fsw%fluxes%bnd_flux_dn_dir, fswc%fluxes, fswc%fluxes%flux_net) &
      !$acc      copy(fns, fcns, fsns, fsnt, soll, sols, solld, solsd)
      !$acc parallel loop collapse(2)
      do k = 1, size(fns, 2)
         do i = 1, size(fns, 1)
            fns(i, k) = 0.0_kind_phys
            fcns(i, k) = 0.0_kind_phys
         end do
      end do
      !$acc parallel loop
      do i = 1, size(fsns)
         fsns(i) = 0.0_kind_phys
         fsnt(i) = 0.0_kind_phys
         soll(i) = 0.0_kind_phys
         sols(i) = 0.0_kind_phys
         solld(i) = 0.0_kind_phys
         solsd(i) = 0.0_kind_phys
      end do

      ! Daytime columns: RRTMGP levels -> host model levels
      !$acc parallel loop collapse(2)
      do k = ktopcam, size(fns, 2)
         do idx = 1, nday
            fns(idxday(idx), k) = fsw%fluxes%flux_net(idx, k - ktopcam + ktoprad)
            fcns(idxday(idx), k) = fswc%fluxes%flux_net(idx, k - ktopcam + ktoprad)
         end do
      end do

      !$acc parallel loop
      do i = 1, size(fsns)
         fsns(i) = fns(i, pverp)
         fsnt(i) = fns(i, ktopcam)
      end do

      ! Export surface fluxes
      ! sols(pcols)      Direct solar rad on surface (< 0.7)
      ! soll(pcols)      Direct solar rad on surface (>= 0.7)
      ! RRTMGP: Near-IR bands (1-10), 820-16000 cm-1, 0.625-12.195 microns
      ! Put half of band 10 in each of the UV/visible and near-IR values,
      ! since this band straddles 0.7 microns:
      ! UV/visible bands 10-13, 16000-50000 cm-1, 0.200-0.625 micron
      ! The diffuse flux is the total minus the direct (bnd_flux_dn - bnd_flux_dn_dir).
      ! The band sums add left to right from zero, as SUM does on the host
      ! (loop seq: not vectorized into a reduction of another order).

      ! The following calculations are hard-coded and will need to be modified
      ! if the RRTMGP flux bands type/array ever changes
      !$acc parallel loop private(s, s2)
      do idx = 1, nday
         s = 0.0_kind_phys
         s2 = 0.0_kind_phys
         !$acc loop seq
         do ib = 1, 9
            s = s + fsw%fluxes%bnd_flux_dn_dir(idx, nlay+1, ib)
            s2 = s2 + (fsw%fluxes%bnd_flux_dn(idx, nlay+1, ib) - fsw%fluxes%bnd_flux_dn_dir(idx, nlay+1, ib))
         end do
         soll(idxday(idx)) = s + 0.5_kind_phys * fsw%fluxes%bnd_flux_dn_dir(idx, nlay+1, 10)
         solld(idxday(idx)) = s2 + 0.5_kind_phys * &
            (fsw%fluxes%bnd_flux_dn(idx, nlay+1, 10) - fsw%fluxes%bnd_flux_dn_dir(idx, nlay+1, 10))

         s = 0.0_kind_phys
         s2 = 0.0_kind_phys
         !$acc loop seq
         do ib = 11, 14
            s = s + fsw%fluxes%bnd_flux_dn_dir(idx, nlay+1, ib)
            s2 = s2 + (fsw%fluxes%bnd_flux_dn(idx, nlay+1, ib) - fsw%fluxes%bnd_flux_dn_dir(idx, nlay+1, ib))
         end do
         sols(idxday(idx)) = 0.5_kind_phys * fsw%fluxes%bnd_flux_dn_dir(idx, nlay+1, 10) + s
         solsd(idxday(idx)) = 0.5_kind_phys * &
            (fsw%fluxes%bnd_flux_dn(idx, nlay+1, 10) - fsw%fluxes%bnd_flux_dn_dir(idx, nlay+1, 10)) + s2
      end do
      !$acc end data

   end subroutine rrtmgp_sw_calculate_fluxes_run

end module rrtmgp_sw_calculate_fluxes
