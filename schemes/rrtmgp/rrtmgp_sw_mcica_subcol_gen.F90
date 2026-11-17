module rrtmgp_sw_mcica_subcol_gen

implicit none
private
save

public :: rrtmgp_sw_mcica_subcol_gen_run

!==================================================================================================
contains
!==================================================================================================

!>
!> \section arg_table_rrtmgp_sw_mcica_subcol_gen_run Argument Table
!! \htmlinclude rrtmgp_sw_mcica_subcol_gen_run.html
subroutine rrtmgp_sw_mcica_subcol_gen_run(dosw, kdist_sw, nswbands, nswgpts, nday, nlay, &
              pver, tiny, idxday, ktopcam, ktoprad, cldfprime, c_cld_tau,       &
              c_cld_tau_w, c_cld_tau_w_g, pmid_day, cloud_sw, errmsg, errflg)
   use ccpp_kinds,              only: kind_phys
   use ccpp_gas_concentrations, only: ty_gas_concs_ccpp
   use ccpp_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp_ccpp
   use ccpp_optical_props,      only: ty_optical_props_2str_ccpp
   use iso_fortran_env,         only: int64
   use atmos_phys_kiss_rng,     only: kiss_rng_seed, kiss_rng_next

   ! Compute combined cloud optical properties.
   ! Create MCICA stochastic arrays for cloud SW optical properties.
   ! Initialize optical properties object (cloud_sw) and load with MCICA columns.

   ! arguments
   class(ty_gas_optics_rrtmgp_ccpp), intent(in)  :: kdist_sw             ! shortwave gas optics object
   integer,                          intent(in)  :: nswbands             ! number of shortwave bands
   integer,                          intent(in)  :: nswgpts              ! number of shortwave g-points
   integer,                          intent(in)  :: nlay                 ! number of layers in radiation calculation (may include "extra layer")
   integer,                          intent(in)  :: nday                 ! number of daylight columns
   integer,                          intent(in)  :: pver                 ! total number of vertical layers
   integer,                          intent(in)  :: ktopcam              ! index in host model arrays of top level (layer or interface) at which RRTMGP is active
   integer,                          intent(in)  :: ktoprad              ! index in RRTMGP array corresponding to top layer or interface of host model arrays
   integer,                          intent(in)  :: idxday(:)            ! indices of daylight columns in the chunk
   real(kind_phys),                  intent(in)  :: tiny                 ! definition of tiny in RRTMGP
   real(kind_phys),                  intent(in)  :: c_cld_tau(:,:,:)     ! combined cloud extinction optical depth
   real(kind_phys),                  intent(in)  :: c_cld_tau_w(:,:,:)   ! combined cloud single scattering albedo * tau
   real(kind_phys),                  intent(in)  :: c_cld_tau_w_g(:,:,:) ! combined cloud asymmetry parameter * w * tau
   real(kind_phys),                  intent(in)  :: cldfprime(:,:)       ! combined cloud fraction
   real(kind_phys),                  intent(in)  :: pmid_day(:,:)        ! air pressure at mid-points [Pa]
   logical,                          intent(in)  :: dosw                 ! Flag to do shortwave radiation this timestep

   type(ty_optical_props_2str_ccpp), intent(out) :: cloud_sw             ! SW cloud optical properties object
   character(len=*),                 intent(out) :: errmsg
   integer,                          intent(out) :: errflg

   ! Local variables

   integer :: i, k, n, icol
   integer :: igpt, isubcol
   integer, parameter :: changeseed = 1
   real(kind_phys), parameter :: cldmin = 1.0e-80_kind_phys  ! min cloud fraction

   integer :: gpt2band(nswgpts)          ! band of each g-point
   integer :: kiss_seed1, kiss_seed2, kiss_seed3, kiss_seed4
   integer(int64) :: s1, s2, s3, s4      ! KISS state of a column
   real(kind_phys) :: rand_num           ! random number (KISS)
   real(kind_phys) :: tau_k, tau_w_k, tau_w_g_k, tauc_k, asmc_k

   ! The daylight columns (1:nday) and the CAM layers in the radiation
   ! calculation (1:pver-ktopcam+1, not the "extra layer")
   real(kind_phys) :: cdf(nswgpts,nday,pver-ktopcam+1)    ! random numbers
   real(kind_phys) :: cldfrac(nday,pver-ktopcam+1)        ! cloud fraction clipped to cldmin
   real(kind_phys) :: tauc(nswbands,nday,pver-ktopcam+1)  ! cloud optical depth
   real(kind_phys) :: ssac(nswbands,nday,pver-ktopcam+1)  ! single scattering albedo
   real(kind_phys) :: asmc(nswbands,nday,pver-ktopcam+1)  ! asymmetry parameter
   integer :: nver
   !--------------------------------------------------------------------------------

   errmsg = ''
   errflg = 0

   ! if no daylight columns the cloud_sw object isn't initialized
   if (nday > 0 .and. dosw) then

      ! Initialize object for SW cloud optical properties.
      errmsg = cloud_sw%optical_props%alloc_2str(nday, nlay, kdist_sw%gas_props)
      if (len_trim(errmsg) > 0) then
         errflg = 1
         return
      end if

      ! number of CAM's layers in radiation calculation.  Does not include the "extra layer".
      nver = pver - ktopcam + 1

      ! Band of each g-point (the gas optics' convert_gpt2band)
      gpt2band(:) = kdist_sw%gas_props%get_gpoint_bands()

      ! The object is only copied in (copying it back could overwrite scalars
      ! set on the host); its tau, ssa and g arrays are written out.
      !$acc data copyin(idxday, cldfprime, c_cld_tau, c_cld_tau_w, c_cld_tau_w_g, pmid_day,   &
      !$acc             gpt2band, cloud_sw%optical_props)                                   &
      !$acc      copyout(cloud_sw%optical_props%tau, cloud_sw%optical_props%ssa,            &
      !$acc              cloud_sw%optical_props%g)                                          &
      !$acc      create(cdf, cldfrac, tauc, ssac, asmc)

      ! Subset data so just the daylight columns, and the number of CAM layers in the
      ! radiation calculation are used by MCICA to produce subcolumns, and compute
      ! the optical properties needed for the 2-stream calculations.  These
      ! calculations are the same as the RRTMG version.
      !$acc parallel loop gang vector collapse(2) default(present) &
      !$acc   private(icol, tau_k, tau_w_k, tau_w_g_k, tauc_k, asmc_k)
      do k = 1, nver
         do i = 1, nday
            icol = idxday(i)
            ! clip cloud fraction
            cldfrac(i,k) = cldfprime(icol,ktopcam+k-1)
            if (cldfrac(i,k) < cldmin) then
               cldfrac(i,k) = 0._kind_phys
            end if
            !$acc loop seq
            do n = 1, nswbands
               tau_k     = c_cld_tau    (n,icol,ktopcam+k-1)
               tau_w_k   = c_cld_tau_w  (n,icol,ktopcam+k-1)
               tau_w_g_k = c_cld_tau_w_g(n,icol,ktopcam+k-1)
               ! set cloud optical depth, clip @ zero
               if (tau_k > 0.0_kind_phys) then
                  tauc_k = tau_k
               else
                  tauc_k = 0.0_kind_phys
               end if
               ! set value of asymmetry
               if (tau_w_k > 0.0_kind_phys) then
                  asmc_k = tau_w_g_k / max(tau_w_k, tiny)
               else
                  asmc_k = 0.0_kind_phys
               end if
               ! set value of single scattering albedo; asymmetry to zero when tauc = 0
               if (tauc_k > 0.0_kind_phys) then
                  ssac(n,i,k) = max(tau_w_k, tiny) / max(tauc_k, tiny)
               else
                  ssac(n,i,k) = 1.0_kind_phys
                  asmc_k = 0.0_kind_phys
               end if
               tauc(n,i,k) = tauc_k
               asmc(n,i,k) = asmc_k
            end do
         end do
      end do

      ! Random numbers: one KISS stream per daylight column, seeded from the
      ! state of the column (pmid_day of the bottom four layers), advanced by
      ! changeseed values, then one number per subcolumn and level (level fastest).
      !$acc parallel loop gang vector default(present) &
      !$acc   private(kiss_seed1, kiss_seed2, kiss_seed3, kiss_seed4, s1, s2, s3, s4, rand_num)
      do i = 1, nday
         kiss_seed1 = (pmid_day(i,nlay)   - int(pmid_day(i,nlay)))    * 1000000000
         kiss_seed2 = (pmid_day(i,nlay-1) - int(pmid_day(i,nlay-1)))  * 1000000000
         kiss_seed3 = (pmid_day(i,nlay-2) - int(pmid_day(i,nlay-2)))  * 1000000000
         kiss_seed4 = (pmid_day(i,nlay-3) - int(pmid_day(i,nlay-3)))  * 1000000000
         s1 = kiss_rng_seed(kiss_seed1)
         s2 = kiss_rng_seed(kiss_seed2)
         s3 = kiss_rng_seed(kiss_seed3)
         s4 = kiss_rng_seed(kiss_seed4)
         !$acc loop seq
         do n = 1, changeseed
            rand_num = kiss_rng_next(s1, s2, s3, s4)
         end do
         !$acc loop seq
         do isubcol = 1, nswgpts
            !$acc loop seq
            do k = 1, nver
               cdf(isubcol,i,k) = kiss_rng_next(s1, s2, s3, s4)
            end do
         end do
      end do

      ! Maximum-Random overlap
      ! i) pick a random number for top layer.
      ! ii) walk down the column:
      !    - if the layer above is cloudy, use the same random number as in the layer above
      !    - if the layer above is clear, use a new random number
      ! (each subcolumn of each column is independent; the levels are sequential)
      !$acc parallel loop gang vector collapse(2) default(present)
      do i = 1, nday
         do isubcol = 1, nswgpts
            !$acc loop seq
            do k = 2, nver
               if (cdf(isubcol,i,k-1) > 1._kind_phys - cldfrac(i,k-1) ) then
                  cdf(isubcol,i,k) = cdf(isubcol,i,k-1)
               else
                  cdf(isubcol,i,k) = cdf(isubcol,i,k) * (1._kind_phys - cldfrac(i,k-1))
               end if
            end do
         end do
      end do

      ! If there is an extra layer in the radiation then this initialization
      ! will provide the optical properties there.
      !$acc parallel loop gang vector collapse(3) default(present)
      do igpt = 1, size(cloud_sw%optical_props%tau, 3)
         do k = 1, size(cloud_sw%optical_props%tau, 2)
            do i = 1, size(cloud_sw%optical_props%tau, 1)
               cloud_sw%optical_props%tau(i,k,igpt) = 0.0_kind_phys
               cloud_sw%optical_props%ssa(i,k,igpt) = 1.0_kind_phys
               cloud_sw%optical_props%g  (i,k,igpt) = 0.0_kind_phys
            end do
         end do
      end do

      ! -- generate subcolumns for homogeneous clouds -----
      ! where there is a cloud (the subcolumn is cloudy), set the subcolumn
      ! cloud properties on the g-points; incoming tauc should be in-cloud
      ! quantites and not grid-averaged quantities
      !$acc parallel loop gang vector collapse(3) default(present) private(n)
      do igpt = 1, nswgpts
         do k = 1, nver
            do i = 1, nday
               if (cdf(igpt,i,k) >= 1._kind_phys - cldfrac(i,k) .and. &
                   cldfrac(i,k) > 0._kind_phys) then
                  n = gpt2band(igpt)
                  cloud_sw%optical_props%g  (i,ktoprad+k-1,igpt) = asmc(n,i,k)
                  cloud_sw%optical_props%ssa(i,ktoprad+k-1,igpt) = ssac(n,i,k)
                  cloud_sw%optical_props%tau(i,ktoprad+k-1,igpt) = tauc(n,i,k)
               else
                  cloud_sw%optical_props%g  (i,ktoprad+k-1,igpt) = 0._kind_phys
                  cloud_sw%optical_props%ssa(i,ktoprad+k-1,igpt) = 1._kind_phys
                  cloud_sw%optical_props%tau(i,ktoprad+k-1,igpt) = 0._kind_phys
               end if
            end do
         end do
      end do

      ! validate checks that: tau > 0, ssa is in range [0,1], and g is in range [-1,1];
      ! delta scaling adjusts for forward scattering (both on the device data)
      errmsg = cloud_sw%optical_props%validate()
      if (len_trim(errmsg) == 0) then
         errmsg = cloud_sw%optical_props%delta_scale()
      end if

      !$acc end data

      if (len_trim(errmsg) > 0) then
         errflg = 1
         return
      end if

   end if

end subroutine rrtmgp_sw_mcica_subcol_gen_run

end module rrtmgp_sw_mcica_subcol_gen
