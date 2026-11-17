module rrtmgp_lw_mcica_subcol_gen

!----------------------------------------------------------------------------------------
! 
! Purpose: Create McICA stochastic arrays for lw cloud optical properties.
! Input cloud optical properties directly: cloud optical depth, single
! scattering albedo and asymmetry parameter.  Output will be stochastic
! arrays of these variables.  (longwave scattering is not yet available)
!
! Original code: From RRTMG, with the following copyright notice,
! based on Raisanen et al., QJRMS, 2004:
!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
! This code is a refactored version of code originally in the files
! rrtmgp_lw_mcica_subcol_gen.F90 and mcica_subcol_gen_sw.F90
! 
! Uses the KISS random number generator.
!
! Overlap assumption: maximum-random.
! 
!----------------------------------------------------------------------------------------

implicit none
private
save

public :: rrtmgp_lw_mcica_subcol_gen_run

!========================================================================================
contains
!========================================================================================

!>
!> \section arg_table_rrtmgp_lw_mcica_subcol_gen_run Argument Table
!! \htmlinclude rrtmgp_lw_mcica_subcol_gen_run.html
subroutine rrtmgp_lw_mcica_subcol_gen_run( &
   dolw, ktoprad, ktopcam, kdist, nbnd, ngpt, ncol, pver, nlay, &
   cldfprime, c_cld_lw_abs, changeseed, pmid, cloud_lw,     &
   errmsg, errflg )
   use ccpp_kinds,             only: kind_phys
   use iso_fortran_env,        only: int64
   use atmos_phys_kiss_rng,    only: kiss_rng_seed, kiss_rng_next
   use ccpp_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp_ccpp
   use ccpp_optical_props,     only: ty_optical_props_1scl_ccpp

   ! Arrays use CAM vertical index convention: index increases from top to bottom.
   ! This index ordering is assumed in the maximum-random overlap algorithm which starts
   ! at the top of a column and marches down, with each layer depending on the state
   ! of the layer above it.
   !
   ! For GCM mode, changeseed must be offset between LW and SW by at least the
   ! number of subcolumns

   ! arguments
   class(ty_gas_optics_rrtmgp_ccpp), intent(in) :: kdist        ! Gas optics object
   logical,                          intent(in) :: dolw         ! Flag for whether to perform longwave calculation
   integer,                          intent(in) :: ktoprad      ! Index in RRTMGP array corresponding to top layer or interface of host model arrays
   integer,                          intent(in) :: ktopcam      ! Index in host model arrays of top level (layer or interface) at which RRTMGP is active
   integer,                          intent(in) :: nbnd         ! Number of spectral bands
   integer,                          intent(in) :: ngpt         ! Number of subcolumns (g-point intervals)
   integer,                          intent(in) :: ncol         ! Number of columns
   integer,                          intent(in) :: pver         ! Number of model layers
   integer,                          intent(in) :: nlay         ! Number of layers in radiation calculation
   integer,                          intent(in) :: changeseed   ! If the subcolumn generator is called multiple times, 
                                                                ! permute the seed between each call.
   real(kind_phys), dimension(:,:),   intent(in)  :: pmid       ! Layer pressures at midpoints (Pa)
   real(kind_phys), dimension(:,:),   intent(in)  :: cldfprime
   real(kind_phys), dimension(:,:,:), intent(in)  :: c_cld_lw_abs
   type(ty_optical_props_1scl_ccpp),  intent(inout) :: cloud_lw ! Cloud optics object
   character(len=*),                  intent(out)   :: errmsg
   integer,                           intent(out)   :: errflg

   ! Local variables

   integer :: idx, isubcol, kdx, ndx, bnd

   real(kind_phys), parameter :: cldmin = 1.0e-80_kind_phys  ! min cloud fraction
   real(kind_phys) :: cldf(ncol,nlay)      ! cloud fraction clipped to cldmin

   integer :: gpt2band(ngpt)               ! band of each g-point
   integer :: kiss_seed1, kiss_seed2, kiss_seed3, kiss_seed4
   integer(int64) :: s1, s2, s3, s4        ! KISS state of a column
   real(kind_phys) :: rand_num             ! random number (KISS)
   real(kind_phys) :: cdf(ngpt,ncol,nlay)  ! random numbers
   real(kind_phys) :: tauc(nbnd,ncol,nlay)
   real(kind_phys) :: tau_k
   !------------------------------------------------------------------------------------------

   ! Set error variables
   errflg = 0
   errmsg = ''

   ! If we're not doing longwave this timestep, no need to proceed
   if (.not. dolw) then
      return
   end if

   ! Band of each g-point (the gas optics' convert_gpt2band)
   gpt2band(:) = kdist%gas_props%get_gpoint_bands()

   ! Subset "chunk" data so just the number of CAM layers in the
   ! radiation calculation are used by MCICA to produce subcolumns;
   ! clip cloud fraction and enforce tauc >= 0.
   do kdx = 1, nlay
      do idx = 1, ncol
         cldf(idx,kdx) = cldfprime(idx,ktopcam+kdx-1)
         if (cldf(idx,kdx) < cldmin) then
            cldf(idx,kdx) = 0._kind_phys
         end if
         do bnd = 1, nbnd
            tau_k = c_cld_lw_abs(bnd,idx,ktopcam+kdx-1)
            if (tau_k > 0.0_kind_phys) then
               tauc(bnd,idx,kdx) = tau_k
            else
               tauc(bnd,idx,kdx) = 0.0_kind_phys
            end if
         end do
      end do
   end do

   ! Random numbers: one KISS stream per column, seeded from the state of
   ! the column (pmid of the bottom four layers), advanced by changeseed
   ! values, then one number per subcolumn and level (level fastest).
   do idx = 1, ncol
      kiss_seed1 = (pmid(idx,pver)   - int(pmid(idx,pver)))    * 1000000000
      kiss_seed2 = (pmid(idx,pver-1) - int(pmid(idx,pver-1)))  * 1000000000
      kiss_seed3 = (pmid(idx,pver-2) - int(pmid(idx,pver-2)))  * 1000000000
      kiss_seed4 = (pmid(idx,pver-3) - int(pmid(idx,pver-3)))  * 1000000000
      s1 = kiss_rng_seed(kiss_seed1)
      s2 = kiss_rng_seed(kiss_seed2)
      s3 = kiss_rng_seed(kiss_seed3)
      s4 = kiss_rng_seed(kiss_seed4)
      do ndx = 1, changeseed
         rand_num = kiss_rng_next(s1, s2, s3, s4)
      end do
      do isubcol = 1, ngpt
         do kdx = 1, nlay
            cdf(isubcol,idx,kdx) = kiss_rng_next(s1, s2, s3, s4)
         end do
      end do
   end do

   ! Maximum-Random overlap
   ! i) pick a random number for top layer.
   ! ii) walk down the column:
   !    - if the layer above is cloudy, use the same random number as in the layer above
   !    - if the layer above is clear, use a new random number

   do kdx = 2, nlay
      do idx = 1, ncol
         do isubcol = 1, ngpt
            if (cdf(isubcol,idx,kdx-1) > 1._kind_phys - cldf(idx,kdx-1) ) then
               cdf(isubcol,idx,kdx) = cdf(isubcol,idx,kdx-1)
            else
               cdf(isubcol,idx,kdx) = cdf(isubcol,idx,kdx) * (1._kind_phys - cldf(idx,kdx-1))
            end if
         end do
      end do
   end do

   ! If there is an extra layer in the radiation then this initialization
   ! will provide zero optical depths there
   cloud_lw%optical_props%tau = 0.0_kind_phys

   ! -- generate subcolumns for homogeneous clouds -----
   ! where there is a cloud (the subcolumn is cloudy), set the subcolumn
   ! cloud properties on the g-points; incoming tauc should be in-cloud
   ! quantites and not grid-averaged quantities
   do isubcol = 1, ngpt
      do kdx = 1, nlay
         do idx = 1, ncol
            if (cdf(isubcol,idx,kdx) >= 1._kind_phys - cldf(idx,kdx) .and. &
                cldf(idx,kdx) > 0._kind_phys) then
               cloud_lw%optical_props%tau(idx,ktoprad+kdx-1,isubcol) = &
                    tauc(gpt2band(isubcol),idx,kdx)
            else
               cloud_lw%optical_props%tau(idx,ktoprad+kdx-1,isubcol) = 0._kind_phys
            end if
         end do
      end do
   end do

   ! validate checks that: tau > 0
   errmsg = cloud_lw%optical_props%validate()
   if (len_trim(errmsg) > 0) then
       errflg = 1
       return
   end if 

end subroutine rrtmgp_lw_mcica_subcol_gen_run


end module rrtmgp_lw_mcica_subcol_gen

