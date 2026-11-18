!> \file rrtmgp_lw_cloud_optics.F90
!!

!> This module contains two routines: The first initializes data and functions
!! needed to compute the longwave cloud radiative properties in RRTMGP. The second routine
!! is a ccpp scheme within the "radiation loop", where the shortwave optical properties
!! (optical-depth, single-scattering albedo, asymmetry parameter) are computed for ALL
!! cloud types visible to RRTMGP.
module rrtmgp_lw_cloud_optics

  implicit none
  private
  public :: rrtmgp_lw_cloud_optics_run

contains

  ! SUBROUTINE rrtmgp_lw_cloud_optics_run()
  ! ######################################################################################
!> \section arg_table_rrtmgp_lw_cloud_optics_run Argument Table
!! \htmlinclude rrtmgp_lw_cloud_optics_run.html
!!
  subroutine rrtmgp_lw_cloud_optics_run(dolw, ncol, nlay, cld, cldfsnow, cldfgrau,      &
             cldfprime, kdist_lw, lamc, pgam, rei, iclwpth, iciwpth, tiny_in, dei, icswpth,  &
             des, icgrauwpth, degrau, nlwbands, do_snow, do_graupel, pver, ktopcam,     &
             cloud_lw, cld_lw_abs, snow_lw_abs, grau_lw_abs, c_cld_lw_abs, errmsg, errflg)
    use ccpp_gas_optics_rrtmgp,    only: ty_gas_optics_rrtmgp_ccpp
    use ccpp_optical_props,        only: ty_optical_props_1scl_ccpp
    use ccpp_kinds,                only: kind_phys
    use mo_rte_util_array,         only: zero_array
    use rrtmgp_cloud_optics_setup, only: g_mu, g_lambda, nmu, nlambda, g_d_eff, n_g_d
    use rrtmgp_cloud_optics_setup, only: abs_lw_liq, abs_lw_ice
    use rrtmgp_cloud_optics_setup, only: liq_cld_optics, ice_cld_optics
    ! Compute combined cloud optical properties
    ! Create MCICA stochastic arrays for cloud LW optical properties
    ! Initialize optical properties object (cloud_lw) and load with MCICA columns

    ! Inputs
    integer,                           intent(in) :: ncol             ! Number of columns
    integer,                           intent(in) :: nlay             ! Number of vertical layers in radiation
    integer,                           intent(in) :: nlwbands         ! Number of longwave bands
    integer,                           intent(in) :: pver             ! Total number of vertical layers
    integer,                           intent(in) :: ktopcam          ! Index in host model arrays of top level (layer or interface) at which RRTMGP is active
    real(kind_phys), dimension(:,:),   intent(in) :: cld              ! Cloud fraction (liq + ice)
    real(kind_phys), dimension(:,:),   intent(in) :: cldfsnow         ! Cloud fraction of just "snow clouds"
    real(kind_phys), dimension(:,:),   intent(in) :: cldfgrau         ! Cloud fraction of just "graupel clouds"
    real(kind_phys), dimension(:,:),   intent(in) :: cldfprime        ! Modified cloud fraction
    real(kind_phys), dimension(:,:),   intent(in) :: lamc             ! Prognosed value of lambda for cloud
    real(kind_phys), dimension(:,:),   intent(in) :: pgam             ! Prognosed value of mu for cloud
    real(kind_phys), dimension(:,:),   intent(in) :: rei              ! Effective radius of stratiform cloud ice crystal [um]
    real(kind_phys), dimension(:,:),   intent(in) :: iclwpth          ! In-cloud liquid water path
    real(kind_phys), dimension(:,:),   intent(in) :: iciwpth          ! In-cloud ice water path 
    real(kind_phys), dimension(:,:),   intent(in) :: icswpth          ! In-cloud snow water path
    real(kind_phys), dimension(:,:),   intent(in) :: icgrauwpth       ! In-cloud graupel water path
    real(kind_phys), dimension(:,:),   intent(in) :: dei              ! Mean effective radius for ice cloud
    real(kind_phys), dimension(:,:),   intent(in) :: des              ! Mean effective radius for snow
    real(kind_phys), dimension(:,:),   intent(in) :: degrau           ! Mean effective radius for graupel
    real(kind_phys),                   intent(in) :: tiny_in          ! Definition of tiny for RRTMGP
    logical,                           intent(in) :: do_snow          ! Flag for whether cldfsnow is present
    logical,                           intent(in) :: do_graupel       ! Flag for whether cldfgrau is present
    logical,                           intent(in) :: dolw             ! Flag for whether to perform longwave calculation
    class(ty_gas_optics_rrtmgp_ccpp),  intent(in) :: kdist_lw         ! Longwave gas optics object

    ! Outputs
    type(ty_optical_props_1scl_ccpp),  intent(out) :: cloud_lw        ! Longwave cloud optics object
    real(kind_phys), dimension(:,:,:), intent(out) :: cld_lw_abs      ! Cloud absorption optics depth (LW)
    real(kind_phys), dimension(:,:,:), intent(out) :: snow_lw_abs     ! Snow absorption optics depth (LW)
    real(kind_phys), dimension(:,:,:), intent(out) :: grau_lw_abs     ! Graupel absorption optics depth (LW)
    real(kind_phys), dimension(:,:,:), intent(out) :: c_cld_lw_abs
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    ! Local variables
    integer :: idx, kdx, lwband

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(kind_phys) :: liq_lw_abs(nlwbands, ncol, pver)   ! liquid absorption optics depth (LW)
    real(kind_phys) :: ice_lw_abs(nlwbands, ncol, pver)   ! ice absorption optics depth (LW)

    character(len=*), parameter :: sub = 'rrtmgp_lw_cloud_optics_run'
    !--------------------------------------------------------------------------------

    ! Set error variables
    errmsg = ''
    errflg = 0

    ! If not doing longwave, no need to proceed
    if (.not. dolw) then
       return
    end if

    ! Check the parameterization selectors on the host, before the data region
    if (trim(liq_cld_optics) /= 'slingo' .and. trim(liq_cld_optics) /= 'gammadist') then
       write(errmsg,'(a,a)') sub, ': liq_cld_optics must be either slingo or gammadist'
       errflg = 1
       return
    end if
    if (trim(ice_cld_optics) /= 'ebertcurry' .and. trim(ice_cld_optics) /= 'mitchell') then
       write(errmsg,'(a,a)') sub, ': ice_cld_optics must be either ebertcurry or mitchell'
       errflg = 1
       return
    end if

    !$acc data copyin(cld, cldfsnow, cldfgrau, cldfprime, lamc, pgam, rei, iclwpth, iciwpth, &
    !$acc             dei, icswpth, des, icgrauwpth, degrau)                               &
    !$acc      copyout(cld_lw_abs, snow_lw_abs, grau_lw_abs, c_cld_lw_abs)                 &
    !$acc      create(liq_lw_abs, ice_lw_abs)

    ! Combine the cloud optical properties.
    ! Note: mix-match slingo/ebertcurry with gammadist/mitchell is NOT supported.
    !       the CAM4-era slingo/ebertcurry longwave absorptions each apply
    !       their coefficient to the TOTAL (liquid+ice) water path weighted by the
    !       ice fraction, so their sum reproduces the CAM4 total only when used as a pair.
    if (trim(liq_cld_optics) == 'slingo') then
       ! Slingo liquid optics (CAM4-era broadband liquid absorption)
       call slingo_liq_get_rad_props_lw(ncol, pver, nlwbands, iclwpth, iciwpth, liq_lw_abs)
    else
       ! gammadist liquid optics
       call liquid_cloud_get_rad_props_lw(ncol, pver, nmu, nlambda, nlwbands, lamc, pgam, g_mu, g_lambda, iclwpth, &
               abs_lw_liq, tiny_in, liq_lw_abs, errmsg, errflg)
    end if

    if (errflg == 0) then
       if (trim(ice_cld_optics) == 'ebertcurry') then
          ! Ebert and Curry (1992) ice optics
          call ec_ice_get_rad_props_lw(ncol, pver, nlwbands, rei, iclwpth, iciwpth, ice_lw_abs)
       else
          ! Mitchell ice optics
          call interpolate_ice_optics_lw(ncol, pver, nlwbands, iciwpth, dei, &
                  n_g_d, g_d_eff, abs_lw_ice, tiny_in, ice_lw_abs, errmsg, errflg)
       end if
    end if

    if (errflg == 0) then
       !$acc parallel loop gang vector collapse(3) default(present)
       do kdx = 1, pver
          do idx = 1, ncol
             do lwband = 1, nlwbands
                cld_lw_abs(lwband,idx,kdx) = liq_lw_abs(lwband,idx,kdx) + ice_lw_abs(lwband,idx,kdx)
             end do
          end do
       end do
    end if

    ! add in snow
    if (errflg == 0) then
       if (do_snow) then
          call interpolate_ice_optics_lw(ncol, pver, nlwbands, icswpth, des, &
                  n_g_d, g_d_eff, abs_lw_ice, tiny_in, snow_lw_abs, errmsg, errflg)
          if (errflg == 0) then
             !$acc parallel loop gang vector collapse(3) default(present)
             do kdx = 1, pver
                do idx = 1, ncol
                   do lwband = 1, nlwbands
                      if (cldfprime(idx,kdx) > 0._kind_phys) then
                         c_cld_lw_abs(lwband,idx,kdx) = ( cldfsnow(idx,kdx)*snow_lw_abs(lwband,idx,kdx) &
                                                        + cld(idx,kdx)*cld_lw_abs(lwband,idx,kdx) )/cldfprime(idx,kdx)
                      else
                         c_cld_lw_abs(lwband,idx,kdx) = 0._kind_phys
                      end if
                   end do
                end do
             end do
          end if
       else
          !$acc parallel loop gang vector collapse(3) default(present)
          do kdx = 1, pver
             do idx = 1, ncol
                do lwband = 1, nlwbands
                   snow_lw_abs(lwband,idx,kdx) = 0._kind_phys
                   c_cld_lw_abs(lwband,idx,kdx) = cld_lw_abs(lwband,idx,kdx)
                end do
             end do
          end do
       end if
    end if

    ! add in graupel
    if (errflg == 0) then
       if (do_graupel) then
          call interpolate_ice_optics_lw(ncol, pver, nlwbands, icgrauwpth, degrau, n_g_d, &
                  g_d_eff, abs_lw_ice, tiny_in, grau_lw_abs, errmsg, errflg)
          if (errflg == 0) then
             !$acc parallel loop gang vector collapse(3) default(present)
             do kdx = 1, pver
                do idx = 1, ncol
                   do lwband = 1, nlwbands
                      if (cldfprime(idx,kdx) > 0._kind_phys) then
                         c_cld_lw_abs(lwband,idx,kdx) = ( cldfgrau(idx,kdx)*grau_lw_abs(lwband,idx,kdx) &
                                                        + cld(idx,kdx)*c_cld_lw_abs(lwband,idx,kdx) )/cldfprime(idx,kdx)
                      else
                         c_cld_lw_abs(lwband,idx,kdx) = 0._kind_phys
                      end if
                   end do
                end do
             end do
          end if
       else
          !$acc parallel loop gang vector collapse(3) default(present)
          do kdx = 1, pver
             do idx = 1, ncol
                do lwband = 1, nlwbands
                   grau_lw_abs(lwband,idx,kdx) = 0._kind_phys
                end do
             end do
          end do
       end if
    end if

    !$acc end data
    if (errflg /= 0) then
       return
    end if

    errmsg =cloud_lw%optical_props%alloc_1scl(ncol, nlay, kdist_lw%gas_props)
    if (len_trim(errmsg) > 0) then
       errflg = 1
       return
    end if
    ! Keep tau on the device for the whole radiation call: created here,
    ! deleted in rrtmgp_post (rrtmgp-rte-object-residency). Defined at
    ! return; rrtmgp_lw_mcica_subcol_gen fills it.
    !$acc enter data copyin(cloud_lw%optical_props)
    !$acc enter data create(cloud_lw%optical_props%tau)
    call zero_array(ncol, nlay, cloud_lw%optical_props%get_ngpt(), cloud_lw%optical_props%tau)

  end subroutine rrtmgp_lw_cloud_optics_run

!==============================================================================

  subroutine liquid_cloud_get_rad_props_lw(ncol, pver, nmu, nlambda, nlwbands, lamc, pgam, &
                  g_mu, g_lambda, iclwpth, abs_lw_liq, tiny, abs_od, errmsg, errflg)
    use ccpp_kinds,              only: kind_phys
    use atmos_phys_interp_point, only: interp_point_wgts, interp_point_wgts_blend, &
                                       interp_point_2d, interp_point_errmsg
    ! Inputs
    integer,                           intent(in) :: ncol
    integer,                           intent(in) :: pver
    integer,                           intent(in) :: nmu
    integer,                           intent(in) :: nlambda
    integer,                           intent(in) :: nlwbands
    real(kind_phys), dimension(:,:),   intent(in) :: lamc
    real(kind_phys), dimension(:,:),   intent(in) :: pgam
    real(kind_phys), dimension(:,:,:), intent(in) :: abs_lw_liq
    real(kind_phys), dimension(:),     intent(in) :: g_mu
    real(kind_phys), dimension(:,:),   intent(in) :: g_lambda
    real(kind_phys), dimension(:,:),   intent(in) :: iclwpth
    real(kind_phys),                   intent(in) :: tiny
    ! Outputs
    real(kind_phys), dimension(:,:,:), intent(out) :: abs_od
    character(len=*),                  intent(out) :: errmsg
    integer,                           intent(out) :: errflg

    integer :: lwband, idx, kdx
    integer :: mu_m, mu_p, lam_m, lam_p, ierr, ierr_max
    real(kind_phys) :: mu_ws, mu_wn, lam_ws, lam_wn

    ! Set error variables
    errflg = 0
    errmsg = ''

    ! For each cell: weights in mu (pgam), then in lambda (lamc) on the
    ! lambda row interpolated to that mu, then the bilinear table lookup
    ! (as get_mu_lambda_weights_ccpp and lininterp do).
    ierr_max = 0
    !$acc parallel loop gang vector collapse(2) default(present) reduction(max:ierr_max) &
    !$acc   private(mu_m, mu_p, lam_m, lam_p, ierr, mu_ws, mu_wn, lam_ws, lam_wn)
    do kdx = 1,pver
       do idx = 1,ncol
          ! lamc <= 0 seems to be the clue for no cloud from microphysics formulation
          if (lamc(idx,kdx) > 0._kind_phys .and. .not. (iclwpth(idx,kdx) < tiny)) then
             call interp_point_wgts(g_mu, nmu, pgam(idx,kdx), mu_m, mu_p, mu_ws, mu_wn, ierr)
             if (ierr == 0) then
                call interp_point_wgts_blend(g_lambda, nlambda, mu_m, mu_p, mu_ws, mu_wn, &
                        lamc(idx,kdx), lam_m, lam_p, lam_ws, lam_wn, ierr)
             end if
             ierr_max = max(ierr_max, ierr)
             if (ierr == 0) then
                !$acc loop seq
                do lwband = 1, nlwbands
                   abs_od(lwband,idx,kdx) = iclwpth(idx,kdx) * interp_point_2d(              &
                        abs_lw_liq(:,:,lwband), mu_m, mu_p, mu_ws, mu_wn, &
                        lam_m, lam_p, lam_ws, lam_wn)
                end do
             else
                !$acc loop seq
                do lwband = 1, nlwbands
                   abs_od(lwband,idx,kdx) = 0._kind_phys
                end do
             end if
          else
             !$acc loop seq
             do lwband = 1, nlwbands
                abs_od(lwband,idx,kdx) = 0._kind_phys
             end do
          endif
       enddo
    enddo
    if (ierr_max /= 0) then
       errmsg = interp_point_errmsg(ierr_max)
       errflg = 1
    end if

  end subroutine liquid_cloud_get_rad_props_lw

!==============================================================================

  subroutine interpolate_ice_optics_lw(ncol, pver, nlwbands, iciwpth, dei, &
                  n_g_d, g_d_eff, abs_lw_ice, tiny, abs_od, errmsg, errflg)
    use atmos_phys_interp_point,  only: interp_point_wgts, interp_point_1d, interp_point_errmsg
    use ccpp_kinds,               only: kind_phys

    integer,           intent(in)                  :: ncol
    integer,           intent(in)                  :: n_g_d
    integer,           intent(in)                  :: pver
    integer,           intent(in)                  :: nlwbands
    real(kind_phys), dimension(:),     intent(in)  :: g_d_eff
    real(kind_phys), dimension(:,:),   intent(in)  :: iciwpth
    real(kind_phys), dimension(:,:),   intent(in)  :: dei
    real(kind_phys), dimension(:,:),   intent(in)  :: abs_lw_ice
    real(kind_phys),                   intent(in)  :: tiny
    real(kind_phys), dimension(:,:,:), intent(out) :: abs_od
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: idx, kdx, lwband
    integer :: jjm, jjp, ierr, ierr_max
    real(kind_phys) :: wgts, wgtn

    ! Set error variables
    errflg = 0
    errmsg = ''

    ierr_max = 0
    !$acc parallel loop gang vector collapse(2) default(present) reduction(max:ierr_max) &
    !$acc   private(jjm, jjp, ierr, wgts, wgtn)
    do kdx = 1,pver
       do idx = 1,ncol
          ! if ice water path is too small, OD := 0
          if( iciwpth(idx,kdx) < tiny .or. dei(idx,kdx) == 0._kind_phys) then
             !$acc loop seq
             do lwband = 1, nlwbands
                abs_od(lwband,idx,kdx) = 0._kind_phys
             end do
          else
             ! for each cell interpolate to find weights in g_d_eff grid.
             call interp_point_wgts(g_d_eff, n_g_d, dei(idx,kdx), jjm, jjp, wgts, wgtn, ierr)
             ierr_max = max(ierr_max, ierr)
             if (ierr == 0) then
                ! interpolate into grid and extract radiative properties
                !$acc loop seq
                do lwband = 1, nlwbands
                   abs_od(lwband,idx,kdx) = iciwpth(idx,kdx) * &
                        interp_point_1d(abs_lw_ice(:,lwband), jjm, jjp, wgts, wgtn)
                   if (abs_od(lwband,idx,kdx) > 50.0_kind_phys) then
                      abs_od(lwband,idx,kdx) = 50.0_kind_phys
                   end if
                enddo
             else
                !$acc loop seq
                do lwband = 1, nlwbands
                   abs_od(lwband,idx,kdx) = 0._kind_phys
                end do
             end if
          endif
       enddo
    enddo
    if (ierr_max /= 0) then
       errmsg = interp_point_errmsg(ierr_max)
       errflg = 1
    end if

  end subroutine interpolate_ice_optics_lw

!==============================================================================

!==============================================================================

  subroutine slingo_liq_get_rad_props_lw(ncol, pver, nlwbands, iclwpth, iciwpth, abs_od)
    ! Kiehl et al. (1998) (https://doi.org/10.1175/1520-0442(1998)011<1131:TNCFAR>2.0.CO;2) longwave liquid
    ! absorption (broadband, CAM4), using the in-cloud water paths.

    use ccpp_kinds, only: kind_phys
    ! Inputs
    integer,                           intent(in) :: ncol
    integer,                           intent(in) :: pver
    integer,                           intent(in) :: nlwbands
    real(kind_phys), dimension(:,:),   intent(in) :: iclwpth   ! In-cloud liquid water path [kg m-2]
    real(kind_phys), dimension(:,:),   intent(in) :: iciwpth   ! In-cloud ice water path [kg m-2]
    ! Outputs
    real(kind_phys), dimension(:,:,:), intent(out) :: abs_od

    ! Local variables
    real(kind_phys) :: ficemr
    real(kind_phys) :: cwp
    real(kind_phys) :: cldtau
    real(kind_phys) :: kabs
    integer :: lwband, i, k

    real(kind_phys), parameter :: kabsl = 0.090361_kind_phys ! longwave liquid absorption coeff (m**2/g)

    !$acc parallel loop gang vector collapse(2) default(present) private(cwp, ficemr, kabs, cldtau)
    do k=1, pver
       do i = 1, ncol
          cwp    = 1000.0_kind_phys * iclwpth(i,k) + 1000.0_kind_phys * iciwpth(i, k)
          ficemr = 1000.0_kind_phys * iciwpth(i,k)/(max(1.e-18_kind_phys, cwp))
          kabs = kabsl*(1._kind_phys-ficemr)
          cldtau = kabs*cwp
          !$acc loop seq
          do lwband = 1, nlwbands
             abs_od(lwband,i,k) = cldtau
          end do
       end do
    end do

  end subroutine slingo_liq_get_rad_props_lw

!==============================================================================

  subroutine ec_ice_get_rad_props_lw(ncol, pver, nlwbands, rei, iclwpth, iciwpth, abs_od)
    ! Ebert and Curry (1992) (https://doi.org/10.1029/91JD02472) longwave ice absorption (broadband, CAM4),
    ! using the in-cloud water paths.

    use ccpp_kinds, only: kind_phys
    ! Inputs
    integer,                           intent(in) :: ncol
    integer,                           intent(in) :: pver
    integer,                           intent(in) :: nlwbands
    real(kind_phys), dimension(:,:),   intent(in) :: rei       ! Effective radius of stratiform cloud ice crystal [um]
    real(kind_phys), dimension(:,:),   intent(in) :: iclwpth   ! In-cloud liquid water path [kg m-2]
    real(kind_phys), dimension(:,:),   intent(in) :: iciwpth   ! In-cloud ice water path [kg m-2]
    ! Outputs
    real(kind_phys), dimension(:,:,:), intent(out) :: abs_od

    ! Local variables
    real(kind_phys) :: ficemr
    real(kind_phys) :: cwp
    real(kind_phys) :: cldtau
    real(kind_phys) :: kabs, kabsi
    integer :: lwband, i, k

    real(kind_phys), parameter :: scalefactor = 1._kind_phys

    !$acc parallel loop gang vector collapse(2) default(present) private(cwp, ficemr, kabsi, kabs, cldtau)
    do k=1,pver
       do i = 1,ncol
          cwp    = 1000.0_kind_phys * iciwpth(i,k) + 1000.0_kind_phys * iclwpth(i,k)
          ficemr = 1000.0_kind_phys * iciwpth(i,k)/(max(1.e-18_kind_phys, cwp))
          ! note that optical properties for ice valid only
          ! in range of 13 > rei > 130 micron
          kabsi = 0.005_kind_phys + 1._kind_phys/min(max(13._kind_phys,scalefactor*rei(i,k)),130._kind_phys)
          kabs =  kabsi*ficemr ! kabsl*(1._r8-ficemr(i,k)) + kabsi*ficemr(i,k)
          cldtau = kabs*cwp
          !$acc loop seq
          do lwband = 1,nlwbands
             abs_od(lwband,i,k) = cldtau
          end do
       end do
    end do

  end subroutine ec_ice_get_rad_props_lw

!==============================================================================

end module rrtmgp_lw_cloud_optics
