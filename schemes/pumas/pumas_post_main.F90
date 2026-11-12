! Post-processing interstitial for PUMAS.
module pumas_post_main
  implicit none
  private

  public :: pumas_post_main_run

contains

  !> \section arg_table_pumas_post_main_run Argument Table
  !! \htmlinclude pumas_post_main_run.html
  subroutine pumas_post_main_run(ncol, nlev, trop_cloud_top_lev, rair, pmid, temp, &
                cldliq, numliq, cldice, numice, snow, numsnow, graupel, numgraup,  &
                strat_cldfrc, effi, dei, pgam, lamc, rel, des, degrau, errmsg, errcode)
    use ccpp_kinds,        only: kind_phys
    use pumas_kinds,       only: pumas_r8=>kind_r8
    use micro_pumas_utils, only: size_dist_param_basic, size_dist_param_liq, &
                                 mg_ice_props, mg_liq_props, avg_diameter,   &
                                 mincld, qsmall, rhog, rhosn, rhows, rhoi

    integer,                 intent(in)  :: ncol
    integer,                 intent(in)  :: nlev
    integer,                 intent(in)  :: trop_cloud_top_lev  ! Index of the top model level for which
                                                                ! tropospheric cloud physics is run (index)
    real(kind_phys),                 intent(in)  :: rair        ! gas constant of dry air (J kg-1 K-1)
    real(kind_phys), dimension(:,:), intent(in)  :: pmid        ! layer midpoint pressure (Pa)
    real(kind_phys), dimension(:,:), intent(in)  :: temp        ! air temperature, before microphysics heating is applied (K)
    real(kind_phys), dimension(:,:), intent(in)  :: cldliq      ! updated cloud liquid water mixing ratio (kg/kg)
    real(kind_phys), dimension(:,:), intent(in)  :: numliq      ! updated cloud droplet number concentration (kg-1)
    real(kind_phys), dimension(:,:), intent(in)  :: cldice      ! updated cloud ice mixing ratio (kg/kg)
    real(kind_phys), dimension(:,:), intent(in)  :: numice      ! updated cloud ice number concentration (kg-1)
    real(kind_phys), dimension(:,:), intent(in)  :: snow        ! updated snow mixing ratio (kg/kg)
    real(kind_phys), dimension(:,:), intent(in)  :: numsnow     ! updated snow number concentration (kg-1)
    real(kind_phys), dimension(:,:), intent(in)  :: graupel     ! updated graupel mixing ratio (kg/kg)
    real(kind_phys), dimension(:,:), intent(in)  :: numgraup    ! updated graupel number concentration (kg-1)
    real(kind_phys), dimension(:,:), intent(in)  :: strat_cldfrc! total stratiform cloud area fraction (= ast)
    real(kind_phys), dimension(:,:), intent(out) :: effi       ! ice effective radius (micron)
    real(kind_phys), dimension(:,:), intent(out) :: dei        ! ice effective diameter for radiation (micron)
    real(kind_phys), dimension(:,:), intent(out) :: pgam       ! droplet size distribution shape parameter for radiation (1)
    real(kind_phys), dimension(:,:), intent(out) :: lamc       ! droplet size distribution slope for radiation (1)
    real(kind_phys), dimension(:,:), intent(out) :: rel        ! liquid effective radius for radiation (micron)
    real(kind_phys), dimension(:,:), intent(out) :: des        ! snow effective diameter for radiation (micron)
    real(kind_phys), dimension(:,:), intent(out) :: degrau     ! graupel effective diameter for radiation (m)
    character(len=*),          intent(out) :: errmsg
    integer,                 intent(out) :: errcode

    integer        :: i, k
    real(pumas_r8) :: icimrst(ncol, nlev)  ! in-cloud (grid-mean) ice mixing ratio
    real(pumas_r8) :: niic(ncol, nlev)     ! in-cloud (grid-mean) ice number conc
    real(pumas_r8) :: rei(ncol, nlev)      ! ice slope param, then effective radius
    real(pumas_r8) :: icwmrst(ncol, nlev)  ! in-cloud (grid-mean) liquid mixing ratio
    real(pumas_r8) :: ncic(ncol, nlev)     ! in-cloud (grid-mean) droplet number conc
    real(pumas_r8) :: rho(ncol, nlev)      ! air density
    real(pumas_r8) :: mu(ncol, nlev)       ! droplet size distribution shape parameter
    real(pumas_r8) :: lambdac(ncol, nlev)  ! droplet size distribution slope
    real(pumas_r8) :: qs(ncol, nlev)       ! grid-mean snow mixing ratio
    real(pumas_r8) :: ns(ncol, nlev)       ! grid-mean snow number conc
    real(pumas_r8) :: qg(ncol, nlev)       ! grid-mean graupel mixing ratio
    real(pumas_r8) :: ng(ncol, nlev)       ! grid-mean graupel number conc
    real(pumas_r8) :: dsout2               ! mean snow particle diameter
    real(pumas_r8) :: desm                 ! snow effective diameter (m)
    real(pumas_r8) :: dgout2               ! mean graupel particle diameter

    errmsg = ' '
    errcode = 0

    ! One data region for the whole scheme. PUMAS's size_dist_param_basic and
    ! size_dist_param_liq run OpenACC kernels (default(present)) when PUMAS is
    ! built for the GPU, so their arguments must be on the device; the rest of
    ! the scheme runs in the two pointwise kernels below. Without -acc the
    ! directives are comments and everything runs on the host as before.
    !$acc data copyin(pmid, temp, cldliq, numliq, cldice, numice, snow, numsnow, &
    !$acc             graupel, numgraup, strat_cldfrc, mg_ice_props, mg_liq_props) &
    !$acc      copyout(effi, dei, pgam, lamc, rel, des, degrau) &
    !$acc      create(icimrst, niic, rei, icwmrst, ncic, rho, mu, lambdac, &
    !$acc             qs, ns, qg, ng)

    ! Inputs of the size distribution calculations, on every level.
    !$acc parallel loop gang vector collapse(2)
    do k = 1, nlev
      do i = 1, ncol
        ! Ice effective radius is recomputed here from the post-microphysics
        ! grid-mean in-cloud ice, NOT from the raw per-substep value PUMAS returns.
        ! PUMAS uses the total stratiform fraction for ice cloud (icecldf = ast).
        icimrst(i,k) = min(real(cldice(i,k), pumas_r8) / &
                    max(mincld, real(strat_cldfrc(i,k), pumas_r8)), 0.005_pumas_r8)
        niic(i,k)   =    real(numice(i,k), pumas_r8) / &
                    max(mincld, real(strat_cldfrc(i,k), pumas_r8))
        rei(i,k) = 25._pumas_r8

        ! Droplet size distribution parameters for radiation are likewise recomputed
        ! from the post-microphysics grid-mean in-cloud liquid.
        ! Air density uses the temperature BEFORE the microphysics heating is applied
        ! The original comment in CAM noted:
        !   "State instead of state_loc to preserve answers for MG1 (and in any
        !    case, it is unlikely to make much difference)."
        ! PUMAS uses the total stratiform fraction for liquid cloud (liqcldf = ast)
        rho(i,k)     = real(pmid(i,k), pumas_r8) / &
                    (real(rair, pumas_r8) * real(temp(i,k), pumas_r8))
        icwmrst(i,k) = min(real(cldliq(i,k), pumas_r8) / &
                    max(mincld, real(strat_cldfrc(i,k), pumas_r8)), 0.005_pumas_r8)
        ncic(i,k)   =    real(numliq(i,k), pumas_r8) / &
                    max(mincld, real(strat_cldfrc(i,k), pumas_r8))
        mu(i,k)      = 0._pumas_r8
        lambdac(i,k) = 0._pumas_r8

        ! Grid-mean snow and graupel
        qs(i,k) = real(snow(i,k), pumas_r8)
        ns(i,k) = real(numsnow(i,k), pumas_r8)
        qg(i,k) = real(graupel(i,k), pumas_r8)
        ng(i,k) = real(numgraup(i,k), pumas_r8)
      end do
    end do

    ! PUMAS's size distribution routines, one level at a time (each launches
    ! its own kernel on the GPU).
    do k = trop_cloud_top_lev, nlev
      call size_dist_param_basic(mg_ice_props, icimrst(:,k), niic(:,k), rei(:,k), ncol)
    end do
    do k = trop_cloud_top_lev, nlev
      call size_dist_param_liq(mg_liq_props, icwmrst(:,k), ncic(:,k), rho(:,k), &
                               mu(:,k), lambdac(:,k), ncol)
    end do

    ! Outputs for radiation.
    !$acc parallel loop gang vector collapse(2) private(dsout2, desm, dgout2)
    do k = 1, nlev
      do i = 1, ncol
        if (icimrst(i,k) >= qsmall) then
          rei(i,k) = 1.5_pumas_r8 / rei(i,k) * 1.e6_pumas_r8
        else
          rei(i,k) = 25._pumas_r8
        end if

        effi(i,k) = real(rei(i,k), kind_phys)

        ! Ice effective diameter for radiation follows from the recomputed
        ! effective radius.
        dei(i,k) = real(rei(i,k) * rhoi/rhows * 2._pumas_r8, kind_phys)

        ! size_dist_param_liq flags no-cloud points with mu = -100 and reset to zero.
        if (k >= trop_cloud_top_lev .and. icwmrst(i,k) < qsmall) then
          mu(i,k) = 0._pumas_r8
        end if

        pgam(i,k) = real(mu(i,k), kind_phys)
        lamc(i,k) = real(lambdac(i,k), kind_phys)

        ! Derive liquid effective radius for radiation from the recomputed size
        ! distribution and overwrite the per-substep value from PUMAS.
        ! 10 micron where there is no cloud.
        rel(i,k) = 10._kind_phys
        if (k >= trop_cloud_top_lev .and. icwmrst(i,k) >= qsmall) then
          rel(i,k) = real((mu(i,k) + 3._pumas_r8) / &
               lambdac(i,k)/2._pumas_r8 * 1.e6_pumas_r8, kind_phys)
        end if

        ! Snow effective diameter for radiation, from the post-microphysics grid-mean snow
        dsout2 = 0._pumas_r8
        desm   = 0._pumas_r8
        if (k >= trop_cloud_top_lev .and. qs(i,k) >= 1.e-7_pumas_r8) then
          dsout2 = avg_diameter(qs(i,k), ns(i,k) * rho(i,k), rho(i,k), rhosn)
          desm = dsout2 * 3._pumas_r8 * rhosn/rhows
        end if

        des(i,k) = real(desm * 1.e6_pumas_r8, kind_phys)

        ! Graupel effective diameter for radiation, from the post-microphysics grid-mean graupel
        degrau(i,k) = 0._kind_phys
        if (k >= trop_cloud_top_lev .and. qg(i,k) >= 1.e-7_pumas_r8) then
          dgout2 = avg_diameter(qg(i,k), ng(i,k) * rho(i,k), rho(i,k), rhog)
          degrau(i,k) = real(dgout2 * 3._pumas_r8 * rhog/rhows, kind_phys)
        end if
      end do
    end do

    !$acc end data

  end subroutine pumas_post_main_run

end module pumas_post_main
