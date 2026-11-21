! Aerosol optical properties of the bulk aerosols (BAM) for RRTMGP radiation,
! from the tables of bulk_aerosol_optics_tables. Same results as the
! aerosol_optics scheme for a climate list of bulk aerosols only (its
! hygroscopic, insoluble and volcanic_radius optics; aerosol_optics_core and
! the *_aerosol_optics_mod classes), without the host's aerosol objects.
! Every expression keeps the order of the original so the answers are
! bit-for-bit.
!
! Adapted from aerosol_optics.F90 (Haipeng Lin, NSF-NCAR/CGD/AMP, March 2026)
! and CAM-SIMA's src/aerosol optics modules (aerosol_optics_core,
! hygroscopic_, insoluble_ and volcrad_aerosol_optics_mod, aerosol_mmr_host).
module bulk_aerosol_optics
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: bulk_aerosol_optics_run

  ! As in aerosol_optics
  real(kind_phys), parameter :: aodvis_min     = 1.e-10_kind_phys
  real(kind_phys), parameter :: ssavis_default = 0.925_kind_phys

  ! As in bulk_aerosol_optics_tables
  integer, parameter :: optics_hygroscopic = 1
  integer, parameter :: optics_insoluble   = 2
  integer, parameter :: optics_volcanic    = 3
  integer, parameter :: class_dust      = 1
  integer, parameter :: class_black_c   = 2
  integer, parameter :: class_sulfate   = 3
  integer, parameter :: class_p_organic = 4
  integer, parameter :: class_s_organic = 5
  integer, parameter :: class_seasalt   = 6

contains

!> \section arg_table_bulk_aerosol_optics_run Argument Table
!! \htmlinclude bulk_aerosol_optics_run.html
  subroutine bulk_aerosol_optics_run( &
    ncol, pver, nswbands, nlwbands, top_lev, rga, idx_sw_diag, &
    num_bulk_aer, nsamples, &
    relh, pdeldry, constituents, &
    optics_method, idx_mmr, idx_radius, species_class, nsamples_aer, &
    density, hygro, refindex_re, refindex_im, &
    sw_ext, sw_ssa, sw_asm, lw_abs, log_radius, &
    aer_tau, aer_tau_w, aer_tau_w_g, aer_lw_abs, &
    dustaod, sulfaod, bcaod, pomaod, soaaod, ssltaod, &
    aodabsbc, &
    burdendust, burdenso4, burdenbc, burdenpom, burdensoa, burdenseasalt, &
    ssavis, aodvis, &
    odv_col_aod, &
    errmsg, errflg)

    ! Input arguments
    integer,         intent(in)  :: ncol
    integer,         intent(in)  :: pver
    integer,         intent(in)  :: nswbands          ! number of SW bands [count]
    integer,         intent(in)  :: nlwbands          ! number of LW bands [count]
    integer,         intent(in)  :: top_lev           ! top level for aerosol model [index]
    real(kind_phys), intent(in)  :: rga               ! 1/g [s2 m-1]
    integer,         intent(in)  :: idx_sw_diag       ! index of SW diagnostic band [index]
    integer,         intent(in)  :: num_bulk_aer      ! number of bulk aerosols [count]
    integer,         intent(in)  :: nsamples          ! samples in the tables [count]
    real(kind_phys), intent(in)  :: relh(:, :)        ! relative humidity [fraction]
    real(kind_phys), intent(in)  :: pdeldry(:, :)     ! dry air pressure thickness [Pa]
    real(kind_phys), intent(in)  :: constituents(:, :, :)
    ! Tables (bulk_aerosol_optics_tables)
    integer,         intent(in)  :: optics_method(:)
    integer,         intent(in)  :: idx_mmr(:)
    integer,         intent(in)  :: idx_radius(:)
    integer,         intent(in)  :: species_class(:)
    integer,         intent(in)  :: nsamples_aer(:)
    real(kind_phys), intent(in)  :: density(:)
    real(kind_phys), intent(in)  :: hygro(:)
    real(kind_phys), intent(in)  :: refindex_re(:, :)
    real(kind_phys), intent(in)  :: refindex_im(:, :)
    real(kind_phys), intent(in)  :: sw_ext(:, :, :)
    real(kind_phys), intent(in)  :: sw_ssa(:, :, :)
    real(kind_phys), intent(in)  :: sw_asm(:, :, :)
    real(kind_phys), intent(in)  :: lw_abs(:, :, :)
    real(kind_phys), intent(in)  :: log_radius(:, :)

    ! Output arguments
    real(kind_phys), intent(out) :: aer_tau(:, :, :)      ! SW extinction optical depth [1]
    real(kind_phys), intent(out) :: aer_tau_w(:, :, :)    ! SW ssa * tau [1]
    real(kind_phys), intent(out) :: aer_tau_w_g(:, :, :)  ! SW asy * ssa * tau [1]
    real(kind_phys), intent(out) :: aer_lw_abs(:, :, :)   ! LW absorption optical depth [1]
    real(kind_phys), intent(out) :: dustaod(:)       ! Dust AOD at vis [1]
    real(kind_phys), intent(out) :: sulfaod(:)       ! Sulfate AOD at vis [1]
    real(kind_phys), intent(out) :: bcaod(:)         ! BC AOD at vis [1]
    real(kind_phys), intent(out) :: pomaod(:)        ! POM AOD at vis [1]
    real(kind_phys), intent(out) :: soaaod(:)        ! SOA AOD at vis [1]
    real(kind_phys), intent(out) :: ssltaod(:)       ! Seasalt AOD at vis [1]
    real(kind_phys), intent(out) :: aodabsbc(:)      ! BC absorption OD at vis [1]
    real(kind_phys), intent(out) :: burdendust(:)    ! Dust column burden [kg m-2]
    real(kind_phys), intent(out) :: burdenso4(:)     ! Sulfate column burden [kg m-2]
    real(kind_phys), intent(out) :: burdenbc(:)      ! BC column burden [kg m-2]
    real(kind_phys), intent(out) :: burdenpom(:)     ! POM column burden [kg m-2]
    real(kind_phys), intent(out) :: burdensoa(:)     ! SOA column burden [kg m-2]
    real(kind_phys), intent(out) :: burdenseasalt(:) ! Seasalt column burden [kg m-2]
    real(kind_phys), intent(out) :: ssavis(:)        ! Column SSA at vis [1]
    real(kind_phys), intent(out) :: aodvis(:)        ! Column AOD at vis [1]
    real(kind_phys), intent(out) :: odv_col_aod(:, :) ! column visible OD per bulk aerosol [1]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local variables
    integer         :: icol, ilev, iwav, ibin
    real(kind_phys) :: mass(ncol, pver)
    ! SW optics of each aerosol in the diagnostic band
    real(kind_phys) :: tau_vis(ncol, pver, num_bulk_aer)
    real(kind_phys) :: ssa_vis(ncol, pver, num_bulk_aer)
    real(kind_phys) :: tau, ssa, asm, tau_lw
    real(kind_phys) :: odsum, specmmr, radius
    ! Per-species diagnostics
    integer         :: icls
    real(kind_phys) :: vol, wetvol, watervol
    real(kind_phys) :: scat_dust, abs_dust, hyg_dust, scat_bc, abs_bc, hyg_bc
    real(kind_phys) :: scat_sulf, abs_sulf, hyg_sulf, scat_pom, abs_pom, hyg_pom
    real(kind_phys) :: scat_soa, abs_soa, hyg_soa, scat_sslt, abs_sslt, hyg_sslt
    real(kind_phys) :: scath2o, absh2o, sumscat, sumabs, sumhygro
    real(kind_phys) :: dopaer, palb, aodc

    errmsg = ''
    errflg = 0

    ! A table index of idx_sw_diag would be out of range otherwise
    if (num_bulk_aer > 0 .and. (idx_sw_diag < 1 .or. idx_sw_diag > nswbands)) then
      errflg = 1
      errmsg = 'bulk_aerosol_optics_run: idx_sw_diag out of range'
      return
    end if

    do ilev = 1, pver
      do icol = 1, ncol
        mass(icol, ilev) = pdeldry(icol, ilev) * rga
      end do
    end do

    ! SW: sum over the aerosols in list order, at each (column, level, band),
    ! from zero, as aerosol_optics accumulates its per-aerosol optics
    do iwav = 1, nswbands
      do ilev = 1, pver
        do icol = 1, ncol
          aer_tau(icol, ilev, iwav)     = 0._kind_phys
          aer_tau_w(icol, ilev, iwav)   = 0._kind_phys
          aer_tau_w_g(icol, ilev, iwav) = 0._kind_phys
          if (ilev >= top_lev) then
            do ibin = 1, num_bulk_aer
              specmmr = aerosol_mmr(constituents, idx_mmr(ibin), icol, ilev)
              radius  = geometric_radius(constituents, idx_radius(ibin), icol, ilev)
              call sw_optics(optics_method(ibin), nsamples_aer(ibin), ibin, iwav, &
                             sw_ext, sw_ssa, sw_asm, log_radius, &
                             relh(icol, ilev), radius, specmmr, mass(icol, ilev), &
                             tau, ssa, asm)
              aer_tau(icol, ilev, iwav)     = aer_tau(icol, ilev, iwav) + tau
              aer_tau_w(icol, ilev, iwav)   = aer_tau_w(icol, ilev, iwav) + tau * ssa
              aer_tau_w_g(icol, ilev, iwav) = aer_tau_w_g(icol, ilev, iwav) + tau * ssa * asm
              if (iwav == idx_sw_diag) then
                tau_vis(icol, ilev, ibin) = tau
                ssa_vis(icol, ilev, ibin) = ssa
              end if
            end do
          end if
        end do
      end do
    end do

    ! LW: all levels
    do iwav = 1, nlwbands
      do ilev = 1, pver
        do icol = 1, ncol
          aer_lw_abs(icol, ilev, iwav) = 0._kind_phys
          do ibin = 1, num_bulk_aer
            specmmr = aerosol_mmr(constituents, idx_mmr(ibin), icol, ilev)
            radius  = geometric_radius(constituents, idx_radius(ibin), icol, ilev)
            call lw_optics(optics_method(ibin), nsamples_aer(ibin), ibin, iwav, &
                           lw_abs, log_radius, relh(icol, ilev), radius, specmmr, &
                           mass(icol, ilev), tau_lw)
            aer_lw_abs(icol, ilev, iwav) = aer_lw_abs(icol, ilev, iwav) + tau_lw
          end do
        end do
      end do
    end do

    ! Column diagnostics (per-species AOD and burdens): aerosols in list
    ! order, then levels, as aerosol_optics. A bulk aerosol is one species
    ! of one class, so five of the six class terms are zero; they are kept
    ! so every sum has the original terms. Bulk aerosols hold no water
    ! (bulk_aerosol_state water_volume = 0): the wet volume is the dry
    ! volume and the water terms are zero (the water refractive index
    ! multiplies that zero volume).
    do icol = 1, ncol
      dustaod(icol)       = 0._kind_phys
      sulfaod(icol)       = 0._kind_phys
      bcaod(icol)         = 0._kind_phys
      pomaod(icol)        = 0._kind_phys
      soaaod(icol)        = 0._kind_phys
      ssltaod(icol)       = 0._kind_phys
      aodabsbc(icol)      = 0._kind_phys
      burdendust(icol)    = 0._kind_phys
      burdenso4(icol)     = 0._kind_phys
      burdenbc(icol)      = 0._kind_phys
      burdenpom(icol)     = 0._kind_phys
      burdensoa(icol)     = 0._kind_phys
      burdenseasalt(icol) = 0._kind_phys
      do ibin = 1, num_bulk_aer
        icls = species_class(ibin)
        do ilev = top_lev, pver
          specmmr = aerosol_mmr(constituents, idx_mmr(ibin), icol, ilev)
          vol     = specmmr / density(ibin)
          scat_dust = 0._kind_phys; abs_dust = 0._kind_phys; hyg_dust = 0._kind_phys
          scat_bc   = 0._kind_phys; abs_bc   = 0._kind_phys; hyg_bc   = 0._kind_phys
          scat_sulf = 0._kind_phys; abs_sulf = 0._kind_phys; hyg_sulf = 0._kind_phys
          scat_pom  = 0._kind_phys; abs_pom  = 0._kind_phys; hyg_pom  = 0._kind_phys
          scat_soa  = 0._kind_phys; abs_soa  = 0._kind_phys; hyg_soa  = 0._kind_phys
          scat_sslt = 0._kind_phys; abs_sslt = 0._kind_phys; hyg_sslt = 0._kind_phys
          select case (icls)
          case (class_dust)
            burdendust(icol) = burdendust(icol) + specmmr*mass(icol, ilev)
            scat_dust = vol*refindex_re(idx_sw_diag, ibin)
            abs_dust  = -vol*refindex_im(idx_sw_diag, ibin)
            hyg_dust  = vol*hygro(ibin)
          case (class_black_c)
            burdenbc(icol) = burdenbc(icol) + specmmr*mass(icol, ilev)
            scat_bc = vol*refindex_re(idx_sw_diag, ibin)
            abs_bc  = -vol*refindex_im(idx_sw_diag, ibin)
            hyg_bc  = vol*hygro(ibin)
          case (class_sulfate)
            burdenso4(icol) = burdenso4(icol) + specmmr*mass(icol, ilev)
            scat_sulf = vol*refindex_re(idx_sw_diag, ibin)
            abs_sulf  = -vol*refindex_im(idx_sw_diag, ibin)
            hyg_sulf  = vol*hygro(ibin)
          case (class_p_organic)
            burdenpom(icol) = burdenpom(icol) + specmmr*mass(icol, ilev)
            scat_pom = vol*refindex_re(idx_sw_diag, ibin)
            abs_pom  = -vol*refindex_im(idx_sw_diag, ibin)
            hyg_pom  = vol*hygro(ibin)
          case (class_s_organic)
            burdensoa(icol) = burdensoa(icol) + specmmr*mass(icol, ilev)
            scat_soa = vol*refindex_re(idx_sw_diag, ibin)
            abs_soa  = -vol*refindex_im(idx_sw_diag, ibin)
            hyg_soa  = vol*hygro(ibin)
          case (class_seasalt)
            burdenseasalt(icol) = burdenseasalt(icol) + specmmr*mass(icol, ilev)
            scat_sslt = vol*refindex_re(idx_sw_diag, ibin)
            abs_sslt  = -vol*refindex_im(idx_sw_diag, ibin)
            hyg_sslt  = vol*hygro(ibin)
          end select

          watervol = 0._kind_phys
          wetvol   = vol + watervol
          if (wetvol > 1.e-40_kind_phys .and. vol > 0._kind_phys) then
            scath2o  = watervol
            absh2o   = -watervol
            sumscat  = scat_sulf + scat_pom + scat_soa + scat_bc + scat_dust + scat_sslt + scath2o
            sumabs   = abs_sulf + abs_pom + abs_soa + abs_bc + abs_dust + abs_sslt + absh2o
            sumhygro = hyg_sulf + hyg_pom + hyg_soa + hyg_bc + hyg_dust + hyg_sslt
            if (sumscat > 0._kind_phys .and. sumabs > 0._kind_phys .and. sumhygro > 0._kind_phys) then
              scat_dust = (scat_dust + scath2o*hyg_dust/sumhygro)/sumscat
              abs_dust  = (abs_dust + absh2o*hyg_dust/sumhygro)/sumabs
              scat_sulf = (scat_sulf + scath2o*hyg_sulf/sumhygro)/sumscat
              abs_sulf  = (abs_sulf + absh2o*hyg_sulf/sumhygro)/sumabs
              scat_pom  = (scat_pom + scath2o*hyg_pom/sumhygro)/sumscat
              abs_pom   = (abs_pom + absh2o*hyg_pom/sumhygro)/sumabs
              scat_soa  = (scat_soa + scath2o*hyg_soa/sumhygro)/sumscat
              abs_soa   = (abs_soa + absh2o*hyg_soa/sumhygro)/sumabs
              scat_bc   = (scat_bc + scath2o*hyg_bc/sumhygro)/sumscat
              abs_bc    = (abs_bc + absh2o*hyg_bc/sumhygro)/sumabs
              scat_sslt = (scat_sslt + scath2o*hyg_sslt/sumhygro)/sumscat
              abs_sslt  = (abs_sslt + absh2o*hyg_sslt/sumhygro)/sumabs
            end if
            dopaer = tau_vis(icol, ilev, ibin)
            palb   = ssa_vis(icol, ilev, ibin)
            aodabsbc(icol) = aodabsbc(icol) + abs_bc*dopaer*(1.0_kind_phys - palb)
            aodc = (abs_sulf*(1.0_kind_phys - palb) + palb*scat_sulf)*dopaer
            sulfaod(icol) = sulfaod(icol) + aodc
            aodc = (abs_pom*(1.0_kind_phys - palb) + palb*scat_pom)*dopaer
            pomaod(icol) = pomaod(icol) + aodc
            aodc = (abs_soa*(1.0_kind_phys - palb) + palb*scat_soa)*dopaer
            soaaod(icol) = soaaod(icol) + aodc
            aodc = (abs_bc*(1.0_kind_phys - palb) + palb*scat_bc)*dopaer
            bcaod(icol) = bcaod(icol) + aodc
            aodc = (abs_sslt*(1.0_kind_phys - palb) + palb*scat_sslt)*dopaer
            ssltaod(icol) = ssltaod(icol) + aodc
            aodc = (abs_dust*(1.0_kind_phys - palb) + palb*scat_dust)*dopaer
            dustaod(icol) = dustaod(icol) + aodc
          end if
        end do
        odsum = 0._kind_phys
        do ilev = top_lev, pver
          odsum = odsum + tau_vis(icol, ilev, ibin)
        end do
        odv_col_aod(icol, ibin) = odsum
      end do

      ssavis(icol) = ssavis_default
      aodvis(icol) = 0._kind_phys
      do ilev = 1, pver
        aodvis(icol) = aodvis(icol) + aer_tau(icol, ilev, idx_sw_diag)
      end do
      if (aodvis(icol) > aodvis_min) then
        odsum = 0._kind_phys
        do ilev = 1, pver
          odsum = odsum + aer_tau_w(icol, ilev, idx_sw_diag)
        end do
        ssavis(icol) = odsum / aodvis(icol)
      end if
    end do

  end subroutine bulk_aerosol_optics_run

  !=============================================================================

  ! Mass mixing ratio of a bulk aerosol: a constituent, or zero (source Z,
  ! index 0), as the host's rad_cnst_get_aer_mmr [kg kg-1]
  pure real(kind_phys) function aerosol_mmr(constituents, idx, icol, ilev)
    real(kind_phys), intent(in) :: constituents(:, :, :)
    integer,         intent(in) :: idx, icol, ilev
    if (idx > 0) then
      aerosol_mmr = constituents(icol, ilev, idx)
    else
      aerosol_mmr = 0._kind_phys
    end if
  end function aerosol_mmr

  ! Geometric radius of a volcanic aerosol (0 when the aerosol has none) [m]
  pure real(kind_phys) function geometric_radius(constituents, idx, icol, ilev)
    real(kind_phys), intent(in) :: constituents(:, :, :)
    integer,         intent(in) :: idx, icol, ilev
    if (idx > 0) then
      geometric_radius = constituents(icol, ilev, idx)
    else
      geometric_radius = 0._kind_phys
    end if
  end function geometric_radius

  ! Relative-humidity table position (hygroscopic_aerosol_optics constructor)
  pure subroutine rh_position(relh, numrh, krh, wrh)
    real(kind_phys), intent(in)  :: relh
    integer,         intent(in)  :: numrh
    integer,         intent(out) :: krh
    real(kind_phys), intent(out) :: wrh
    real(kind_phys) :: rhtrunc
    rhtrunc = min(relh, 1._kind_phys)
    krh     = min(floor(rhtrunc * numrh) + 1, numrh - 1)
    wrh     = rhtrunc * numrh - krh
  end subroutine rh_position

  ! Radius table position (volcrad_aerosol_optics constructor)
  pure subroutine radius_position(radius, nmu, ibin, log_radius, kmu, wmu)
    real(kind_phys), intent(in)  :: radius
    integer,         intent(in)  :: nmu, ibin
    real(kind_phys), intent(in)  :: log_radius(:, :)
    integer,         intent(out) :: kmu
    real(kind_phys), intent(out) :: wmu
    real(kind_phys) :: r_mu_min, r_mu_max, mutrunc, mu
    r_mu_max = log_radius(nmu, ibin)
    r_mu_min = log_radius(1, ibin)
    if (radius > 0._kind_phys) then
      mu = log(radius)
    else
      mu = 0._kind_phys
    end if
    mutrunc = max(min(mu, r_mu_max), r_mu_min)
    ! Real to integer by assignment (truncation), as in volcrad
    kmu = max(min(1 + (mutrunc - r_mu_min)/(r_mu_max - r_mu_min)*(nmu - 1), nmu - 1._kind_phys), &
              1._kind_phys)
    wmu = max(min((mutrunc - log_radius(kmu, ibin)) / &
                  (log_radius(kmu+1, ibin) - log_radius(kmu, ibin)), 1._kind_phys), 0._kind_phys)
  end subroutine radius_position

  ! SW optical depth, single-scattering albedo and asymmetry parameter of one
  ! bulk aerosol (aerosol_optics_sw_bin with the classes' sw_props)
  pure subroutine sw_optics(method, nsamp, ibin, iwav, sw_ext, sw_ssa, sw_asm, log_radius, &
                            relh, radius, mmr, mass, tau, ssa, asm)
    integer,         intent(in)  :: method, nsamp, ibin, iwav
    real(kind_phys), intent(in)  :: sw_ext(:, :, :), sw_ssa(:, :, :), sw_asm(:, :, :)
    real(kind_phys), intent(in)  :: log_radius(:, :)
    real(kind_phys), intent(in)  :: relh, radius, mmr, mass
    real(kind_phys), intent(out) :: tau, ssa, asm
    integer         :: krh, kmu
    real(kind_phys) :: wrh, wmu, pext, palb, pasm, scat, ascat
    select case (method)
    case (optics_hygroscopic)
      call rh_position(relh, nsamp, krh, wrh)
      pext = (1._kind_phys + wrh) * sw_ext(krh+1, iwav, ibin) &
                           - wrh  * sw_ext(krh,   iwav, ibin)
      palb = (1._kind_phys + wrh) * sw_ssa(krh+1, iwav, ibin) &
                           - wrh  * sw_ssa(krh,   iwav, ibin)
      pasm = (1._kind_phys + wrh) * sw_asm(krh+1, iwav, ibin) &
                           - wrh  * sw_asm(krh,   iwav, ibin)
      pext = pext * mmr
    case (optics_insoluble)
      pext = sw_ext(1, iwav, ibin) * mmr
      palb = sw_ssa(1, iwav, ibin)
      pasm = sw_asm(1, iwav, ibin)
    case (optics_volcanic)
      call radius_position(radius, nsamp, ibin, log_radius, kmu, wmu)
      pext  = ((1._kind_phys - wmu) * sw_ext(kmu,   iwav, ibin) + &
                              (wmu) * sw_ext(kmu+1, iwav, ibin))
      scat  = ((1._kind_phys - wmu) * sw_ssa(kmu,   iwav, ibin) + &
                              (wmu) * sw_ssa(kmu+1, iwav, ibin))
      ascat = ((1._kind_phys - wmu) * sw_asm(kmu,   iwav, ibin) + &
                              (wmu) * sw_asm(kmu+1, iwav, ibin))
      palb = scat / pext
      if (scat > 0._kind_phys) then
        pasm = ascat / scat
      else
        pasm = 0._kind_phys
      end if
      pext = pext * mmr
    case default
      pext = 0._kind_phys
      palb = 0._kind_phys
      pasm = 0._kind_phys
    end select
    tau = pext * mass
    ssa = palb
    asm = pasm
  end subroutine sw_optics

  ! LW absorption optical depth of one bulk aerosol (aerosol_optics_lw_bin)
  pure subroutine lw_optics(method, nsamp, ibin, iwav, lw_abs, log_radius, &
                            relh, radius, mmr, mass, tau_lw)
    integer,         intent(in)  :: method, nsamp, ibin, iwav
    real(kind_phys), intent(in)  :: lw_abs(:, :, :)
    real(kind_phys), intent(in)  :: log_radius(:, :)
    real(kind_phys), intent(in)  :: relh, radius, mmr, mass
    real(kind_phys), intent(out) :: tau_lw
    integer         :: krh, kmu
    real(kind_phys) :: wrh, wmu, pabs
    select case (method)
    case (optics_hygroscopic)
      call rh_position(relh, nsamp, krh, wrh)
      pabs = (1._kind_phys + wrh) * lw_abs(krh+1, iwav, ibin) &
                           - wrh  * lw_abs(krh,   iwav, ibin)
      pabs = pabs * mmr
    case (optics_insoluble)
      pabs = lw_abs(1, iwav, ibin) * mmr
    case (optics_volcanic)
      call radius_position(radius, nsamp, ibin, log_radius, kmu, wmu)
      pabs = ((1._kind_phys - wmu) * lw_abs(kmu,   iwav, ibin) + &
                             (wmu) * lw_abs(kmu+1, iwav, ibin))
      pabs = pabs * mmr
    case default
      pabs = 0._kind_phys
    end select
    tau_lw = pabs * mass
  end subroutine lw_optics

end module bulk_aerosol_optics
