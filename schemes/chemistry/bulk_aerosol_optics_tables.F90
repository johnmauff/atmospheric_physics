! Optics tables and per-aerosol properties of the bulk aerosols (BAM) in the
! climate list, flattened into plain arrays for bulk_aerosol_optics (which
! can then run without the host's aerosol objects, e.g. on a GPU).
!
! The three bulk optics types the BAM physprops files use share one table
! layout, (sample, band, aerosol):
!   hygroscopic      samples = the relative-humidity mesh (phys_prop nrh);
!                    tables = extinction, single-scattering albedo,
!                    asymmetry parameter (SW) and absorption (LW)
!   insoluble        one sample; same quantities as hygroscopic
!   volcanic_radius  samples = log geometric-radius mesh; tables =
!                    extinction, scattering, asymmetric scattering (SW)
!                    and absorption (LW)
! Unused samples are zero.
!
! Adapted from aerosol_optics.F90 (Haipeng Lin, NSF-NCAR/CGD/AMP, March 2026)
! and CAM-SIMA's src/aerosol modules (bulk_aerosol_properties_mod,
! aerosol_mmr_host, radiative_aerosol).
module bulk_aerosol_optics_tables
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: bulk_aerosol_optics_tables_init

  ! Optics types (bulk_aerosol_optics_method_index)
  integer, parameter, public :: optics_hygroscopic = 1
  integer, parameter, public :: optics_insoluble   = 2
  integer, parameter, public :: optics_volcanic    = 3

  ! Species classes of the per-species diagnostics (bulk_aerosol_species_class_index)
  integer, parameter, public :: class_dust      = 1
  integer, parameter, public :: class_black_c   = 2
  integer, parameter, public :: class_sulfate   = 3
  integer, parameter, public :: class_p_organic = 4
  integer, parameter, public :: class_s_organic = 5
  integer, parameter, public :: class_seasalt   = 6

contains

!> \section arg_table_bulk_aerosol_optics_tables_init Argument Table
!! \htmlinclude bulk_aerosol_optics_tables_init.html
  subroutine bulk_aerosol_optics_tables_init(nswbands, nlwbands, num_bulk_aer, nsamples, &
       optics_method, idx_mmr, idx_radius, species_class, nsamples_aer, &
       density, hygro, refindex_re, refindex_im, &
       sw_ext, sw_ssa, sw_asm, lw_abs, log_radius, &
       errmsg, errflg)
    use ccpp_scheme_utils,       only: ccpp_constituent_index
    use atmos_phys_string_utils, only: to_lower
    ! The aerosol lists and optics properties are host data read by the
    ! host (as in aerosol_optics_init); they should be passed through the
    ! CCPP in the future.
    use radiative_aerosol_definitions, only: bulk_aerosol_list, id_climate
    use radiative_aerosol,             only: rad_aer_get_props
    use phys_prop,                     only: ot_length

    integer,          intent(in)  :: nswbands         ! number of SW bands [count]
    integer,          intent(in)  :: nlwbands         ! number of LW bands [count]
    integer,          intent(in)  :: num_bulk_aer     ! number of bulk aerosols [count]
    integer,          intent(in)  :: nsamples         ! samples in the tables [count]
    integer,          intent(out) :: optics_method(:) ! optics type [index]
    integer,          intent(out) :: idx_mmr(:)       ! constituent index of the mass mixing ratio; 0: zero [index]
    integer,          intent(out) :: idx_radius(:)    ! constituent index of the geometric radius; 0: none [index]
    integer,          intent(out) :: species_class(:) ! diagnostic species class [index]
    integer,          intent(out) :: nsamples_aer(:)  ! samples used per aerosol [count]
    real(kind_phys),  intent(out) :: density(:)       ! density [kg m-3]
    real(kind_phys),  intent(out) :: hygro(:)         ! hygroscopicity [1]
    real(kind_phys),  intent(out) :: refindex_re(:,:) ! SW refractive index, real part (0: none) [1]
    real(kind_phys),  intent(out) :: refindex_im(:,:) ! SW refractive index, imaginary part (0: none) [1]
    real(kind_phys),  intent(out) :: sw_ext(:,:,:)    ! SW specific extinction [m2 kg-1]
    real(kind_phys),  intent(out) :: sw_ssa(:,:,:)    ! SW ssa [1], or specific scattering (volcanic) [m2 kg-1]
    real(kind_phys),  intent(out) :: sw_asm(:,:,:)    ! SW asymmetry [1], or asymmetric scattering (volcanic) [m2 kg-1]
    real(kind_phys),  intent(out) :: lw_abs(:,:,:)    ! LW specific absorption [m2 kg-1]
    real(kind_phys),  intent(out) :: log_radius(:,:)  ! log geometric radius samples (volcanic) [ln(m)]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer                  :: ibin, nsamp
    character(len=ot_length) :: opticstype
    character(len=20)        :: aername
    character(len=16)        :: radius_name
    real(kind_phys),    pointer :: tbl_sw_ext(:,:), tbl_sw_ssa(:,:), tbl_sw_asm(:,:)
    real(kind_phys),    pointer :: tbl_lw_abs(:,:)
    real(kind_phys),    pointer :: vec_sw_ext(:), vec_sw_ssa(:), vec_sw_asm(:), vec_lw_abs(:)
    real(kind_phys),    pointer :: mu(:)
    complex(kind_phys), pointer :: refindex_sw(:)

    errmsg = ''
    errflg = 0

    optics_method(:) = 0
    idx_mmr(:)       = 0
    idx_radius(:)    = 0
    species_class(:) = 0
    nsamples_aer(:)  = 0
    density(:)       = 0._kind_phys
    hygro(:)         = 0._kind_phys
    refindex_re(:,:) = 0._kind_phys
    refindex_im(:,:) = 0._kind_phys
    sw_ext(:,:,:)    = 0._kind_phys
    sw_ssa(:,:,:)    = 0._kind_phys
    sw_asm(:,:,:)    = 0._kind_phys
    lw_abs(:,:,:)    = 0._kind_phys
    log_radius(:,:)  = 0._kind_phys

    do ibin = 1, num_bulk_aer
      call rad_aer_get_props(id_climate, ibin, opticstype=opticstype, aername=aername, &
                             density_aer=density(ibin), hygro_aer=hygro(ibin))

      ! Mass mixing ratio: a constituent (sources A and N) or zero (Z), as
      ! in the host's rad_cnst_get_aer_mmr
      select case (bulk_aerosol_list(id_climate)%aer(ibin)%source)
      case ('A', 'N')
        idx_mmr(ibin) = bulk_aerosol_list(id_climate)%aer(ibin)%idx
      case ('Z')
        idx_mmr(ibin) = 0
      case default
        errflg = 1
        errmsg = 'bulk_aerosol_optics_tables_init: unrecognized source for '//trim(aername)
        return
      end select

      ! Diagnostic species class from the name, as bulk_aerosol_properties%get
      select case (to_lower(aername(:4)))
      case ('dust')
        species_class(ibin) = class_dust
      case ('sulf', 'volc')
        species_class(ibin) = class_sulfate
      case ('bcar', 'bcph')
        species_class(ibin) = class_black_c
      case ('ocar', 'ocph')
        species_class(ibin) = class_p_organic
      case ('sslt', 'seas', 'ssam', 'sscm')
        species_class(ibin) = class_seasalt
      case default
        errflg = 1
        errmsg = 'bulk_aerosol_optics_tables_init: aerosol name not recognized: '//trim(aername)
        return
      end select

      nullify(refindex_sw)
      call rad_aer_get_props(id_climate, ibin, refindex_aer_sw=refindex_sw)
      if (associated(refindex_sw)) then
        ! Zero where there is none: the diagnostics then add the same zeros
        refindex_re(:, ibin) = real(refindex_sw(:nswbands), kind_phys)
        refindex_im(:, ibin) = aimag(refindex_sw(:nswbands))
      end if

      select case (trim(opticstype))
      case ('hygroscopic')
        nullify(tbl_sw_ext, tbl_sw_ssa, tbl_sw_asm, tbl_lw_abs)
        call rad_aer_get_props(id_climate, ibin, sw_hygro_ext=tbl_sw_ext, &
                               sw_hygro_ssa=tbl_sw_ssa, sw_hygro_asm=tbl_sw_asm, &
                               lw_hygro_ext=tbl_lw_abs)
        nsamp = size(tbl_sw_ext, 1)
        optics_method(ibin) = optics_hygroscopic
        nsamples_aer(ibin)  = nsamp
        sw_ext(:nsamp, :, ibin) = tbl_sw_ext(:, :nswbands)
        sw_ssa(:nsamp, :, ibin) = tbl_sw_ssa(:, :nswbands)
        sw_asm(:nsamp, :, ibin) = tbl_sw_asm(:, :nswbands)
        lw_abs(:nsamp, :, ibin) = tbl_lw_abs(:, :nlwbands)

      case ('insoluble', 'nonhygro')
        nullify(vec_sw_ext, vec_sw_ssa, vec_sw_asm, vec_lw_abs)
        call rad_aer_get_props(id_climate, ibin, sw_nonhygro_ext=vec_sw_ext, &
                               sw_nonhygro_ssa=vec_sw_ssa, sw_nonhygro_asm=vec_sw_asm, &
                               lw_ext=vec_lw_abs)
        optics_method(ibin) = optics_insoluble
        nsamples_aer(ibin)  = 1
        sw_ext(1, :, ibin) = vec_sw_ext(:nswbands)
        sw_ssa(1, :, ibin) = vec_sw_ssa(:nswbands)
        sw_asm(1, :, ibin) = vec_sw_asm(:nswbands)
        lw_abs(1, :, ibin) = vec_lw_abs(:nlwbands)

      case ('volcanic_radius', 'volcanic_radius1', 'volcanic_radius2', &
            'volcanic_radius3', 'volcanic_radius5')
        nullify(tbl_sw_ext, tbl_sw_ssa, tbl_sw_asm, tbl_lw_abs, mu)
        call rad_aer_get_props(id_climate, ibin, r_sw_ext=tbl_sw_ext, &
                               r_sw_scat=tbl_sw_ssa, r_sw_ascat=tbl_sw_asm, &
                               r_lw_abs=tbl_lw_abs, mu=mu)
        nsamp = size(mu)
        optics_method(ibin) = optics_volcanic
        nsamples_aer(ibin)  = nsamp
        ! The files store (band, sample)
        sw_ext(:nsamp, :, ibin) = transpose(tbl_sw_ext(:nswbands, :))
        sw_ssa(:nsamp, :, ibin) = transpose(tbl_sw_ssa(:nswbands, :))
        sw_asm(:nsamp, :, ibin) = transpose(tbl_sw_asm(:nswbands, :))
        lw_abs(:nsamp, :, ibin) = transpose(tbl_lw_abs(:nlwbands, :))
        log_radius(:nsamp, ibin) = mu(:)

        ! The geometric radius comes from a constituent named after the type
        radius_name = 'VOLC_RAD_GEOM '
        if (len_trim(opticstype) > 15) then
          radius_name = trim(radius_name)//opticstype(16:16)
        end if
        call ccpp_constituent_index(trim(radius_name), idx_radius(ibin), errflg, errmsg)
        if (errflg /= 0) return
        if (idx_radius(ibin) < 1) then
          errflg = 1
          errmsg = trim(radius_name)//' constituent required for '// &
                   trim(opticstype)//' optics but not found'
          return
        end if

      case default
        errflg = 1
        errmsg = 'bulk_aerosol_optics_tables_init: optics type '//trim(opticstype)// &
                 ' not supported; use the aerosol_optics scheme'
        return
      end select

      if (nsamples_aer(ibin) > nsamples) then
        errflg = 1
        errmsg = 'bulk_aerosol_optics_tables_init: table larger than nsamples for '// &
                 trim(aername)
        return
      end if
    end do

  end subroutine bulk_aerosol_optics_tables_init

end module bulk_aerosol_optics_tables
