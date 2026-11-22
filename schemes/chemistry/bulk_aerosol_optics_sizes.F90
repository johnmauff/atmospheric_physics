! Sizes of the bulk aerosol (BAM) optics tables built by
! bulk_aerosol_optics_tables and used by bulk_aerosol_optics.
! A separate scheme so that the tables, dimensioned by these sizes, can be
! allocated before bulk_aerosol_optics_tables fills them. The tables hold at
! least one (zero) aerosol: with no bulk aerosols (e.g. aquaplanet) a
! zero-size array kept on the GPU for the whole run can collide with later
! device data.
!
! Adapted from aerosol_optics.F90 (Haipeng Lin, NSF-NCAR/CGD/AMP, March 2026).
module bulk_aerosol_optics_sizes
  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: bulk_aerosol_optics_sizes_init

contains

!> \section arg_table_bulk_aerosol_optics_sizes_init Argument Table
!! \htmlinclude bulk_aerosol_optics_sizes_init.html
  subroutine bulk_aerosol_optics_sizes_init(num_bulk_aer, nsamples, naer_tables, errmsg, errflg)
    ! The aerosol lists and optics properties are host data read by the
    ! host (as in aerosol_optics_init); they should be passed through the
    ! CCPP in the future.
    use radiative_aerosol_definitions, only: bulk_aerosol_list, id_climate
    use radiative_aerosol,             only: rad_aer_get_props
    use aerosol_instances_mod,         only: aerosol_instances_is_active
    use phys_prop,                     only: ot_length

    integer,          intent(out) :: num_bulk_aer ! number of bulk aerosols [count]
    integer,          intent(out) :: nsamples     ! samples in the optics tables [count]
    integer,          intent(out) :: naer_tables  ! aerosols in the optics tables, at least 1 [count]
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer                  :: ibin
    character(len=ot_length) :: opticstype
    real(kind_phys), pointer :: sw_hygro_ext(:, :)
    real(kind_phys), pointer :: mu(:)

    errmsg = ''
    errflg = 0

    if (aerosol_instances_is_active('modal') .or. &
        aerosol_instances_is_active('carma')) then
      errflg = 1
      errmsg = 'bulk_aerosol_optics_sizes_init: modal or sectional aerosols are '// &
               'active; use the aerosol_optics scheme'
      return
    end if

    num_bulk_aer = bulk_aerosol_list(id_climate)%numaerosols
    naer_tables  = max(num_bulk_aer, 1)

    ! One sample at least: the insoluble tables have one per band.
    nsamples = 1
    do ibin = 1, num_bulk_aer
      call rad_aer_get_props(id_climate, ibin, opticstype=opticstype)
      select case (trim(opticstype))
      case ('hygroscopic')
        nullify(sw_hygro_ext)
        call rad_aer_get_props(id_climate, ibin, sw_hygro_ext=sw_hygro_ext)
        nsamples = max(nsamples, size(sw_hygro_ext, 1))
      case ('insoluble', 'nonhygro')
        ! one sample per band
      case ('volcanic_radius', 'volcanic_radius1', 'volcanic_radius2', &
            'volcanic_radius3', 'volcanic_radius5')
        nullify(mu)
        call rad_aer_get_props(id_climate, ibin, mu=mu)
        nsamples = max(nsamples, size(mu))
      case default
        errflg = 1
        errmsg = 'bulk_aerosol_optics_sizes_init: optics type '//trim(opticstype)// &
                 ' not supported; use the aerosol_optics scheme'
        return
      end select
    end do

  end subroutine bulk_aerosol_optics_sizes_init

end module bulk_aerosol_optics_sizes
