! Aerosol optics diagnostics scheme for CAM-SIMA.
! Outputs aggregate SW/LW optical property diagnostics and per-species
! AODs/burdens from the aerosol_optics CCPP scheme.
!
! Day-only fields use fillvalue for night columns;
! day+night (dn) fields carry data for all columns.
!
! Author: Haipeng Lin, NSF-NCAR/CGD/AMP, March 2026
module aerosol_optics_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: aerosol_optics_diagnostics_init
   public :: aerosol_optics_diagnostics_run

   ! ODV_ aerosol diagnostics for bulk aerosol
   integer :: num_odv = 0
   character(len=64), allocatable :: odv_names(:)

contains

   !> \section arg_table_aerosol_optics_diagnostics_init  Argument Table
   !! \htmlinclude aerosol_optics_diagnostics_init.html
   subroutine aerosol_optics_diagnostics_init(errmsg, errflg)
      use cam_history,                   only: history_add_field
      use cam_history_support,           only: horiz_only
      use radiative_aerosol_definitions, only: bulk_aerosol_list, id_climate

      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      integer :: i

      errmsg = ''
      errflg = 0

      !-----------------------------------------------------------------
      ! SW column diagnostics (horiz_only, avg)
      !-----------------------------------------------------------------
      call history_add_field('AODVIS',      'Aerosol optical depth at 550 nm',                    horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODVISdn',    'Aerosol optical depth at 550 nm, day and night',      horiz_only, 'avg', '1')
      call history_add_field('AODUV',       'Aerosol optical depth in UV band',                    horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODUVdn',     'Aerosol optical depth in UV band, day and night',     horiz_only, 'avg', '1')
      call history_add_field('AODNIR',      'Aerosol optical depth in near-IR band',               horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODNIRdn',    'Aerosol optical depth in near-IR band, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODABS',      'Aerosol absorption optical depth at 550 nm',          horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODABSdn',    'Aerosol absorption optical depth at 550 nm, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODTOT',      'Total aerosol optical depth (all SW bands)',          horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODTOTdn',    'Total aerosol optical depth (all SW bands), day and night', horiz_only, 'avg', '1')
      call history_add_field('SSAVIS',      'Aerosol single scatter albedo at 550 nm',             horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('SSAVISdn',    'Aerosol single scatter albedo at 550 nm, day and night', horiz_only, 'avg', '1')

      ! Stratospheric column AODs
      call history_add_field('AODVISst',    'Stratospheric aerosol optical depth at 550 nm',       horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODVISstdn',  'Stratospheric aerosol optical depth at 550 nm, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODUVst',     'Stratospheric aerosol optical depth in UV band',      horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODUVstdn',   'Stratospheric aerosol optical depth in UV band, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODNIRst',    'Stratospheric aerosol optical depth in near-IR band', horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODNIRstdn',  'Stratospheric aerosol optical depth in near-IR band, day and night', horiz_only, 'avg', '1')

      !-----------------------------------------------------------------
      ! SW per-level diagnostics (lev, avg)
      !-----------------------------------------------------------------
      call history_add_field('EXTINCT',      'Aerosol extinction at 550 nm',                       'lev', 'avg', 'm-1', flag_xyfill=.true.)
      call history_add_field('EXTINCTdn',    'Aerosol extinction at 550 nm, day and night',         'lev', 'avg', 'm-1')
      call history_add_field('EXTINCTUV',    'Aerosol extinction in UV band',                       'lev', 'avg', 'm-1', flag_xyfill=.true.)
      call history_add_field('EXTINCTUVdn',  'Aerosol extinction in UV band, day and night',        'lev', 'avg', 'm-1')
      call history_add_field('EXTINCTNIR',   'Aerosol extinction in near-IR band',                  'lev', 'avg', 'm-1', flag_xyfill=.true.)
      call history_add_field('EXTINCTNIRdn', 'Aerosol extinction in near-IR band, day and night',   'lev', 'avg', 'm-1')
      call history_add_field('ABSORB',       'Aerosol absorption at 550 nm',                        'lev', 'avg', 'm-1', flag_xyfill=.true.)
      call history_add_field('ABSORBdn',     'Aerosol absorption at 550 nm, day and night',         'lev', 'avg', 'm-1')

      !-----------------------------------------------------------------
      ! Per-species column diagnostics (horiz_only, avg)
      !-----------------------------------------------------------------
      call history_add_field('AODDUST',        'Dust aerosol optical depth',          horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODDUSTdn',      'Dust aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODSO4',         'Sulfate aerosol optical depth',        horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODSO4dn',       'Sulfate aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODBC',          'Black carbon aerosol optical depth',   horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODBCdn',        'Black carbon aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODPOM',         'POM aerosol optical depth',            horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODPOMdn',       'POM aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODSOA',         'SOA aerosol optical depth',            horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODSOAdn',       'SOA aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODSS',          'Seasalt aerosol optical depth',        horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODSSdn',        'Seasalt aerosol optical depth, day and night', horiz_only, 'avg', '1')
      call history_add_field('AODABSBC',       'Black carbon absorption optical depth', horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODABSBCdn',     'Black carbon absorption optical depth, day and night', horiz_only, 'avg', '1')

      call history_add_field('BURDENDUST',     'Dust column burden',                   horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENDUSTdn',   'Dust column burden, day and night',    horiz_only, 'avg', 'kg m-2')
      call history_add_field('BURDENSO4',      'Sulfate column burden',                horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENSO4dn',    'Sulfate column burden, day and night',  horiz_only, 'avg', 'kg m-2')
      call history_add_field('BURDENBC',       'Black carbon column burden',           horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENBCdn',     'Black carbon column burden, day and night', horiz_only, 'avg', 'kg m-2')
      call history_add_field('BURDENPOM',      'POM column burden',                    horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENPOMdn',    'POM column burden, day and night',      horiz_only, 'avg', 'kg m-2')
      call history_add_field('BURDENSOA',      'SOA column burden',                    horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENSOAdn',    'SOA column burden, day and night',      horiz_only, 'avg', 'kg m-2')
      call history_add_field('BURDENSEASALT',  'Seasalt column burden',                horiz_only, 'avg', 'kg m-2', flag_xyfill=.true.)
      call history_add_field('BURDENSEASALTdn','Seasalt column burden, day and night',  horiz_only, 'avg', 'kg m-2')

      !-----------------------------------------------------------------
      ! LW per-level diagnostics (lev, avg)
      !-----------------------------------------------------------------
      call history_add_field('TOTABSLW',  'Total LW aerosol absorption optical depth', 'lev', 'avg', '1')
      call history_add_field('AODABSLW',  'LW aerosol absorption optical depth at diagnostic band', 'lev', 'avg', '1')

      !-----------------------------------------------------------------
      ! Additional CAM-standard column diagnostics (aliases)
      !-----------------------------------------------------------------
      call history_add_field('AEROD_v',  'Total aerosol optical depth in visible band', horiz_only, 'avg', '1', flag_xyfill=.true.)
      call history_add_field('AODvstrt', 'Stratospheric aerosol optical depth in visible band', horiz_only, 'avg', '1', flag_xyfill=.true.)

      !-----------------------------------------------------------------
      ! Per-aerosol visible OD diagnostics (ODV_<name>) for bulk aerosols
      !-----------------------------------------------------------------
      num_odv = bulk_aerosol_list(id_climate)%numaerosols
      if (num_odv > 0) then
        allocate(odv_names(num_odv), stat=errflg, errmsg=errmsg)
        if (errflg /= 0) then
          errmsg = 'aerosol_optics_diagnostics_init: ' // trim(errmsg)
          return
        end if
        do i = 1, num_odv
          odv_names(i) = 'ODV_' // trim(bulk_aerosol_list(id_climate)%aer(i)%camname)
          call history_add_field(trim(odv_names(i)), &
               trim(bulk_aerosol_list(id_climate)%aer(i)%camname)//' optical depth in visible band', &
               horiz_only, 'avg', '1', flag_xyfill=.true.)
        end do
      end if

   end subroutine aerosol_optics_diagnostics_init

   !> \section arg_table_aerosol_optics_diagnostics_run  Argument Table
   !! \htmlinclude aerosol_optics_diagnostics_run.html
   subroutine aerosol_optics_diagnostics_run( &
      ncol, pver, &
      N_DIAG, active_calls, &
      constituents, &
      aer_tau, aer_tau_w, aer_lw_abs, &
      idx_sw_diag, idx_uv_diag, idx_nir_diag, idx_lw_diag, &
      pdeldry, pmid, t, rga, rair, &
      nnite, idxnite, troplev, &
      dustaod, sulfaod, bcaod, pomaod, soaaod, ssltaod, aodabsbc, &
      burdendust, burdenso4, burdenbc, burdenpom, burdensoa, burdenseasalt, &
      ssavis, aodvis, &
      num_bulk_aer, odv_col_aod, &
      errmsg, errflg)

      use cam_history,         only: history_out_field
      use cam_history_support, only: fillvalue

      ! host-dependent diagnostic output subroutine for radiatively active aerosol:
      use aerosol_mmr_host,    only: rad_aer_diag_out

      !-----------------------------------------------------------------
      ! Input arguments
      !-----------------------------------------------------------------
      integer,         intent(in) :: ncol
      integer,         intent(in) :: pver
      integer,         intent(in) :: N_DIAG                ! max number of diagnostic lists
      logical,         intent(in) :: active_calls(0:)      ! flag for active diagnostic list calls
      ! n.b.: the lower bound of the active_calls array has to be set here.
      ! in radiative_aerosol_definitions, active_calls is 0-indexed because
      ! 0 is the climate list index.
      ! The Fortran 90 standard (and higher) states assumed-shape arrays
      ! lower bound default to 1 unless specified (ISO/IEC 1539 5.1.2.4.2)

      real(kind_phys), intent(in) :: constituents(:, :, :)

      real(kind_phys), intent(in) :: aer_tau(:,:,:)        ! SW extinction OD (ncol,pver,nswbands)
      real(kind_phys), intent(in) :: aer_tau_w(:,:,:)      ! SW ssa*tau (ncol,pver,nswbands)
      real(kind_phys), intent(in) :: aer_lw_abs(:,:,:)     ! LW absorption OD (ncol,pver,nlwbands)
      integer,         intent(in) :: idx_sw_diag           ! Index of SW diagnostic (vis) band
      integer,         intent(in) :: idx_uv_diag           ! Index of UV band
      integer,         intent(in) :: idx_nir_diag          ! Index of near-IR band
      integer,         intent(in) :: idx_lw_diag           ! Index of LW diagnostic band
      real(kind_phys), intent(in) :: pdeldry(:,:)          ! Dry air pressure thickness (ncol,pver) [Pa]
      real(kind_phys), intent(in) :: pmid(:,:)             ! Air pressure (ncol,pver) [Pa]
      real(kind_phys), intent(in) :: t(:,:)                ! Air temperature (ncol,pver) [K]
      real(kind_phys), intent(in) :: rga                   ! 1/g [s2 m-1]
      real(kind_phys), intent(in) :: rair                  ! Gas constant of dry air [J kg-1 K-1]
      integer,         intent(in) :: nnite                 ! Number of night columns
      integer,         intent(in) :: idxnite(:)            ! Indices of night columns
      integer,         intent(in) :: troplev(:)            ! Tropopause level index (ncol)

      ! Per-species inputs from aerosol_optics
      real(kind_phys), intent(in) :: dustaod(:)
      real(kind_phys), intent(in) :: sulfaod(:)
      real(kind_phys), intent(in) :: bcaod(:)
      real(kind_phys), intent(in) :: pomaod(:)
      real(kind_phys), intent(in) :: soaaod(:)
      real(kind_phys), intent(in) :: ssltaod(:)
      real(kind_phys), intent(in) :: aodabsbc(:)
      real(kind_phys), intent(in) :: burdendust(:)
      real(kind_phys), intent(in) :: burdenso4(:)
      real(kind_phys), intent(in) :: burdenbc(:)
      real(kind_phys), intent(in) :: burdenpom(:)
      real(kind_phys), intent(in) :: burdensoa(:)
      real(kind_phys), intent(in) :: burdenseasalt(:)
      real(kind_phys), intent(in) :: ssavis(:)
      real(kind_phys), intent(in) :: aodvis(:)

      ! Per-constituent visible OD inputs from aerosol_optics
      integer,         intent(in) :: num_bulk_aer
      real(kind_phys), intent(in) :: odv_col_aod(:,:)      ! (ncol, num_bulk_aer)

      ! CCPP error handling
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------
      integer :: icol, i, iaer, k, ib
      integer :: ilist

      ! Derived column quantities
      real(kind_phys) :: aodvis_l(ncol)
      real(kind_phys) :: aoduv(ncol)
      real(kind_phys) :: aodnir(ncol)
      real(kind_phys) :: aodabs(ncol)
      real(kind_phys) :: aodtot(ncol)
      real(kind_phys) :: ssavis_l(ncol)

      ! Stratospheric column AODs
      real(kind_phys) :: aodvisst(ncol)
      real(kind_phys) :: aoduvst(ncol)
      real(kind_phys) :: aodnirst(ncol)

      ! Per-level quantities
      real(kind_phys) :: mass(ncol, pver)
      real(kind_phys) :: air_density(ncol, pver)
      real(kind_phys) :: extinct(ncol, pver)
      real(kind_phys) :: extinctuv(ncol, pver)
      real(kind_phys) :: extinctnir(ncol, pver)
      real(kind_phys) :: absorb(ncol, pver)

      ! LW diagnostics
      real(kind_phys) :: totabslw(ncol, pver)
      real(kind_phys) :: aodabslw(ncol, pver)

      ! Local copy for ODV_ output with fillvalue
      real(kind_phys) :: odv_tmp(ncol)

      ! Local copies of per-species fields (for fillvalue modification)
      real(kind_phys) :: dustaod_l(ncol), sulfaod_l(ncol), bcaod_l(ncol)
      real(kind_phys) :: pomaod_l(ncol), soaaod_l(ncol), ssltaod_l(ncol)
      real(kind_phys) :: aodabsbc_l(ncol)
      real(kind_phys) :: burdendust_l(ncol), burdenso4_l(ncol), burdenbc_l(ncol)
      real(kind_phys) :: burdenpom_l(ncol), burdensoa_l(ncol), burdenseasalt_l(ncol)

      errmsg = ''
      errflg = 0

      ! Runs on the device: history_out_field accumulates device data in
      ! place (gpu-history-accumulate). Same arithmetic as before; the column
      ! sums are sequential loops over levels (and bands) in the order of
      ! the former sum() intrinsics. rad_aer_diag_out (aerosol_mmr_host) also
      ! runs on the device, on the constituents and pdeldry mapped here.
      !$acc data present(aer_tau, aer_tau_w, aer_lw_abs, pmid, t, troplev,    &
      !$acc              idxnite, dustaod, sulfaod, bcaod, pomaod, soaaod,      &
      !$acc              ssltaod, aodabsbc, burdendust, burdenso4, burdenbc,    &
      !$acc              burdenpom, burdensoa, burdenseasalt, ssavis, aodvis,   &
      !$acc              odv_col_aod, pdeldry, constituents)                    &
      !$acc      create(aodvis_l, aoduv, aodnir, aodabs, aodtot, ssavis_l,      &
      !$acc             aodvisst, aoduvst, aodnirst, mass, air_density,         &
      !$acc             extinct, extinctuv, extinctnir, absorb, totabslw,       &
      !$acc             aodabslw, odv_tmp, dustaod_l, sulfaod_l, bcaod_l,       &
      !$acc             pomaod_l, soaaod_l, ssltaod_l, aodabsbc_l,              &
      !$acc             burdendust_l, burdenso4_l, burdenbc_l, burdenpom_l,     &
      !$acc             burdensoa_l, burdenseasalt_l)

      !-----------------------------------------------------------------
      ! Derived quantities, per-level extinction [m-1] and absorption
      ! [m-1], LW diagnostics
      !-----------------------------------------------------------------
      !$acc parallel loop collapse(2)
      do k = 1, pver
        do icol = 1, ncol
          mass(icol,k) = pdeldry(icol,k) * rga
          air_density(icol,k) = pmid(icol,k) / (rair * t(icol,k))
          extinct(icol,k)    = aer_tau(icol,k,idx_sw_diag) * air_density(icol,k) / mass(icol,k)
          extinctuv(icol,k)  = aer_tau(icol,k,idx_uv_diag) * air_density(icol,k) / mass(icol,k)
          extinctnir(icol,k) = aer_tau(icol,k,idx_nir_diag) * air_density(icol,k) / mass(icol,k)
          absorb(icol,k)     = (aer_tau(icol,k,idx_sw_diag) - aer_tau_w(icol,k,idx_sw_diag)) &
                               * air_density(icol,k) / mass(icol,k)
          totabslw(icol,k) = 0._kind_phys
          !$acc loop seq
          do ib = 1, size(aer_lw_abs, 3)
            totabslw(icol,k) = totabslw(icol,k) + aer_lw_abs(icol,k,ib)
          end do
          aodabslw(icol,k) = aer_lw_abs(icol,k,idx_lw_diag)
        end do
      end do

      !-----------------------------------------------------------------
      ! Column AODs (sum over levels), stratospheric column AODs (levels
      ! above the tropopause), copies of the per-species fields for
      ! fillvalue modification
      !-----------------------------------------------------------------
      !$acc parallel loop
      do icol = 1, ncol
        ! aodvis and ssavis are already computed in aerosol_optics_run
        aodvis_l(icol) = aodvis(icol)
        ssavis_l(icol) = ssavis(icol)

        aoduv(icol)  = 0._kind_phys
        aodnir(icol) = 0._kind_phys
        aodabs(icol) = 0._kind_phys
        !$acc loop seq
        do k = 1, size(aer_tau, 2)
          aoduv(icol)  = aoduv(icol)  + aer_tau(icol,k,idx_uv_diag)
          aodnir(icol) = aodnir(icol) + aer_tau(icol,k,idx_nir_diag)
          aodabs(icol) = aodabs(icol) + (aer_tau(icol,k,idx_sw_diag) - aer_tau_w(icol,k,idx_sw_diag))
        end do
        aodtot(icol) = 0._kind_phys
        !$acc loop seq
        do ib = 1, size(aer_tau, 3)
          !$acc loop seq
          do k = 1, size(aer_tau, 2)
            aodtot(icol) = aodtot(icol) + aer_tau(icol,k,ib)
          end do
        end do

        aodvisst(icol) = 0._kind_phys
        aoduvst(icol)  = 0._kind_phys
        aodnirst(icol) = 0._kind_phys
        !$acc loop seq
        do k = 1, troplev(icol)
          aodvisst(icol) = aodvisst(icol) + aer_tau(icol,k,idx_sw_diag)
          aoduvst(icol)  = aoduvst(icol)  + aer_tau(icol,k,idx_uv_diag)
          aodnirst(icol) = aodnirst(icol) + aer_tau(icol,k,idx_nir_diag)
        end do

        dustaod_l(icol)       = dustaod(icol)
        sulfaod_l(icol)       = sulfaod(icol)
        bcaod_l(icol)         = bcaod(icol)
        pomaod_l(icol)        = pomaod(icol)
        soaaod_l(icol)        = soaaod(icol)
        ssltaod_l(icol)       = ssltaod(icol)
        aodabsbc_l(icol)      = aodabsbc(icol)
        burdendust_l(icol)    = burdendust(icol)
        burdenso4_l(icol)     = burdenso4(icol)
        burdenbc_l(icol)      = burdenbc(icol)
        burdenpom_l(icol)     = burdenpom(icol)
        burdensoa_l(icol)     = burdensoa(icol)
        burdenseasalt_l(icol) = burdenseasalt(icol)
      end do

      !-----------------------------------------------------------------
      ! Output day+night (dn) variants first (all columns valid)
      !-----------------------------------------------------------------

      ! SW column
      call history_out_field('AODVISdn',   aodvis_l)
      call history_out_field('AODUVdn',    aoduv)
      call history_out_field('AODNIRdn',   aodnir)
      call history_out_field('AODABSdn',   aodabs)
      call history_out_field('AODTOTdn',   aodtot)
      call history_out_field('SSAVISdn',   ssavis_l)

      ! Stratospheric
      call history_out_field('AODVISstdn', aodvisst)
      call history_out_field('AODUVstdn',  aoduvst)
      call history_out_field('AODNIRstdn', aodnirst)

      ! SW per-level
      call history_out_field('EXTINCTdn',    extinct)
      call history_out_field('EXTINCTUVdn',  extinctuv)
      call history_out_field('EXTINCTNIRdn', extinctnir)
      call history_out_field('ABSORBdn',     absorb)

      ! Per-species column (dn)
      call history_out_field('AODDUSTdn',       dustaod_l)
      call history_out_field('AODSO4dn',        sulfaod_l)
      call history_out_field('AODBCdn',         bcaod_l)
      call history_out_field('AODPOMdn',        pomaod_l)
      call history_out_field('AODSOAdn',        soaaod_l)
      call history_out_field('AODSSdn',         ssltaod_l)
      call history_out_field('AODABSBCdn',      aodabsbc_l)

      call history_out_field('BURDENDUSTdn',    burdendust_l)
      call history_out_field('BURDENSO4dn',     burdenso4_l)
      call history_out_field('BURDENBCdn',      burdenbc_l)
      call history_out_field('BURDENPOMdn',     burdenpom_l)
      call history_out_field('BURDENSOAdn',     burdensoa_l)
      call history_out_field('BURDENSEASALTdn', burdenseasalt_l)

      ! LW per-level (no day/night distinction)
      call history_out_field('TOTABSLW',  totabslw)
      call history_out_field('AODABSLW',  aodabslw)

      !-----------------------------------------------------------------
      ! Apply fillvalue to night columns for day-only variants
      ! (idxnite(1:nnite) are distinct columns)
      !-----------------------------------------------------------------
      !$acc parallel loop
      do i = 1, nnite
        aodvis_l(idxnite(i))  = fillvalue
        aoduv(idxnite(i))     = fillvalue
        aodnir(idxnite(i))    = fillvalue
        aodabs(idxnite(i))    = fillvalue
        aodtot(idxnite(i))    = fillvalue
        ssavis_l(idxnite(i))  = fillvalue
        aodvisst(idxnite(i))  = fillvalue
        aoduvst(idxnite(i))   = fillvalue
        aodnirst(idxnite(i))  = fillvalue

        dustaod_l(idxnite(i))       = fillvalue
        sulfaod_l(idxnite(i))       = fillvalue
        bcaod_l(idxnite(i))         = fillvalue
        pomaod_l(idxnite(i))        = fillvalue
        soaaod_l(idxnite(i))        = fillvalue
        ssltaod_l(idxnite(i))       = fillvalue
        aodabsbc_l(idxnite(i))      = fillvalue

        burdendust_l(idxnite(i))    = fillvalue
        burdenso4_l(idxnite(i))     = fillvalue
        burdenbc_l(idxnite(i))      = fillvalue
        burdenpom_l(idxnite(i))     = fillvalue
        burdensoa_l(idxnite(i))     = fillvalue
        burdenseasalt_l(idxnite(i)) = fillvalue
      end do
      !$acc parallel loop collapse(2)
      do k = 1, pver
        do i = 1, nnite
          extinct(idxnite(i), k)    = fillvalue
          extinctuv(idxnite(i), k)  = fillvalue
          extinctnir(idxnite(i), k) = fillvalue
          absorb(idxnite(i), k)     = fillvalue
        end do
      end do

      !-----------------------------------------------------------------
      ! Output day-only variants (night columns have fillvalue)
      !-----------------------------------------------------------------

      ! SW column
      call history_out_field('AODVIS',   aodvis_l)
      call history_out_field('AODUV',    aoduv)
      call history_out_field('AODNIR',   aodnir)
      call history_out_field('AODABS',   aodabs)
      call history_out_field('AODTOT',   aodtot)
      call history_out_field('SSAVIS',   ssavis_l)

      ! Stratospheric
      call history_out_field('AODVISst', aodvisst)
      call history_out_field('AODUVst',  aoduvst)
      call history_out_field('AODNIRst', aodnirst)

      ! SW per-level
      call history_out_field('EXTINCT',    extinct)
      call history_out_field('EXTINCTUV',  extinctuv)
      call history_out_field('EXTINCTNIR', extinctnir)
      call history_out_field('ABSORB',     absorb)

      ! Per-species column
      call history_out_field('AODDUST',       dustaod_l)
      call history_out_field('AODSO4',        sulfaod_l)
      call history_out_field('AODBC',         bcaod_l)
      call history_out_field('AODPOM',        pomaod_l)
      call history_out_field('AODSOA',        soaaod_l)
      call history_out_field('AODSS',         ssltaod_l)
      call history_out_field('AODABSBC',      aodabsbc_l)

      call history_out_field('BURDENDUST',    burdendust_l)
      call history_out_field('BURDENSO4',     burdenso4_l)
      call history_out_field('BURDENBC',      burdenbc_l)
      call history_out_field('BURDENPOM',     burdenpom_l)
      call history_out_field('BURDENSOA',     burdensoa_l)
      call history_out_field('BURDENSEASALT', burdenseasalt_l)

      ! AEROD_v and AODvstrt (aliases for AODVIS/AODVISst, day-only)
      call history_out_field('AEROD_v',  aodvis_l)
      call history_out_field('AODvstrt', aodvisst)

      ! Per-aerosol visible OD (ODV_<name>), day-only with fillvalue
      do iaer = 1, min(num_bulk_aer, num_odv)
        !$acc parallel loop
        do icol = 1, ncol
          odv_tmp(icol) = odv_col_aod(icol, iaer)
        end do
        !$acc parallel loop
        do i = 1, nnite
          odv_tmp(idxnite(i)) = fillvalue
        end do
        call history_out_field(trim(odv_names(iaer)), odv_tmp)
      end do

      ! Output diagnostics for active climate and diagnostic lists
      do ilist = 0, N_DIAG
        if (active_calls(ilist)) then
          call rad_aer_diag_out(ilist, constituents, pdeldry, ncol)
        end if
      end do

      !$acc end data

   end subroutine aerosol_optics_diagnostics_run

end module aerosol_optics_diagnostics
