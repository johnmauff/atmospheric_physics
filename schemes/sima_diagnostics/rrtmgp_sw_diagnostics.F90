module rrtmgp_sw_diagnostics
   use ccpp_kinds, only:  kind_phys

   implicit none
   private

   public :: rrtmgp_sw_diagnostics_init ! init routine
   public :: rrtmgp_sw_diagnostics_run  ! main routine

   integer, parameter :: N_DIAG=10
   character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ',&
                                      '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

CONTAINS

   !> \section arg_table_rrtmgp_sw_diagnostics_init  Argument Table
   !! \htmlinclude rrtmgp_sw_diagnostics_init.html
   subroutine rrtmgp_sw_diagnostics_init(num_diag_subcycles, active_calls, errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      integer,             intent(in) :: num_diag_subcycles  ! Number of diagnostic subcycles
      logical,             intent(in) :: active_calls(:)     ! Logical array of flags for whether a specified subcycle is active

      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables:
      integer :: icall

      errmsg = ''
      errflg = 0

      ! Loop over number of diagnostics subcycles
      !  and add the relevant fields for each cycle if it's active
      do icall = 1, num_diag_subcycles
         if (active_calls(icall)) then
            call history_add_field('SOLIN'//diag(icall-1),   'Solar isolation',                                   horiz_only, 'avg',  'W m-2')
            call history_add_field('QRS'//diag(icall-1),     'Solar heating rate',                                'lev',      'avg',  'K s-1')
            call history_add_field('QRSC'//diag(icall-1),    'Clearsky solar heating rate',                       'lev',      'avg',  'K s-1')
            call history_add_field('FSNT'//diag(icall-1),    'Net solar flux at top of model',                    horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNTC'//diag(icall-1),   'Clearky net solar flux at top of model',            horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNTOA'//diag(icall-1),  'Net solar flux at top of atmosphere',               horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNTOAC'//diag(icall-1), 'Clearsky net solar flux at top of atmosphere',      horiz_only, 'avg',  'W m-2')
            call history_add_field('SWCF'//diag(icall-1),    'Shortwave cloud forcing',                           horiz_only, 'avg',  'W m-2')
            call history_add_field('FSUTOA'//diag(icall-1),  'Upwelling solar flux at top of atmospehre',         horiz_only, 'avg',  'W m-2')
            call history_add_field('FSN200'//diag(icall-1),  'Net shortwave flux at 200 mb',                      horiz_only, 'avg',  'W m-2')
            call history_add_field('FSN200C'//diag(icall-1), 'Clearsky net shortwave flux at 200 mb',             horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNR'//diag(icall-1),    'Net solar flux at tropopause',                      horiz_only, 'avg',  'W m-2')
            call history_add_field('SOLL'//diag(icall-1),    'Solar downward near infrared direct to surface',    horiz_only, 'avg',  'W m-2')
            call history_add_field('SOLS'//diag(icall-1),    'Solar downward visible direct to surface',          horiz_only, 'avg',  'W m-2')
            call history_add_field('SOLLD'//diag(icall-1),   'Solar downward near infrared diffuse to surface',   horiz_only, 'avg',  'W m-2')
            call history_add_field('SOLSD'//diag(icall-1),   'Solar downward visible diffuse to surface',         horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNS'//diag(icall-1),    'Net solar flux at surface',                         horiz_only, 'avg',  'W m-2')
            call history_add_field('FSNSC'//diag(icall-1),   'Clearsky net solar flux at surface',                horiz_only, 'avg',  'W m-2')
            call history_add_field('FSDS'//diag(icall-1),    'Downwelling solar flux at surface',                 horiz_only, 'avg',  'W m-2')
            call history_add_field('FSDSC'//diag(icall-1),   'Clearsky downwelling solar flux at surface',        horiz_only, 'avg',  'W m-2')

            ! Fluxes on CAM grid
            call history_add_field('FUS'//diag(icall-1),     'Shortwave upward flux',                             'ilev',     'inst', 'W m-2')
            call history_add_field('FDS'//diag(icall-1),     'Shortwave downward flux',                           'ilev',     'inst', 'W m-2')
            call history_add_field('FUSC'//diag(icall-1),    'Shortwave clear-sky upward flux',                   'ilev',     'inst', 'W m-2')
            call history_add_field('FDSC'//diag(icall-1),    'Shortwave clear-sky downward flux',                 'ilev',     'inst', 'W m-2')
         end if
      end do

   end subroutine rrtmgp_sw_diagnostics_init

   !> \section arg_table_rrtmgp_sw_diagnostics_run  Argument Table
   !! \htmlinclude rrtmgp_sw_diagnostics_run.html
   subroutine rrtmgp_sw_diagnostics_run(num_diag_subcycles, icall, active_calls, fsw, fswc, rpdel, ncol, nday, idxday, &
                  nlay, pver, pverp, pint, gravit, cpair, p_trop, fns, fcns, qrs, qrsc, fsnt, fsns, sols, soll, solsd, &
                  solld, ktopcam, ktoprad, write_output, dosw, errmsg, errflg)

      use cam_history,        only: history_out_field
      use ccpp_fluxes,        only: ty_fluxes_broadband_ccpp
      use ccpp_fluxes_byband, only: ty_fluxes_byband_ccpp
      use atmos_phys_interp_utils, only: vertinterp_column
      !------------------------------------------------
      !   Input / output parameters
      !------------------------------------------------
      integer,                        intent(in) :: num_diag_subcycles  ! Number of diagnostics subcycles
      integer,                        intent(in) :: icall               ! Current diagnostic subcycle
      integer,                        intent(in) :: ncol                ! Number of horizontal points
      integer,                        intent(in) :: nday                ! Number of daytime points
      integer,                        intent(in) :: idxday(:)           ! Daytime points
      integer,                        intent(in) :: nlay                ! Number of vertical layers used in radiation calcluation
      integer,                        intent(in) :: pver                ! Number of vertical layers
      integer,                        intent(in) :: pverp               ! Number of vertical layer interfaces
      integer,                        intent(in) :: ktopcam             ! Index in CAM arrays of top level (layer or interface) at which RRTMGP is active
      integer,                        intent(in) :: ktoprad             ! Index in RRTMGP array corresponding to top layer or interface of CAM arrays
      logical,                        intent(in) :: active_calls(:)     ! Logical array of flags for whether a specified subcycle is active
      logical,                        intent(in) :: write_output        ! Flag to write output for radiation
      logical,                        intent(in) :: dosw                ! Flag for whether shortwave radiation was calculated this timestep
      real(kind_phys),                intent(in) :: gravit              ! Standard gravitiational acceleration
      real(kind_phys),                intent(in) :: cpair
      real(kind_phys),                intent(in) :: pint(:,:)           ! Air pressure at layer interfaces [Pa]
      real(kind_phys),                intent(in) :: p_trop(:)           ! Tropopause air pressure [Pa]
      real(kind_phys),                intent(in) :: rpdel(:,:)          ! Reciprocal of layer thickness [Pa-1]
      real(kind_phys),                intent(in) :: fns(:,:)            ! Net shortwave all-sky flux [W m-2]
      real(kind_phys),                intent(in) :: fcns(:,:)           ! Net shortwave clear-sky flux [W m-2]
      real(kind_phys),                intent(in) :: qrs(:,:)            ! Heating rate (all-sky, shortwave) [J kg-1 s-1]
      real(kind_phys),                intent(in) :: qrsc(:,:)           ! Heating rate (clear-sky, shortwave) [J kg-1 s-1]
      real(kind_phys),                intent(in) :: fsnt(:)             ! Shortwave outgoing flux at model top [W m-2]
      real(kind_phys),                intent(in) :: fsns(:)             ! Shortwave upward flux at surface [W m-2]
      real(kind_phys),                intent(in) :: sols(:)
      real(kind_phys),                intent(in) :: soll(:)
      real(kind_phys),                intent(in) :: solsd(:)
      real(kind_phys),                intent(in) :: solld(:)
      type(ty_fluxes_byband_ccpp),    intent(in) :: fsw                 ! Shortwave all-sky flux object
      type(ty_fluxes_broadband_ccpp), intent(in) :: fswc                ! Shortwave clear-sky flux object
      
      ! CCPP error handling variables
      character(len=512), intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      integer :: diag_index, i, k, nbad
      logical :: ok200, ok200c, okr
      real(kind_phys) :: swcf(ncol)
      real(kind_phys) :: qrs_cp(size(qrs,1), size(qrs,2))   ! qrs/cpair
      real(kind_phys) :: qrsc_cp(size(qrsc,1), size(qrsc,2)) ! qrsc/cpair
      real(kind_phys) :: fsds(ncol)
      real(kind_phys) :: solin(ncol)
      real(kind_phys) :: fsntoa(ncol)
      real(kind_phys) :: fsntoac(ncol)
      real(kind_phys) :: fsutoa(ncol)
      real(kind_phys) :: fsdsc(ncol)
      real(kind_phys) :: flux_sw_up(ncol,pverp)
      real(kind_phys) :: flux_sw_dn(ncol,pverp)
      real(kind_phys) :: flux_sw_clr_up(ncol,pverp)
      real(kind_phys) :: flux_sw_clr_dn(ncol,pverp)
      real(kind_phys) :: fsn200(ncol)
      real(kind_phys) :: fsn200c(ncol)
      real(kind_phys) :: fsnr(ncol)

      errmsg = ''
      errflg = 0

      ! Diagnostic indices are reversed
      diag_index = num_diag_subcycles - icall

      ! Don't do anything if this subcycle is inactive,
      !  we're not configured to write radiation output, or
      !  radiation was not run this timestep.
      !
      ! The flux objects are zeroed on non-radiation timesteps so they should not be
      ! written to accumulated history fields on non-radiation timesteps.
      if ((.not. dosw) .or. (.not. active_calls(diag_index+1)) .or. (.not. write_output)) then
         return
      end if

      ! Everything below runs on the device: the inputs are device data, the
      ! flux object arrays too (rrtmgp_pre; only when nday > 0), and every
      ! field goes to the history buffer as a device array (a local created
      ! here or a contiguous section of an input). The flux arrays are only
      ! read through plain-array dummies (scatter_daylit), not as components
      ! in a kernel (nvfortran "partially present").
      !$acc data present(idxday, pint, p_trop, fns, fcns, qrs, qrsc, fsnt, fsns, &
      !$acc              sols, soll, solsd, solld) &
      !$acc      create(solin, fsntoa, fsntoac, fsutoa, fsdsc, fsds, flux_sw_up, &
      !$acc             flux_sw_dn, flux_sw_clr_up, flux_sw_clr_dn, swcf, fsn200, &
      !$acc             fsn200c, fsnr, qrs_cp, qrsc_cp)

      !$acc parallel loop
      do i = 1, ncol
         solin(i) = 0.0_kind_phys
         fsntoa(i) = 0.0_kind_phys
         fsntoac(i) = 0.0_kind_phys
         fsutoa(i) = 0.0_kind_phys
         fsdsc(i) = 0.0_kind_phys
         fsds(i) = 0.0_kind_phys
      end do
      !$acc parallel loop collapse(2)
      do k = 1, pverp
         do i = 1, ncol
            flux_sw_up(i,k) = 0.0_kind_phys
            flux_sw_dn(i,k) = 0.0_kind_phys
            flux_sw_clr_up(i,k) = 0.0_kind_phys
            flux_sw_clr_dn(i,k) = 0.0_kind_phys
         end do
      end do

      ! Daylit columns from the flux objects (nday columns) to the CAM grid.
      ! With no daylit columns the flux arrays have no device copies.
      if (nday > 0) then
         call scatter_daylit(nday, idxday, nlay, ktopcam, ktoprad,          &
              fsw%fluxes%flux_up, fsw%fluxes%flux_dn, fsw%fluxes%flux_net,  &
              fswc%fluxes%flux_up, fswc%fluxes%flux_dn, fswc%fluxes%flux_net, &
              fsds, solin, fsntoa, fsntoac, fsutoa, fsdsc,                  &
              flux_sw_up, flux_sw_dn, flux_sw_clr_up, flux_sw_clr_dn)
      end if

      !$acc parallel loop
      do i = 1, ncol
         swcf(i) = fsntoa(i) - fsntoac(i)
      end do
      !$acc parallel loop collapse(2)
      do k = 1, size(qrs, 2)
         do i = 1, size(qrs, 1)
            qrs_cp(i,k) = qrs(i,k)/cpair
            qrsc_cp(i,k) = qrsc(i,k)/cpair
         end do
      end do

      ! FSN200 and FSNR (vertinterp, default method, one column at a time)
      nbad = 0
      !$acc parallel loop reduction(+:nbad)
      do i = 1, ncol
         call vertinterp_column(i, pverp, pint, 20000._kind_phys, fns,  fsn200(i),  ok200)
         call vertinterp_column(i, pverp, pint, 20000._kind_phys, fcns, fsn200c(i), ok200c)
         call vertinterp_column(i, pverp, pint, p_trop(i),        fns,  fsnr(i),    okr)
         if (.not. (ok200 .and. ok200c .and. okr)) nbad = nbad + 1
      end do

      if (nbad == 0) then
         call history_out_field('SOLIN'//diag(diag_index),   solin)
         call history_out_field('QRS'//diag(diag_index),     qrs_cp)
         call history_out_field('QRSC'//diag(diag_index),    qrsc_cp)

         call history_out_field('FSNT'//diag(diag_index),    fsnt)
         call history_out_field('FSNTC'//diag(diag_index),   fcns(:, ktopcam)) ! net sw clearsky flux at top
         call history_out_field('FSNTOA'//diag(diag_index),  fsntoa)
         call history_out_field('FSNTOAC'//diag(diag_index), fsntoac)

         call history_out_field('SWCF'//diag(diag_index),    swcf)

         call history_out_field('FSUTOA'//diag(diag_index),  fsutoa)

         call history_out_field('FSN200'//diag(diag_index),  fsn200)
         call history_out_field('FSN200C'//diag(diag_index), fsn200c)

         call history_out_field('FSNR'//diag(diag_index),    fsnr)

         call history_out_field('SOLS'//diag(diag_index),    sols)
         call history_out_field('SOLL'//diag(diag_index),    soll)
         call history_out_field('SOLSD'//diag(diag_index),   solsd)
         call history_out_field('SOLLD'//diag(diag_index),   solld)

         call history_out_field('FSNS'//diag(diag_index),    fsns)
         call history_out_field('FSNSC'//diag(diag_index),   fcns(:, pverp))   ! net sw clearsky flux at surface

         call history_out_field('FSDS'//diag(diag_index),    fsds)
         call history_out_field('FSDSC'//diag(diag_index),   fsdsc)

         ! Fluxes on the CAM grid
         call history_out_field('FDS'//diag(diag_index),     flux_sw_dn)
         call history_out_field('FDSC'//diag(diag_index),    flux_sw_clr_dn)
         call history_out_field('FUS'//diag(diag_index),     flux_sw_up)
         call history_out_field('FUSC'//diag(diag_index),    flux_sw_clr_up)
      end if

      !$acc end data

      if (nbad > 0) then
         errmsg = 'VERTINTERP: output pressure not bracketed by the input levels'
         errflg = 1
      end if

   end subroutine rrtmgp_sw_diagnostics_run

   !=========================================================================================

   subroutine scatter_daylit(nday, idxday, nlay, ktopcam, ktoprad,      &
                             up, dn, net, clr_up, clr_dn, clr_net,      &
                             fsds, solin, fsntoa, fsntoac, fsutoa, fsdsc, &
                             flux_up, flux_dn, flux_clr_up, flux_clr_dn)
      ! Copy the daylit columns (1:nday) of the shortwave flux arrays to
      ! their CAM columns idxday(1:nday), on the device (all present). The
      ! idxday values are distinct, so the columns are independent.
      integer,         intent(in)    :: nday, nlay, ktopcam, ktoprad
      integer,         intent(in)    :: idxday(:)
      real(kind_phys), intent(in)    :: up(:,:), dn(:,:), net(:,:)            ! all-sky (nday, nlay+1)
      real(kind_phys), intent(in)    :: clr_up(:,:), clr_dn(:,:), clr_net(:,:) ! clear-sky
      real(kind_phys), intent(inout) :: fsds(:), solin(:), fsntoa(:), fsntoac(:), fsutoa(:), fsdsc(:)
      real(kind_phys), intent(inout) :: flux_up(:,:), flux_dn(:,:), flux_clr_up(:,:), flux_clr_dn(:,:)

      integer :: idx, k

      !$acc data present(idxday, up, dn, net, clr_up, clr_dn, clr_net, fsds, solin, &
      !$acc              fsntoa, fsntoac, fsutoa, fsdsc, flux_up, flux_dn, flux_clr_up, &
      !$acc              flux_clr_dn)
      !$acc parallel loop
      do idx = 1, nday
         fsds(idxday(idx)) = dn(idx, nlay+1)
         solin(idxday(idx)) = clr_dn(idx, 1)
         fsntoa(idxday(idx)) = net(idx, 1)
         fsntoac(idxday(idx)) = clr_net(idx, 1)
         fsutoa(idxday(idx)) = up(idx, 1)
         fsdsc(idxday(idx))  = clr_dn(idx, nlay+1)
      end do
      !$acc parallel loop collapse(2)
      do k = ktopcam, size(flux_up, 2)
         do idx = 1, nday
            flux_up(idxday(idx), k) = up(idx, k - ktopcam + ktoprad)
            flux_dn(idxday(idx), k) = dn(idx, k - ktopcam + ktoprad)
            flux_clr_up(idxday(idx), k) = clr_up(idx, k - ktopcam + ktoprad)
            flux_clr_dn(idxday(idx), k) = clr_dn(idx, k - ktopcam + ktoprad)
         end do
      end do
      !$acc end data

   end subroutine scatter_daylit

end module rrtmgp_sw_diagnostics
