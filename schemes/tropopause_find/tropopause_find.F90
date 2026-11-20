! This module is used to diagnose the location of the tropopause. Multiple
! algorithms are provided, some of which may not be able to identify a
! tropopause in all situations. To handle these cases, an analytic
! definition and a climatology are provided that can be used to fill in
! when the original algorithm fails. The tropopause temperature and
! pressure are determined and can be output to the history file.
!
! Author: Charles Bardeen
! Created: April, 2009
!
! CCPP-ized: Haipeng Lin, August 2024
module tropopause_find
  !---------------------------------------------------------------
  ! ... variables for the tropopause module
  !---------------------------------------------------------------

  use ccpp_kinds,           only : kind_phys

  implicit none

  private

  ! CCPP-compliant subroutines
  public :: tropopause_find_init
  public :: tropopause_find_run

  ! "Wrapped" routine for use in old CAM for backward compatibility.
  ! Also called by tropopause_find_run driver routine.
  public :: tropopause_findWithBackup

  ! Switches for tropopause method
  public :: TROP_ALG_NONE, TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE
  public :: TROP_ALG_STOBIE, TROP_ALG_HYBSTOB, TROP_ALG_TWMO, TROP_ALG_WMO
  public :: TROP_ALG_CPP
  public :: NOTFOUND

  save

  ! These parameters define an enumeration to be used to define the primary
  ! and backup algorithms to be used with the tropopause_find() method. The
  ! backup algorithm is meant to provide a solution when the primary algorithm
  ! fails. The algorithms that can't fail (i.e., always find a tropopause) are:
  ! TROP_ALG_ANALYTIC, TROP_ALG_CLIMATE, and TROP_ALG_STOBIE.
  integer, parameter    :: TROP_ALG_NONE      = 1    ! Don't evaluate
  integer, parameter    :: TROP_ALG_ANALYTIC  = 2    ! Analytic Expression
  integer, parameter    :: TROP_ALG_CLIMATE   = 3    ! Climatology
  integer, parameter    :: TROP_ALG_STOBIE    = 4    ! Stobie Algorithm
  integer, parameter    :: TROP_ALG_TWMO      = 5    ! WMO Definition, Reichler et al. [2003]
  integer, parameter    :: TROP_ALG_WMO       = 6    ! WMO Definition
  integer, parameter    :: TROP_ALG_HYBSTOB   = 7    ! Hybrid Stobie Algorithm
  integer, parameter    :: TROP_ALG_CPP       = 8    ! Cold Point Parabolic
  integer, parameter    :: TROP_ALG_CHEMTROP  = 9    ! Chemical tropopause

  integer, parameter    :: default_primary    = TROP_ALG_TWMO        ! default primary algorithm
  integer, parameter    :: default_backup     = TROP_ALG_CLIMATE     ! default backup algorithm

  integer, parameter    :: NOTFOUND = -1

  real(kind_phys), parameter :: ALPHA  = 0.03_kind_phys

  ! physical constants
  ! These constants are set in module variables rather than as parameters
  ! to support the aquaplanet mode in which the constants have values determined
  ! by the experiment protocol
  real(kind_phys) :: cnst_kap     ! = cappa
  real(kind_phys) :: cnst_faktor  ! = -gravit/rair
  real(kind_phys) :: cnst_rga     ! = 1/gravit
  real(kind_phys) :: cnst_ka1     ! = cnst_kap - 1._kind_phys
  real(kind_phys) :: cnst_rad2deg ! = 180/pi
  ! Device copies for the kernels and acc routines below; set in _init.
  !$acc declare create(cnst_kap, cnst_faktor, cnst_rga, cnst_ka1, cnst_rad2deg)

!================================================================================================
contains
!================================================================================================

!> \section arg_table_tropopause_find_init Argument Table
!! \htmlinclude tropopause_find_init.html
  subroutine tropopause_find_init(cappa, rair, gravit, pi, errmsg, errflg)

    real(kind_phys), intent(in)         :: cappa   ! R/Cp
    real(kind_phys), intent(in)         :: rair    ! Dry air gas constant (J K-1 kg-1)
    real(kind_phys), intent(in)         :: gravit  ! Gravitational acceleration (m s-2)
    real(kind_phys), intent(in)         :: pi      ! Pi

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    ! define physical constants
    cnst_kap     = cappa
    cnst_faktor  = -gravit/rair
    cnst_rga     = 1._kind_phys/gravit              ! Reciprocal of gravit (s2 m-1)
    cnst_ka1     = cnst_kap - 1._kind_phys
    cnst_rad2deg = 180._kind_phys/pi               ! radians to degrees conversion factor
    !$acc update device(cnst_kap, cnst_faktor, cnst_rga, cnst_ka1, cnst_rad2deg)

  end subroutine tropopause_find_init

  ! "Driver" routine for tropopause_find. Identifies the tropopause using several methods
  ! and populates them into the model state.
  ! Most methods use climatological tropopause as a backup, and as such is guaranteed to
  ! find a tropopause in all columns. Others are explicitly single-method and used by
  ! other parameterizations as-is with NOTFOUND values being intentional.
!> \section arg_table_tropopause_find_run Argument Table
!! \htmlinclude tropopause_find_run.html
  subroutine tropopause_find_run(ncol, pver, fillvalue, lat, pint, pmid, t, zi, zm, phis, &
                                 calday, tropp_p_loc, tropp_days, &
                                 tropLev, tropP, tropT, tropZ, & ! Default primary+backup (twmo+climate)
                                 tropLev_twmo, tropP_twmo, tropT_twmo, tropZ_twmo, & ! Primary only (twmo)
                                 tropLev_clim, tropP_clim, tropT_clim, tropZ_clim, & ! Climate-only
                                 tropLev_hybstob, tropP_hybstob, tropT_hybstob, tropZ_hybstob, & ! Hybridstobie + climate backup
                                 tropLev_cpp, tropP_cpp, tropT_cpp, tropZ_cpp, & ! Cold point only
                                 tropLev_chem, tropP_chem, tropT_chem, tropZ_chem, & ! Chemical tropopause only
                                 hstobie_trop, hstobie_linoz, hstobie_tropop, & ! Hybridstobie only for chemistry diagnostics
                                 scheme_name, errmsg, errflg)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: fillvalue     ! Fill value for diagnostic outputs
    real(kind_phys), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    real(kind_phys), intent(in)         :: calday        ! Day of year including fraction of day

    ! Climatological tropopause pressures (Pa), (ncol,ntimes=12).
    real(kind_phys), intent(in)         :: tropp_p_loc(:,:)
    real(kind_phys), intent(in)         :: tropp_days(:) ! Day-of-year for climo data, 12

    integer,         intent(out)     :: tropLev(:)         ! tropopause level index
    real(kind_phys), intent(out)     :: tropP(:)           ! tropopause pressure (Pa)
    real(kind_phys), intent(out)     :: tropT(:)           ! tropopause temperature (K)
    real(kind_phys), intent(out)     :: tropZ(:)           ! tropopause height (m)

    integer,         intent(out)     :: tropLev_twmo(:)    ! lapse-rate tropopause level index
    real(kind_phys), intent(out)     :: tropP_twmo(:)      ! lapse-rate tropopause pressure (Pa)
    real(kind_phys), intent(out)     :: tropT_twmo(:)      ! lapse-rate tropopause temperature (K)
    real(kind_phys), intent(out)     :: tropZ_twmo(:)      ! lapse-rate tropopause height (m)

    integer,         intent(out)     :: tropLev_clim(:)    ! climatology-backed tropopause level index
    real(kind_phys), intent(out)     :: tropP_clim(:)      ! climatology-backed tropopause pressure (Pa)
    real(kind_phys), intent(out)     :: tropT_clim(:)      ! climatology-backed tropopause temperature (K)
    real(kind_phys), intent(out)     :: tropZ_clim(:)      ! climatology-backed tropopause height (m)

    integer,         intent(out)     :: tropLev_hybstob(:) ! hybridstobie climatology-backed tropopause level index
    real(kind_phys), intent(out)     :: tropP_hybstob(:)   ! hybridstobie climatology-backed tropopause pressure (Pa)
    real(kind_phys), intent(out)     :: tropT_hybstob(:)   ! hybridstobie climatology-backed tropopause temperature (K)
    real(kind_phys), intent(out)     :: tropZ_hybstob(:)   ! hybridstobie climatology-backed tropopause height (m)

    integer,         intent(out)     :: tropLev_cpp(:)     ! cold point tropopause level index
    real(kind_phys), intent(out)     :: tropP_cpp(:)       ! cold point tropopause pressure (Pa)
    real(kind_phys), intent(out)     :: tropT_cpp(:)       ! cold point tropopause temperature (K)
    real(kind_phys), intent(out)     :: tropZ_cpp(:)       ! cold point tropopause height (m)

    integer,         intent(out)     :: tropLev_chem(:)    ! chemical tropopause level index
    real(kind_phys), intent(out)     :: tropP_chem(:)      ! chemical tropopause pressure (Pa)
    real(kind_phys), intent(out)     :: tropT_chem(:)      ! chemical tropopause temperature (K)
    real(kind_phys), intent(out)     :: tropZ_chem(:)      ! chemical tropopause height (m)

    ! Optional output arguments for hybridstobie with chemistry
    real(kind_phys), intent(out)   :: hstobie_trop(:,:)   ! Lowest level with strat. chem
    real(kind_phys), intent(out)   :: hstobie_linoz(:,:)  ! Lowest possible Linoz level
    real(kind_phys), intent(out)   :: hstobie_tropop(:,:) ! Troposphere boundary calculated in chem.

    character(len=64),  intent(out) :: scheme_name
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    scheme_name = 'tropopause_find'
    errmsg = ' '
    errflg = 0

    ! The column work below runs on the device (OpenACC); the routines it
    ! calls find these arrays present. tropp_days and the scalars stay on
    ! the host (the climatology's time weights are computed there).
    !$acc data present(lat, pint, pmid, t, zi, zm, phis, tropp_p_loc, &
    !$acc              tropLev, tropP, tropT, tropZ, &
    !$acc              tropLev_twmo, tropP_twmo, tropT_twmo, tropZ_twmo, &
    !$acc              tropLev_clim, tropP_clim, tropT_clim, tropZ_clim, &
    !$acc              tropLev_hybstob, tropP_hybstob, tropT_hybstob, tropZ_hybstob, &
    !$acc              tropLev_cpp, tropP_cpp, tropT_cpp, tropZ_cpp, &
    !$acc              tropLev_chem, tropP_chem, tropT_chem, tropZ_chem, &
    !$acc              hstobie_trop, hstobie_linoz, hstobie_tropop)

    ! Obtain the primary output, which is TWMO + climate
    call tropopause_findWithBackup( &
         ncol           = ncol, &
         pver           = pver, &
         fillvalue      = fillvalue, &
         lat            = lat, &
         pint           = pint, &
         pmid           = pmid, &
         t              = t, &
         zi             = zi, &
         zm             = zm, &
         phis           = phis, &
         calday         = calday, &
         tropp_p_loc    = tropp_p_loc, &
         tropp_days     = tropp_days, &
         tropLev        = tropLev, &
         tropP          = tropP, &
         tropT          = tropT, &
         tropZ          = tropZ, &
         primary        = default_primary, &
         backup         = default_backup, &
         errmsg         = errmsg, &
         errflg         = errflg &
    )

    ! Any other intended outputs
    ! Primary (TWMO) only
    call tropopause_findWithBackup( &
         ncol           = ncol, &
         pver           = pver, &
         fillvalue      = fillvalue, &
         lat            = lat, &
         pint           = pint, &
         pmid           = pmid, &
         t              = t, &
         zi             = zi, &
         zm             = zm, &
         phis           = phis, &
         calday         = calday, &
         tropp_p_loc    = tropp_p_loc, &
         tropp_days     = tropp_days, &
         tropLev        = tropLev_twmo, &
         tropP          = tropP_twmo, &
         tropT          = tropT_twmo, &
         tropZ          = tropZ_twmo, &
         primary        = default_primary, &
         backup         = TROP_ALG_NONE, &
         errmsg         = errmsg, &
         errflg         = errflg &
    )

    ! Climatology only
    call tropopause_findWithBackup( &
         ncol           = ncol, &
         pver           = pver, &
         fillvalue      = fillvalue, &
         lat            = lat, &
         pint           = pint, &
         pmid           = pmid, &
         t              = t, &
         zi             = zi, &
         zm             = zm, &
         phis           = phis, &
         calday         = calday, &
         tropp_p_loc    = tropp_p_loc, &
         tropp_days     = tropp_days, &
         tropLev        = tropLev_clim, &
         tropP          = tropP_clim, &
         tropT          = tropT_clim, &
         tropZ          = tropZ_clim, &
         primary        = TROP_ALG_CLIMATE, &
         backup         = TROP_ALG_NONE, &
         errmsg         = errmsg, &
         errflg         = errflg &
    )

    ! Cold point (CPP) only
    call tropopause_findWithBackup( &
         ncol           = ncol, &
         pver           = pver, &
         fillvalue      = fillvalue, &
         lat            = lat, &
         pint           = pint, &
         pmid           = pmid, &
         t              = t, &
         zi             = zi, &
         zm             = zm, &
         phis           = phis, &
         calday         = calday, &
         tropp_p_loc    = tropp_p_loc, &
         tropp_days     = tropp_days, &
         tropLev        = tropLev_cpp, &
         tropP          = tropP_cpp, &
         tropT          = tropT_cpp, &
         tropZ          = tropZ_cpp, &
         primary        = TROP_ALG_CPP, &
         backup         = TROP_ALG_NONE, &
         errmsg         = errmsg, &
         errflg         = errflg &
    )

    ! Hybridstobie with climatology-backed
    call tropopause_findWithBackup( &
         ncol           = ncol, &
         pver           = pver, &
         fillvalue      = fillvalue, &
         lat            = lat, &
         pint           = pint, &
         pmid           = pmid, &
         t              = t, &
         zi             = zi, &
         zm             = zm, &
         phis           = phis, &
         calday         = calday, &
         tropp_p_loc    = tropp_p_loc, &
         tropp_days     = tropp_days, &
         tropLev        = tropLev_hybstob, &
         tropP          = tropP_hybstob, &
         tropT          = tropT_hybstob, &
         tropZ          = tropZ_hybstob, &
         primary        = TROP_ALG_HYBSTOB, &
         backup         = TROP_ALG_CLIMATE, &
         hstobie_trop   = hstobie_trop, &    ! Only used if TROP_ALG_HYBSTOB
         hstobie_linoz  = hstobie_linoz, &   ! Only used if TROP_ALG_HYBSTOB
         hstobie_tropop = hstobie_tropop, &  ! Only used if TROP_ALG_HYBSTOB
         errmsg         = errmsg, &
         errflg         = errflg &
    )

    ! Chemical tropopause (used for chemistry)
    call tropopause_findWithBackup( &
         ncol           = ncol, &
         pver           = pver, &
         fillvalue      = fillvalue, &
         lat            = lat, &
         pint           = pint, &
         pmid           = pmid, &
         t              = t, &
         zi             = zi, &
         zm             = zm, &
         phis           = phis, &
         calday         = calday, &
         tropp_p_loc    = tropp_p_loc, &
         tropp_days     = tropp_days, &
         tropLev        = tropLev_chem, &
         tropP          = tropP_chem, &
         tropT          = tropT_chem, &
         tropZ          = tropZ_chem, &
         primary        = TROP_ALG_CHEMTROP, &
         backup         = TROP_ALG_CLIMATE, &
         errmsg         = errmsg, &
         errflg         = errflg &
    )

    !$acc end data

  end subroutine tropopause_find_run

  ! Searches all the columns and attempts to identify the tropopause.
  ! Two routines can be specifed, a primary routine which is tried first and a
  ! backup routine which will be tried only if the first routine fails. If the
  ! tropopause can not be identified by either routine, then a NOTFOUND is returned
  ! for the tropopause level, temperature and pressure.
  subroutine tropopause_findWithBackup(ncol, pver, fillvalue, lat, pint, pmid, t, zi, zm, phis, &
                                       calday, tropp_p_loc, tropp_days, &
                                       tropLev, tropP, tropT, tropZ, &
                                       hstobie_trop, hstobie_linoz, hstobie_tropop, &
                                       primary, backup, &
                                       errmsg, errflg)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: fillvalue     ! Fill value for diagnostic outputs
    real(kind_phys), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    real(kind_phys), intent(in)         :: calday        ! Day of year including fraction of day

    ! Climatological tropopause pressures (Pa), (ncol,ntimes=12).
    real(kind_phys), intent(in)         :: tropp_p_loc(:,:)
    real(kind_phys), intent(in)         :: tropp_days(:) ! Day-of-year for climo data, 12

    integer,                   intent(out)     :: tropLev(:)          ! tropopause level index
    real(kind_phys), optional, intent(out)     :: tropP(:)            ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(out)     :: tropT(:)            ! tropopause temperature (K)
    real(kind_phys), optional, intent(out)     :: tropZ(:)            ! tropopause height (m)

    ! Optional output arguments for hybridstobie with chemistry
    real(kind_phys), optional, intent(inout)   :: hstobie_trop(:,:)   ! Lowest level with strat. chem
    real(kind_phys), optional, intent(inout)   :: hstobie_linoz(:,:)  ! Lowest possible Linoz level
    real(kind_phys), optional, intent(inout)   :: hstobie_tropop(:,:) ! Troposphere boundary calculated in chem.

    ! primary and backup are no longer optional arguments for CCPP-compliance.
    ! specify defaults when calling (TWMO, CLIMO)
    integer, intent(in)       :: primary                   ! primary detection algorithm
    integer, intent(in)       :: backup                    ! backup detection algorithm

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local variables
    integer :: i
    logical :: do_p, do_t, do_z
    logical :: do_hstobie       ! pass the hybrid-Stobie chemistry outputs (all three present)
    integer :: last, next       ! climatology months bracketing calday
    real(kind_phys) :: dels     ! climatology time weight

    errmsg = ' '
    errflg = 0

    ! Fused path: when both algorithms have a column form (the ones
    ! tropopause_find_run uses), one kernel does, per column, the fill, the
    ! primary algorithm and the backup algorithm if this column was not
    ! found. Every algorithm only acts on columns still NOTFOUND (and the
    ! hybrid-Stobie outputs are written per column), so this is what the
    ! separate kernels and the any-NOTFOUND test before the backup did,
    ! column by column. CHEMTROP as a backup keeps the separate path: its
    ! latitude check acts on columns already found.
    if (column_algorithm(primary) .and. column_algorithm(backup) .and. &
        backup /= TROP_ALG_CHEMTROP) then
      do_p = present(tropP)
      do_t = present(tropT)
      do_z = present(tropZ)
      do_hstobie = present(hstobie_trop) .and. present(hstobie_linoz) .and. &
                   present(hstobie_tropop)
      call climate_weights(calday, tropp_days, last, next, dels)

      !$acc parallel loop gang vector &
      !$acc   copyin(lat, pint, pmid, t, zi, zm, phis, tropp_p_loc) &
      !$acc   copy(tropLev, tropP, tropT, tropZ, hstobie_trop, hstobie_linoz, hstobie_tropop)
      do i = 1, size(tropLev)
        ! Initialize the results to a missing value, so that the algorithms
        ! will attempt to find the tropopause for all of them.
        tropLev(i) = NOTFOUND
        if (do_p) tropP(i) = fillvalue
        if (do_t) tropT(i) = fillvalue
        if (do_z) tropZ(i) = fillvalue
        if (i <= ncol) then
          ! Try to find the tropopause using the primary algorithm.
          if (do_hstobie) then
            call tropopause_column_using(primary, i, pver, lat, pint, pmid, t, zi, zm, phis, &
                 tropp_p_loc, last, next, dels, tropLev, tropP, tropT, tropZ, &
                 hstobie_trop, hstobie_linoz, hstobie_tropop)
          else
            call tropopause_column_using(primary, i, pver, lat, pint, pmid, t, zi, zm, phis, &
                 tropp_p_loc, last, next, dels, tropLev, tropP, tropT, tropZ)
          end if
          ! Then the backup algorithm, for this column if it was not found.
          if (backup /= TROP_ALG_NONE .and. tropLev(i) == NOTFOUND) then
            call tropopause_column_using(backup, i, pver, lat, pint, pmid, t, zi, zm, phis, &
                 tropp_p_loc, last, next, dels, tropLev, tropP, tropT, tropZ)
          end if
        end if
      end do
      return
    end if

    ! Separate path (an algorithm without a column form: ANALYTIC, STOBIE,
    ! WMO, which are host code).
    ! Initialize the results to a missing value, so that the algorithms will
    ! attempt to find the tropopause for all of them.
    ! (OpenACC: whether an optional output is present is decided here on the
    ! host; a data clause on an absent optional argument is ignored.)
    do_p = present(tropP)
    do_t = present(tropT)
    do_z = present(tropZ)
    !$acc parallel loop gang vector copy(tropLev)
    do i = 1, size(tropLev)
      tropLev(i) = NOTFOUND
    end do
    if (do_p) then
      !$acc parallel loop gang vector copy(tropP)
      do i = 1, size(tropP)
        tropP(i) = fillvalue
      end do
    end if
    if (do_t) then
      !$acc parallel loop gang vector copy(tropT)
      do i = 1, size(tropT)
        tropT(i) = fillvalue
      end do
    end if
    if (do_z) then
      !$acc parallel loop gang vector copy(tropZ)
      do i = 1, size(tropZ)
        tropZ(i) = fillvalue
      end do
    end if

    ! Try to find the tropopause using the primary algorithm.
    if (primary /= TROP_ALG_NONE) then
      call tropopause_findUsing(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                calday, tropp_p_loc, tropp_days, &
                                primary, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, &
                                hstobie_trop=hstobie_trop, hstobie_linoz=hstobie_linoz, hstobie_tropop=hstobie_tropop, & ! only for HYBSTOB
                                errmsg=errmsg, errflg=errflg)
    end if

    ! (any_notfound: the any(tropLev(:) == NOTFOUND) test, evaluated where
    ! tropLev lives -- on the device in GPU builds.)
    if (backup /= TROP_ALG_NONE) then
      if (any_notfound(tropLev)) then
        call tropopause_findUsing(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                  calday, tropp_p_loc, tropp_days, &
                                  backup, tropLev, tropP=tropP, tropT=tropT, tropZ=tropZ, &
                                  errmsg=errmsg, errflg=errflg)
      end if
    end if

  end subroutine tropopause_findWithBackup

  ! True if any column has no tropopause yet: any(tropLev(:) == NOTFOUND),
  ! as an OpenACC reduction so that in GPU builds only the result, not
  ! tropLev, comes back to the host.
  logical function any_notfound(tropLev)

    integer, intent(in) :: tropLev(:)          ! tropopause level index

    integer :: i
    logical :: found_missing

    found_missing = .false.
    !$acc parallel loop gang vector reduction(.or.:found_missing) copyin(tropLev)
    do i = 1, size(tropLev)
      found_missing = found_missing .or. (tropLev(i) == NOTFOUND)
    end do
    any_notfound = found_missing

  end function any_notfound

  ! True if the algorithm has a column form (tropopause_column_using).
  pure logical function column_algorithm(algorithm)

    integer, intent(in) :: algorithm         ! detection algorithm

    column_algorithm = algorithm == TROP_ALG_NONE .or. algorithm == TROP_ALG_CLIMATE .or. &
                       algorithm == TROP_ALG_TWMO .or. algorithm == TROP_ALG_CPP .or.     &
                       algorithm == TROP_ALG_HYBSTOB .or. algorithm == TROP_ALG_CHEMTROP

  end function column_algorithm

  ! Time interpolation of the tropopause climatology: the months bracketing
  ! calday and the weight of the later one.
  subroutine climate_weights(calday, tropp_days, last, next, dels)

    real(kind_phys), intent(in)  :: calday        ! Day of year including fraction of day
    real(kind_phys), intent(in)  :: tropp_days(:) ! Day-of-year for climo data, 12
    integer,         intent(out) :: last, next    ! bracketing months
    real(kind_phys), intent(out) :: dels          ! weight of month next

    integer :: m

    !--------------------------------------------------------
    ! ... setup the time interpolation
    !--------------------------------------------------------
    if( calday < tropp_days(1) ) then
      next = 1
      last = 12
      dels = (365._kind_phys + calday - tropp_days(12)) / (365._kind_phys + tropp_days(1) - tropp_days(12))
    else if( calday >= tropp_days(12) ) then
      next = 1
      last = 12
      dels = (calday - tropp_days(12)) / (365._kind_phys + tropp_days(1) - tropp_days(12))
    else
      do m = 11,1,-1
         if( calday >= tropp_days(m) ) then
            exit
         end if
      end do
      last = m
      next = m + 1
      dels = (calday - tropp_days(m)) / (tropp_days(m+1) - tropp_days(m))
    end if

    dels = max( min( 1._kind_phys,dels ),0._kind_phys )

  end subroutine climate_weights

  ! Column form of tropopause_findUsing for the algorithms with one (see
  ! column_algorithm): applies the algorithm to column icol, which it
  ! changes only if tropLev(icol) is NOTFOUND (hybrid Stobie also writes its
  ! chemistry outputs of the column, if present).
  subroutine tropopause_column_using(algorithm, icol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                     tropp_p_loc, last, next, dels, tropLev, tropP, tropT, tropZ, &
                                     hstobie_trop, hstobie_linoz, hstobie_tropop)
    !$acc routine seq

    integer,         intent(in)         :: algorithm     ! detection algorithm
    integer,         intent(in)         :: icol          ! column index
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    real(kind_phys), intent(in)         :: tropp_p_loc(:,:) ! Climatological tropopause pressures (Pa)
    integer,         intent(in)         :: last, next    ! climatology months (climate_weights)
    real(kind_phys), intent(in)         :: dels          ! climatology time weight
    integer,                   intent(inout) :: tropLev(:)  ! tropopause level index
    real(kind_phys), optional, intent(inout) :: tropP(:)    ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout) :: tropT(:)    ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout) :: tropZ(:)    ! tropopause height (m)
    real(kind_phys), optional, intent(inout) :: hstobie_trop(:,:)   ! Lowest level with strat. chem
    real(kind_phys), optional, intent(inout) :: hstobie_linoz(:,:)  ! Lowest possible Linoz level
    real(kind_phys), optional, intent(inout) :: hstobie_tropop(:,:) ! Troposphere boundary calculated in chem.

    real(kind_phys) :: dlats                ! latitude (degrees), this column

    select case (algorithm)
    case (TROP_ALG_CLIMATE)
      call tropopause_climate_column(icol, pver, pint, pmid, t, zi, zm, phis, &
                                     tropp_p_loc, last, next, dels, tropLev, tropP, tropT, tropZ)
    case (TROP_ALG_HYBSTOB)
      call tropopause_hybridstobie_column(icol, pver, pmid, t, zm, tropLev, tropP, tropT, tropZ, &
                                          hstobie_trop, hstobie_linoz, hstobie_tropop)
    case (TROP_ALG_TWMO)
      call tropopause_twmo_column(icol, pver, pint, pmid, t, zi, zm, phis, &
                                  tropLev, tropP, tropT, tropZ)
    case (TROP_ALG_CPP)
      call tropopause_cpp_column(icol, pver, pmid, t, zm, tropLev, tropP, tropT, tropZ)
    case (TROP_ALG_CHEMTROP)
      ! tropopause_findChemTrop for one column: the lapse rate tropopause,
      ! then at high latitudes (poleward of 50) the level is reset if found
      ! at P <= 125 hPa, then the climatology if not found.
      call tropopause_twmo_column(icol, pver, pint, pmid, t, zi, zm, phis, &
                                  tropLev, tropP, tropT, tropZ)
      dlats = lat(icol) * cnst_rad2deg ! convert to degrees
      if (abs(dlats) > 50._kind_phys) then
        if (tropLev(icol) .ne. NOTFOUND) then
          if (pmid(icol, tropLev(icol)) <= 12500._kind_phys) then
            tropLev(icol) = NOTFOUND
          end if
        end if
      end if
      call tropopause_climate_column(icol, pver, pint, pmid, t, zi, zm, phis, &
                                     tropp_p_loc, last, next, dels, tropLev, tropP, tropT, tropZ)
    end select

  end subroutine tropopause_column_using

  ! Call the appropriate tropopause detection routine based upon the algorithm
  ! specifed.
  !
  ! NOTE: It is assumed that the output fields have been initialized by the
  ! caller, and only output values set to fillvalue will be detected.
  subroutine tropopause_findUsing(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                  calday, tropp_p_loc, tropp_days, &
                                  algorithm, tropLev, tropP, tropT, tropZ, &
                                  hstobie_trop, hstobie_linoz, hstobie_tropop, &
                                  errmsg, errflg)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    real(kind_phys), intent(in)         :: calday        ! Day of year including fraction of day

    ! Climatological tropopause pressures (Pa), (ncol,ntimes=12).
    real(kind_phys), intent(in)         :: tropp_p_loc(:,:)
    real(kind_phys), intent(in)         :: tropp_days(:) ! Day-of-year for climo data, 12

    integer,                   intent(in)      :: algorithm             ! detection algorithm
    integer,                   intent(inout)   :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)   :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)   :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)   :: tropZ(:)              ! tropopause height (m)

    ! Optional output arguments for hybridstobie with chemistry
    real(kind_phys), optional, intent(inout)   :: hstobie_trop(:,:)   ! Lowest level with strat. chem
    real(kind_phys), optional, intent(inout)   :: hstobie_linoz(:,:)  ! Lowest possible Linoz level
    real(kind_phys), optional, intent(inout)   :: hstobie_tropop(:,:) ! Troposphere boundary calculated in chem.

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ' '
    errflg = 0

    ! Dispatch the request to the appropriate routine.
    ! OpenACC: the algorithms tropopause_find_run uses (CLIMATE, HYBSTOB,
    ! TWMO, CPP, CHEMTROP) run on the device; ANALYTIC, STOBIE and WMO are
    ! host code, for callers of tropopause_findWithBackup with host data.
    select case(algorithm)
      case(TROP_ALG_ANALYTIC)
        call tropopause_analytic(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_CLIMATE)
        call tropopause_climate(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                calday, tropp_p_loc, tropp_days, &
                                tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_STOBIE)
        call tropopause_stobie(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_HYBSTOB)
        if(present(hstobie_trop) .and. present(hstobie_linoz) .and. present(hstobie_tropop)) then
          call tropopause_hybridstobie(ncol, pver, pmid, t, zm, &
                                       tropLev, tropP, tropT, tropZ, &
                                       hstobie_trop, hstobie_linoz, hstobie_tropop)
        else
          call tropopause_hybridstobie(ncol, pver, pmid, t, zm, &
                                       tropLev, tropP, tropT, tropZ)
        endif

      case(TROP_ALG_TWMO)
        call tropopause_twmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_WMO)
        call tropopause_wmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_CPP)
        call tropopause_cpp(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

      case(TROP_ALG_CHEMTROP)
        ! hplin: needs climatological arguments as calling tropopause_climate from within findChemTrop
        call tropopause_findChemTrop(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                     calday, tropp_p_loc, tropp_days, &
                                     tropLev, tropP, tropT, tropZ, &
                                     errmsg, errflg)

      case default
        errflg = 1
        write(errmsg,*) 'tropopause: Invalid detection algorithm (',  algorithm, ') specified.'
    end select

  end subroutine tropopause_findUsing

  ! This analytic expression closely matches the mean tropopause determined
  ! by the NCEP reanalysis and has been used by the radiation code.
  subroutine tropopause_analytic(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                 tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer,                   intent(inout)   :: tropLev(:)             ! tropopause level index
    real(kind_phys), optional, intent(inout)   :: tropP(:)               ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)   :: tropT(:)               ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)   :: tropZ(:)               ! tropopause height (m)

    ! Local Variables
    integer       :: i
    integer       :: k
    real(kind_phys)      :: tP                       ! tropopause pressure (Pa)

    ! Iterate over all of the columns.
    do i = 1, ncol

      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then

        ! Calculate the pressure of the tropopause.
        tP = (25000.0_kind_phys - 15000.0_kind_phys * (cos(lat(i)))**2)

        ! Find the level that contains the tropopause.
        do k = pver, 2, -1
          if (tP >= pint(i, k)) then
            tropLev(i) = k
            exit
          end if
        end do

        ! Return the optional outputs
        if (present(tropP)) tropP(i) = tP

        if (present(tropT)) then
          tropT(i) = tropopause_interpolateT(pver, pmid, t, i, tropLev(i), tP)
        end if

        if (present(tropZ)) then
          tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
        end if
      end if
    end do

  end subroutine tropopause_analytic

  ! Determine the tropopause pressure from a climatology,
  ! interpolated to the current day of year and latitude.
  subroutine tropopause_climate(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                calday, tropp_p_loc, tropp_days, tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    real(kind_phys), intent(in)         :: calday        ! Day of year including fraction of day

    ! Climatological tropopause pressures (Pa), (ncol,ntimes=12).
    real(kind_phys), intent(in)         :: tropp_p_loc(:,:)
    real(kind_phys), intent(in)         :: tropp_days(:) ! Day-of-year for climo data, 12

    integer,                   intent(inout)  :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    ! Local Variables
    integer       :: i
    real(kind_phys)      :: dels
    integer       :: last
    integer       :: next

    ! If any columns remain to be indentified, then get the current
    ! day from the calendar.

    if (any_notfound(tropLev)) then

      call climate_weights(calday, tropp_days, last, next, dels)

      ! Iterate over all of the columns.
      !$acc parallel loop gang vector &
      !$acc   copyin(tropp_p_loc, pint, pmid, t, zi, zm, phis) &
      !$acc   copy(tropLev, tropP, tropT, tropZ)
      do i = 1, ncol
        call tropopause_climate_column(i, pver, pint, pmid, t, zi, zm, phis, &
                                       tropp_p_loc, last, next, dels, tropLev, tropP, tropT, tropZ)
      end do
    end if

  end subroutine tropopause_climate

  ! tropopause_climate for column icol (months and weight from
  ! climate_weights).
  subroutine tropopause_climate_column(icol, pver, pint, pmid, t, zi, zm, phis, &
                                       tropp_p_loc, last, next, dels, tropLev, tropP, tropT, tropZ)
    !$acc routine seq

    integer,         intent(in)         :: icol          ! column index
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa), pverp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    real(kind_phys), intent(in)         :: tropp_p_loc(:,:) ! Climatological tropopause pressures (Pa), (ncol,12)
    integer,         intent(in)         :: last, next    ! bracketing months
    real(kind_phys), intent(in)         :: dels          ! weight of month next
    integer,                   intent(inout)  :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    integer         :: k
    real(kind_phys) :: tP                       ! tropopause pressure (Pa)

    ! Skip column in which the tropopause has already been found.
    if (tropLev(icol) == NOTFOUND) then

    !--------------------------------------------------------
    ! ... get tropopause level from climatology
    !--------------------------------------------------------
      ! Interpolate the tropopause pressure.
      tP = tropp_p_loc(icol,last) &
        + dels * (tropp_p_loc(icol,next) - tropp_p_loc(icol,last))

      ! Find the associated level.
      do k = pver, 2, -1
        if (tP >= pint(icol, k)) then
          tropLev(icol) = k
          exit
        end if
      end do

      ! Return the optional outputs
      if (present(tropP)) tropP(icol) = tP

      if (present(tropT)) then
        tropT(icol) = tropopause_interpolateT(pver, pmid, t, icol, tropLev(icol), tP)
      end if

      if (present(tropZ)) then
        tropZ(icol) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, icol, tropLev(icol), tP)
      end if
    end if

  end subroutine tropopause_climate_column

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine tropopause_hybridstobie(ncol, pver, pmid, t, zm, &
                                     tropLev, tropP, tropT, tropZ, &
                                     hstobie_trop, hstobie_linoz, hstobie_tropop)

    !-----------------------------------------------------------------------
    ! Originally written by Philip Cameron-Smith, LLNL
    !
    !   Stobie-Linoz hybrid: the highest altitude of
    !          a) Stobie algorithm, or
    !          b) minimum Linoz pressure.
    !
    ! NOTE: the ltrop(i) gridbox itself is assumed to be a STRATOSPHERIC gridbox.
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !        ... Local variables
    !-----------------------------------------------------------------------
    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levelserp
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)

    integer,            intent(inout)   :: tropLev(:)             ! tropopause level index
    real(kind_phys), optional, intent(inout)   :: tropP(:)               ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)   :: tropT(:)               ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)   :: tropZ(:)               ! tropopause height (m)

    ! Optional output arguments for hybridstobie with chemistry
    real(kind_phys), optional, intent(inout)   :: hstobie_trop(:,:)   ! Lowest level with strat. chem
    real(kind_phys), optional, intent(inout)   :: hstobie_linoz(:,:)  ! Lowest possible Linoz level
    real(kind_phys), optional, intent(inout)   :: hstobie_tropop(:,:) ! Troposphere boundary calculated in chem.

    integer      :: i

    ! One thread per column (tropopause_hybridstobie_column), which writes
    ! its own rows of the (ncol, pver) 0/1 outputs.
    !$acc parallel loop gang vector &
    !$acc   copyin(pmid, t, zm) &
    !$acc   copy(tropLev, tropP, tropT, tropZ, hstobie_trop, hstobie_linoz, hstobie_tropop)
    do i=1,ncol
      call tropopause_hybridstobie_column(i, pver, pmid, t, zm, tropLev, tropP, tropT, tropZ, &
                                          hstobie_trop, hstobie_linoz, hstobie_tropop)
    end do

  end subroutine tropopause_hybridstobie

  ! tropopause_hybridstobie for column icol.
  subroutine tropopause_hybridstobie_column(icol, pver, pmid, t, zm, &
                                            tropLev, tropP, tropT, tropZ, &
                                            hstobie_trop, hstobie_linoz, hstobie_tropop)
    !$acc routine seq

    integer,         intent(in)         :: icol          ! column index
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)

    integer,            intent(inout)   :: tropLev(:)             ! tropopause level index
    real(kind_phys), optional, intent(inout)   :: tropP(:)               ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)   :: tropT(:)               ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)   :: tropZ(:)               ! tropopause height (m)

    ! Optional output arguments for hybridstobie with chemistry
    real(kind_phys), optional, intent(inout)   :: hstobie_trop(:,:)   ! Lowest level with strat. chem
    real(kind_phys), optional, intent(inout)   :: hstobie_linoz(:,:)  ! Lowest possible Linoz level
    real(kind_phys), optional, intent(inout)   :: hstobie_tropop(:,:) ! Troposphere boundary calculated in chem.

    real(kind_phys),parameter  ::  min_Stobie_Pressure= 40.E2_kind_phys !For case 2 & 4.  [Pa]
    real(kind_phys),parameter  ::  max_Linoz_Pressure =208.E2_kind_phys !For case     4.  [Pa]

    integer      :: k
    real(kind_phys)     :: stobie_min, shybrid_temp      !temporary variable for case 2 & 3.
    integer      :: ltrop_linoz                  !Lowest possible Linoz vertical level (this column)
    integer      :: ltrop_trop                   !Tropopause level for hybrid case (this column)
    logical      :: ltrop_linoz_set               !Flag that lowest linoz level already found.

    ltrop_linoz = 1  ! Initialize to default value.
    ltrop_trop = 1   ! Initialize to default value.

    ! Skip column in which the tropopause has already been found.
    not_found: if (tropLev(icol) == NOTFOUND) then

       stobie_min = 1.e10_kind_phys    ! An impossibly large number
       ltrop_linoz_set = .FALSE.
       LOOP_LEV: do k=pver,1,-1
          IF (pmid(icol,k) < min_stobie_pressure) cycle
          shybrid_temp = ALPHA * t(icol,k) - Log10(pmid(icol,k))
          !PJC_NOTE: the units of pmid won't matter, because it is just an additive offset.
          IF (shybrid_temp<stobie_min) then
             ltrop_trop=k
             stobie_min = shybrid_temp
          ENDIF
          IF (pmid(icol,k) < max_Linoz_pressure .AND. .NOT. ltrop_linoz_set) THEN
             ltrop_linoz = k
             ltrop_linoz_set = .TRUE.
          ENDIF
       enddo LOOP_LEV

       tropLev(icol) = MIN(ltrop_trop,ltrop_linoz)

       if (present(tropP)) then
          tropP(icol) = pmid(icol,tropLev(icol))
       endif
       if (present(tropT)) then
          tropT(icol) = t(icol,tropLev(icol))
       endif
       if (present(tropZ)) then
          tropZ(icol) = zm(icol,tropLev(icol))
       endif

    endif not_found

    ! For output purposes only.
    if (present(hstobie_trop)) then
       do k = 1, pver
          hstobie_trop(icol,k) = 0._kind_phys
       end do
       if (tropLev(icol)>0) hstobie_trop(icol,tropLev(icol)) = 1._kind_phys
    endif
    if (present(hstobie_linoz)) then
       do k = 1, pver
          hstobie_linoz(icol,k) = 0._kind_phys
       end do
       if (tropLev(icol)>0) hstobie_linoz(icol,ltrop_linoz) = 1._kind_phys
    endif
    if (present(hstobie_tropop)) then
       do k = 1, pver
          hstobie_tropop(icol,k) = 0._kind_phys
       end do
       if (tropLev(icol)>0) hstobie_tropop(icol,ltrop_trop) = 1._kind_phys
    endif

  end subroutine tropopause_hybridstobie_column

  ! This routine originates with Stobie at NASA Goddard, but does not have a
  ! known reference. It was supplied by Philip Cameron-Smith of LLNL.
  !
  subroutine tropopause_stobie(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                               tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    integer,            intent(inout)   :: tropLev(:)             ! tropopause level index
    real(kind_phys), optional, intent(inout)   :: tropP(:)               ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)   :: tropT(:)               ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)   :: tropZ(:)               ! tropopause height (m)

    ! Local Variables
    integer       :: i
    integer       :: k
    integer       :: tLev                     ! tropopause level
    real(kind_phys)      :: tP                       ! tropopause pressure (Pa)
    real(kind_phys)      :: stobie(pver)             ! stobie weighted temperature
    real(kind_phys)      :: sTrop                    ! stobie value at the tropopause

    ! Iterate over all of the columns.
    do i = 1, ncol

      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then

        ! Caclulate a pressure weighted temperature.
        stobie(:) = ALPHA * t(i,:) - log10(pmid(i, :))

        ! Search from the bottom up, looking for the first minimum.
        tLev  = -1

        do k = pver-1, 1, -1

          if (pmid(i, k) <= 4000._kind_phys) then
            exit
          end if

          if (pmid(i, k) >= 55000._kind_phys) then
            cycle
          end if

          if ((tLev == -1) .or. (stobie(k) < sTrop)) then
            tLev  = k
            tP    = pmid(i, k)
            sTrop = stobie(k)
          end if
        end do

        if (tLev /= -1) then
          tropLev(i) = tLev

          ! Return the optional outputs
          if (present(tropP)) tropP(i) = tP

          if (present(tropT)) then
            tropT(i) = tropopause_interpolateT(pver, pmid, t, i, tropLev(i), tP)
          end if

          if (present(tropZ)) then
            tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
          end if
        end if
      end if
    end do

  end subroutine tropopause_stobie


  ! This routine is an implementation of Reichler et al. [2003] done by
  ! Reichler and downloaded from his web site. Minimal modifications were
  ! made to have the routine work within the CAM framework (i.e. using
  ! CAM constants and types).
  !
  ! NOTE: I am not a big fan of the goto's and multiple returns in this
  ! code, but for the moment I have left them to preserve as much of the
  ! original and presumably well tested code as possible.
  ! UPDATE: The most "obvious" substitutions have been made to replace
  ! goto/return statements with cycle/exit. The structure is still
  ! somewhat tangled.
  ! UPDATE 2: "gamma" renamed to "gam" in order to avoid confusion
  ! with the Fortran 2008 intrinsic. "level" argument removed because
  ! a physics column is not contiguous, so using explicit dimensions
  ! will cause the data to be needlessly copied.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! determination of tropopause height from gridded temperature data
  !
  ! Reichler, T., M. Dameris, and R. Sausen (2003),
  ! Determining the tropopause height from gridded data,
  ! Geophys. Res. Lett., 30, 2042, doi:10.1029/2003GL018240, 20.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine twmo(t, p, icol, plimu, pliml, gam, trp)
    !$acc routine seq

    ! t, p: the (ncol, pver) arrays; this column is icol. (Formerly the row
    ! sections t(i,:), p(i,:) were passed: strided, so copied per call, and
    ! on the device a temporary per thread. Same arithmetic.)
    real(kind_phys), intent(in), dimension(:,:)    :: t, p
    integer,         intent(in)                    :: icol
    real(kind_phys), intent(in)                    :: plimu, pliml, gam
    real(kind_phys), intent(out)                   :: trp

    real(kind_phys), parameter                     :: deltaz = 2000.0_kind_phys

    real(kind_phys)                                :: pmk, pm, a, b, tm, dtdp, dtdz
    real(kind_phys)                                :: ag, bg, ptph
    real(kind_phys)                                :: pm0, pmk0, dtdz0
    real(kind_phys)                                :: p2km, asum, aquer
    real(kind_phys)                                :: pmk2, pm2, a2, b2, tm2, dtdp2, dtdz2
    integer                                 :: level
    integer                                 :: icount, jj
    integer                                 :: j


    trp=-99.0_kind_phys                           ! negative means not valid

    ! initialize start level
    ! dt/dz
    level = size(t, 2)
    pmk= .5_kind_phys * (p(icol,level-1)**cnst_kap+p(icol,level)**cnst_kap)
    pm = pmk**(1/cnst_kap)
    a = (t(icol,level-1)-t(icol,level))/(p(icol,level-1)**cnst_kap-p(icol,level)**cnst_kap)
    b = t(icol,level)-(a*p(icol,level)**cnst_kap)
    tm = a * pmk + b
    dtdp = a * cnst_kap * (pm**cnst_ka1)
    dtdz = cnst_faktor*dtdp*pm/tm

    main_loop: do j=level-1,2,-1
      pm0 = pm
      pmk0 = pmk
      dtdz0  = dtdz

      ! dt/dz
      pmk= .5_kind_phys * (p(icol,j-1)**cnst_kap+p(icol,j)**cnst_kap)
      pm = pmk**(1/cnst_kap)
      a = (t(icol,j-1)-t(icol,j))/(p(icol,j-1)**cnst_kap-p(icol,j)**cnst_kap)
      b = t(icol,j)-(a*p(icol,j)**cnst_kap)
      tm = a * pmk + b
      dtdp = a * cnst_kap * (pm**cnst_ka1)
      dtdz = cnst_faktor*dtdp*pm/tm
      ! dt/dz valid?
      if (dtdz.le.gam) cycle main_loop    ! no, dt/dz < -2 K/km
      if (pm.gt.plimu)   cycle main_loop    ! no, too low

      ! dtdz is valid, calculate tropopause pressure
      if (dtdz0.lt.gam) then
        ag = (dtdz-dtdz0) / (pmk-pmk0)
        bg = dtdz0 - (ag * pmk0)
        ptph = exp(log((gam-bg)/ag)/cnst_kap)
      else
        ptph = pm
      endif

      if (ptph.lt.pliml) cycle main_loop
      if (ptph.gt.plimu) cycle main_loop

      ! 2nd test: dtdz above 2 km must not exceed gam
      p2km = ptph + deltaz*(pm/tm)*cnst_faktor     ! p at ptph + 2km
      asum = 0.0_kind_phys                                ! dtdz above
      icount = 0                                   ! number of levels above

      ! test until apm < p2km
      in_loop: do jj=j,2,-1

        pmk2 = .5_kind_phys * (p(icol,jj-1)**cnst_kap+p(icol,jj)**cnst_kap) ! p mean ^kappa
        pm2 = pmk2**(1/cnst_kap)                           ! p mean
        if(pm2.gt.ptph) cycle in_loop            ! doesn't happen
        if(pm2.lt.p2km) exit in_loop             ! ptropo is valid

        a2 = (t(icol,jj-1)-t(icol,jj))                     ! a
        a2 = a2/(p(icol,jj-1)**cnst_kap-p(icol,jj)**cnst_kap)
        b2 = t(icol,jj)-(a2*p(icol,jj)**cnst_kap)          ! b
        tm2 = a2 * pmk2 + b2                     ! T mean
        dtdp2 = a2 * cnst_kap * (pm2**(cnst_kap-1))  ! dt/dp
        dtdz2 = cnst_faktor*dtdp2*pm2/tm2
        asum = asum+dtdz2
        icount = icount+1
        aquer = asum/float(icount)               ! dt/dz mean

        ! discard ptropo ?
        if (aquer.le.gam) cycle main_loop      ! dt/dz above < gam

      enddo in_loop  ! test next level

      trp = ptph
      exit main_loop
    enddo main_loop

  end subroutine twmo


  ! This routine uses an implementation of Reichler et al. [2003] done by
  ! Reichler and downloaded from his web site. This is similar to the WMO
  ! routines, but is designed for GCMs with a coarse vertical grid.
  subroutine tropopause_twmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                             tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    integer                 :: i

    ! Iterate over all of the columns.
    !$acc parallel loop gang vector &
    !$acc   copyin(pint, pmid, t, zi, zm, phis) &
    !$acc   copy(tropLev, tropP, tropT, tropZ)
    do i = 1, ncol
      call tropopause_twmo_column(i, pver, pint, pmid, t, zi, zm, phis, &
                                  tropLev, tropP, tropT, tropZ)
    end do

  end subroutine tropopause_twmo

  ! tropopause_twmo for column icol.
  subroutine tropopause_twmo_column(icol, pver, pint, pmid, t, zi, zm, phis, &
                                    tropLev, tropP, tropT, tropZ)
    !$acc routine seq

    integer,         intent(in)         :: icol          ! column index
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    ! Local Variables
    real(kind_phys), parameter     :: gam    = -0.002_kind_phys         ! K/m
    real(kind_phys), parameter     :: plimu    = 45000._kind_phys         ! Pa
    real(kind_phys), parameter     :: pliml    = 7500._kind_phys          ! Pa

    integer                 :: k
    real(kind_phys)                :: tP                       ! tropopause pressure (Pa)

    ! Skip column in which the tropopause has already been found.
    if (tropLev(icol) == NOTFOUND) then

      ! Use the routine from Reichler.
      call twmo(t, pmid, icol, plimu, pliml, gam, tP)

      ! if successful, store of the results and find the level and temperature.
      if (tP > 0) then

        ! Find the associated level.
        do k = pver, 2, -1
          if (tP >= pint(icol, k)) then
            tropLev(icol) = k
            exit
          end if
        end do

        ! Return the optional outputs
        if (present(tropP)) tropP(icol) = tP

        if (present(tropT)) then
          tropT(icol) = tropopause_interpolateT(pver, pmid, t, icol, tropLev(icol), tP)
        end if

        if (present(tropZ)) then
          tropZ(icol) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, icol, tropLev(icol), tP)
        end if
      end if
    end if

  end subroutine tropopause_twmo_column

  ! This routine implements the WMO definition of the tropopause (WMO, 1957; Seidel and Randel, 2006).
  ! Seidel, D. J., and W. J. Randel (2006),
  ! Variability and trends in the global tropopause estimated from radiosonde data,
  ! J. Geophys. Res., 111, D21101, doi:10.1029/2006JD007363.
  !
  ! This requires that the lapse rate be less than 2 K/km for an altitude range
  ! of 2 km. The search starts at the surface and stops the first time this
  ! criteria is met.
  subroutine tropopause_wmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                            tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    ! Local Variables
    real(kind_phys), parameter    :: ztrop_low   = 5000._kind_phys        ! lowest tropopause level allowed (m)
    real(kind_phys), parameter    :: ztrop_high  = 20000._kind_phys       ! highest tropopause level allowed (m)
    real(kind_phys), parameter    :: max_dtdz    = 0.002_kind_phys        ! max dt/dz for tropopause level (K/m)
    real(kind_phys), parameter    :: min_trop_dz = 2000._kind_phys        ! min tropopause thickness (m)

    integer                 :: i
    integer                 :: k
    integer                 :: k2
    real(kind_phys)                :: tP                           ! tropopause pressure (Pa)
    real(kind_phys)                :: dt

    ! Iterate over all of the columns.
    do i = 1, ncol

      ! Skip column in which the tropopause has already been found.
      if (tropLev(i) == NOTFOUND) then

        kloop: do k = pver-1, 2, -1

          ! Skip levels below the minimum and stop if nothing is found
          ! before the maximum.
          if (zm(i, k) < ztrop_low) then
            cycle kloop
          else if (zm(i, k) > ztrop_high) then
            exit kloop
          end if

          ! Compare the actual lapse rate to the threshold
          dt = t(i, k) - t(i, k-1)

          if (dt <= (max_dtdz * (zm(i, k-1) - zm(i, k)))) then

            ! Make sure that the lapse rate stays below the threshold for the
            ! specified range.
            k2loop: do k2 = k-1, 2, -1
              if ((zm(i, k2) - zm(i, k)) >= min_trop_dz) then
                tP = pmid(i, k)
                tropLev(i) = k
                exit k2loop
              end if

              dt = t(i, k) - t(i, k2)
              if (dt > (max_dtdz * (zm(i, k2) - zm(i, k)))) then
                exit k2loop
              end if
           end do k2loop

           if (tropLev(i) == NOTFOUND) then
              cycle kloop
           else

              ! Return the optional outputs
              if (present(tropP)) tropP(i) = tP

              if (present(tropT)) then
                tropT(i) = tropopause_interpolateT(pver, pmid, t, i, tropLev(i), tP)
              end if

              if (present(tropZ)) then
                tropZ(i) = tropopause_interpolateZ(pint, pmid, zi, zm, phis, i, tropLev(i), tP)
              end if

              exit kloop
            end if
          end if
        end do kloop
      end if
    end do

  end subroutine tropopause_wmo


  ! This routine searches for the cold point tropopause, and uses a parabolic
  ! fit of the coldest point and two adjacent points to interpolate the cold point
  ! between model levels.
  subroutine tropopause_cpp(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                            tropLev, tropP, tropT, tropZ)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    integer                 :: i

    ! Iterate over all of the columns.
    !$acc parallel loop gang vector &
    !$acc   copyin(pmid, t, zm) copy(tropLev, tropP, tropT, tropZ)
    do i = 1, ncol
      call tropopause_cpp_column(i, pver, pmid, t, zm, tropLev, tropP, tropT, tropZ)
    end do

  end subroutine tropopause_cpp

  ! tropopause_cpp for column icol.
  subroutine tropopause_cpp_column(icol, pver, pmid, t, zm, tropLev, tropP, tropT, tropZ)
    !$acc routine seq

    integer,         intent(in)         :: icol          ! column index
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    integer,            intent(inout)  :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)  :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)  :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)  :: tropZ(:)              ! tropopause height (m)

    ! Local Variables
    real(kind_phys), parameter    :: ztrop_low   = 5000._kind_phys        ! lowest tropopause level allowed (m)
    real(kind_phys), parameter    :: ztrop_high  = 25000._kind_phys       ! highest tropopause level allowed (m)

    integer                 :: k, firstk, lastk
    real(kind_phys)                :: tZ                           ! tropopause height (m)
    real(kind_phys)                :: tmin
    real(kind_phys)                :: f0, f1, f2
    real(kind_phys)                :: x0, x1, x2
    real(kind_phys)                :: c0, c1, c2
    real(kind_phys)                :: a, b, c

    firstk = 0
    lastk  = pver+1

    ! Skip column in which the tropopause has already been found.
    if (tropLev(icol) == NOTFOUND) then
      tmin = 1e6_kind_phys

      kloop: do k = pver-1, 2, -1

        ! Skip levels below the minimum and stop if nothing is found
        ! before the maximum.
        if (zm(icol, k) < ztrop_low) then
          firstk = k
          cycle kloop
        else if (zm(icol, k) > ztrop_high) then
          lastk = k
          exit kloop
        end if

        ! Find the coldest point
        if (t(icol, k) < tmin) then
          tropLev(icol) = k
          tmin = t(icol,k)
        end if
      end do kloop


      ! If the minimum is at the edge of the search range, then don't
      ! consider this to be a minima
      if ((tropLev(icol) >= (firstk-1)) .or. (tropLev(icol) <= (lastk+1))) then
        tropLev(icol) = NOTFOUND
      else

        ! If returning P, Z, or T, then do a parabolic fit using the
        ! cold point and it its 2 surrounding points to interpolate
        ! between model levels.
        if (present(tropP) .or. present(tropZ) .or. present(tropT)) then
          f0 = t(icol, tropLev(icol)-1)
          f1 = t(icol, tropLev(icol))
          f2 = t(icol, tropLev(icol)+1)

          x0 = zm(icol, tropLev(icol)-1)
          x1 = zm(icol, tropLev(icol))
          x2 = zm(icol, tropLev(icol)+1)

          c0 = (x0-x1)*(x0-x2)
          c1 = (x1-x0)*(x1-x2)
          c2 = (x2-x0)*(x2-x1)

          ! Determine the quadratic coefficients of:
          !   T = a * z^2 - b*z + c
          a = (f0/c0 + f1/c1 + f2/c2)
          b = (f0/c0*(x1+x2) + f1/c1*(x0+x2) + f2/c2*(x0+x1))
          c = f0/c0*x1*x2 + f1/c1*x0*x2 + f2/c2*x0*x1

          ! Find the altitude of the minimum temperature
          if (a == 0) then
             tropLev(icol) = NOTFOUND
          else
             tZ = 0.5_kind_phys * b / a

             ! The fit should be between the upper and lower points,
             ! so skip the point if the fit fails.
             if ((tZ >= x0) .or. (tZ <= x2)) then
               tropLev(icol) = NOTFOUND
             else

               ! Return the optional outputs
               if (present(tropP)) then
                 tropP(icol) = tropopause_interpolateP(pver, pmid, zm, icol, tropLev(icol), tZ)
               end if

               if (present(tropT)) then
                 tropT(icol) = a * tZ*tZ - b*tZ + c
               end if

               if (present(tropZ)) then
                 tropZ(icol) = tZ
               end if
             end if
          end if
        end if
      end if
    end if

  end subroutine tropopause_cpp_column

  ! Searches all the columns and attempts to identify the "chemical"
  ! tropopause. This is the lapse rate tropopause, backed up by the climatology
  ! if the lapse rate fails to find the tropopause at pressures higher than a certain
  ! threshold. This pressure threshold depends on latitude. Between 50S and 50N,
  ! the climatology is used if the lapse rate tropopause is not found at P > 75 hPa.
  ! At high latitude (poleward of 50), the threshold is increased to 125 hPa to
  ! eliminate false events that are sometimes detected in the cold polar stratosphere.
  !
  ! NOTE: This routine was adapted from code in chemistry.F90 and mo_gasphase_chemdr.F90.
  ! During the CCPP-ization, findChemTrop is now called from tropopause_find_run using method CHEMTROP
  ! and now also returns the standard tropLev, tropP, tropT, tropZ outputs (optional).
  ! The "backup" option is dropped as it is not used anywhere in current CAM.
  subroutine tropopause_findChemTrop(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                                     calday, tropp_p_loc, tropp_days, &
                                     tropLev, tropP, tropT, tropZ, errmsg, errflg)

    integer,         intent(in)         :: ncol          ! Number of atmospheric columns
    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: lat(:)        ! Latitudes (radians)
    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)

    real(kind_phys), intent(in)         :: calday        ! Day of year including fraction of day

    ! Climatological tropopause pressures (Pa), (ncol,ntimes=12).
    real(kind_phys), intent(in)         :: tropp_p_loc(:,:)
    real(kind_phys), intent(in)         :: tropp_days(:) ! Day-of-year for climo data, 12

    integer,                   intent(out)     :: tropLev(:)            ! tropopause level index
    real(kind_phys), optional, intent(inout)   :: tropP(:)              ! tropopause pressure (Pa)
    real(kind_phys), optional, intent(inout)   :: tropT(:)              ! tropopause temperature (K)
    real(kind_phys), optional, intent(inout)   :: tropZ(:)              ! tropopause height (m)

    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    ! Local Variable
    real(kind_phys)     :: dlats                ! latitude (degrees), this column
    integer             :: i

    errmsg = ' '
    errflg = 0

    ! First use the lapse rate tropopause.
    call tropopause_twmo(ncol, pver, lat, pint, pmid, t, zi, zm, phis, tropLev, tropP, tropT, tropZ)

    ! Now check high latitudes (poleward of 50) and set the level to the
    ! climatology if the level was not found or is at P <= 125 hPa.
    !$acc parallel loop gang vector private(dlats) copyin(lat, pmid) copy(tropLev)
    do i = 1, ncol
      dlats = lat(i) * cnst_rad2deg ! convert to degrees
      if (abs(dlats) > 50._kind_phys) then
        if (tropLev(i) .ne. NOTFOUND) then
          if (pmid(i, tropLev(i)) <= 12500._kind_phys) then
            tropLev(i) = NOTFOUND
          end if
        end if
      end if
    end do

    ! Now use the climatology backup
    if (any_notfound(tropLev)) then
      call tropopause_climate(ncol, pver, lat, pint, pmid, t, zi, zm, phis, &
                              calday, tropp_p_loc, tropp_days, &
                              tropLev, tropP, tropT, tropZ)
    end if

  end subroutine tropopause_findChemTrop

  ! This routine interpolates the pressures in the physics state to
  ! find the pressure at the specified tropopause altitude.
  function tropopause_interpolateP(pver, pmid, zm, icol, tropLev, tropZ)
    !$acc routine seq

    implicit none

    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    integer, intent(in)                 :: icol               ! column being processed
    integer, intent(in)                 :: tropLev            ! tropopause level index
    real(kind_phys), optional, intent(in)      :: tropZ              ! tropopause pressure (m)
    real(kind_phys)                            :: tropopause_interpolateP

    ! Local Variables
    real(kind_phys)   :: tropP              ! tropopause pressure (Pa)
    real(kind_phys)   :: dlogPdZ            ! dlog(p)/dZ

    ! Interpolate the temperature linearly against log(P)

    ! Is the tropopause at the midpoint?
    if (tropZ == zm(icol, tropLev)) then
      tropP = pmid(icol, tropLev)

    else if (tropZ > zm(icol, tropLev)) then

      ! It is above the midpoint? Make sure we aren't at the top.
      if (tropLev > 1) then
        dlogPdZ = (log(pmid(icol, tropLev)) - log(pmid(icol, tropLev - 1))) / &
          (zm(icol, tropLev) - zm(icol, tropLev - 1))
        tropP = pmid(icol, tropLev) + exp((tropZ - zm(icol, tropLev)) * dlogPdZ)
      end if
    else

      ! It is below the midpoint. Make sure we aren't at the bottom.
      if (tropLev < pver) then
        dlogPdZ =  (log(pmid(icol, tropLev + 1)) - log(pmid(icol, tropLev))) / &
          (zm(icol, tropLev + 1) - zm(icol, tropLev))
        tropP = pmid(icol, tropLev) + exp((tropZ - zm(icol, tropLev)) * dlogPdZ)
      end if
    end if

    tropopause_interpolateP = tropP
  end function tropopause_interpolateP


  ! This routine interpolates the temperatures in the physics state to
  ! find the temperature at the specified tropopause pressure.
  function tropopause_interpolateT(pver, pmid, t, icol, tropLev, tropP)
    !$acc routine seq

    implicit none

    integer,         intent(in)         :: pver          ! Number of vertical levels
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: t(:,:)        ! Temperature (K)
    integer, intent(in)                 :: icol               ! column being processed
    integer, intent(in)                 :: tropLev            ! tropopause level index
    real(kind_phys), optional, intent(in)      :: tropP              ! tropopause pressure (Pa)
    real(kind_phys)                            :: tropopause_interpolateT

    ! Local Variables
    real(kind_phys)   :: tropT              ! tropopause temperature (K)
    real(kind_phys)   :: dTdlogP            ! dT/dlog(P)

    ! Interpolate the temperature linearly against log(P)

    ! Is the tropopause at the midpoint?
    if (tropP == pmid(icol, tropLev)) then
      tropT = t(icol, tropLev)

    else if (tropP < pmid(icol, tropLev)) then

      ! It is above the midpoint? Make sure we aren't at the top.
      if (tropLev > 1) then
        dTdlogP = (t(icol, tropLev) - t(icol, tropLev - 1)) / &
          (log(pmid(icol, tropLev)) - log(pmid(icol, tropLev - 1)))
        tropT = t(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dTdlogP
      end if
    else

      ! It is below the midpoint. Make sure we aren't at the bottom.
      if (tropLev < pver) then
        dTdlogP = (t(icol, tropLev + 1) - t(icol, tropLev)) / &
          (log(pmid(icol, tropLev + 1)) - log(pmid(icol, tropLev)))
        tropT = t(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dTdlogP
      end if
    end if

    tropopause_interpolateT = tropT
  end function tropopause_interpolateT


  ! This routine interpolates the geopotential height in the physics state to
  ! find the geopotential height at the specified tropopause pressure.
  function tropopause_interpolateZ(pint, pmid, zi, zm, phis, icol, tropLev, tropP)
    !$acc routine seq

    real(kind_phys), intent(in)         :: pint(:,:)     ! Interface pressures (Pa)
    real(kind_phys), intent(in)         :: pmid(:,:)     ! Midpoint pressures (Pa)
    real(kind_phys), intent(in)         :: zi(:,:)       ! Geopotential height above surface at interfaces (m)
    real(kind_phys), intent(in)         :: zm(:,:)       ! Geopotential height above surface at midpoints (m)
    real(kind_phys), intent(in)         :: phis(:)       ! Surface geopotential (m2 s-2)
    integer, intent(in)                 :: icol               ! column being processed
    integer, intent(in)                 :: tropLev            ! tropopause level index
    real(kind_phys), optional, intent(in)      :: tropP              ! tropopause pressure (Pa)
    real(kind_phys)                            :: tropopause_interpolateZ

    ! Local Variables
    real(kind_phys)   :: tropZ              ! tropopause geopotential height (m)
    real(kind_phys)   :: dZdlogP            ! dZ/dlog(P)

    ! Interpolate the geopotential height linearly against log(P)

    ! Is the tropopause at the midpoint?
    if (tropP == pmid(icol, tropLev)) then
      tropZ = zm(icol, tropLev)

    else if (tropP < pmid(icol, tropLev)) then

      ! It is above the midpoint? Make sure we aren't at the top.
      dZdlogP = (zm(icol, tropLev) - zi(icol, tropLev)) / &
        (log(pmid(icol, tropLev)) - log(pint(icol, tropLev)))
      tropZ = zm(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dZdlogP
    else

      ! It is below the midpoint. Make sure we aren't at the bottom.
      dZdlogP = (zm(icol, tropLev) - zi(icol, tropLev+1)) / &
        (log(pmid(icol, tropLev)) - log(pint(icol, tropLev+1)))
      tropZ = zm(icol, tropLev) + (log(tropP) - log(pmid(icol, tropLev))) * dZdlogP
    end if

    tropopause_interpolateZ = tropZ + phis(icol)*cnst_rga
  end function tropopause_interpolateZ
end module tropopause_find
