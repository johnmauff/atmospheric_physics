! Cloud cover diagnostics
! total/low/mid/high cloud cover assuming maximum-random overlap and 3D cloud fraction.
module cloud_cover_diagnostics
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: cloud_cover_diagnostics_init
   public :: cloud_cover_diagnostics_run

   ! Pressure bounds of the low/mid/high cloud cover ranges [Pa]
   real(kind_phys), parameter :: plowmax = 120000._kind_phys ! Max pressure for low cloud cover range
   real(kind_phys), parameter :: plowmin =  70000._kind_phys ! Min pressure for low cloud cover range
   real(kind_phys), parameter :: pmedmax =  70000._kind_phys ! Max pressure for mid cloud cover range
   real(kind_phys), parameter :: pmedmin =  40000._kind_phys ! Min pressure for mid cloud cover range
   real(kind_phys), parameter :: phghmax =  40000._kind_phys ! Max pressure for hgh cloud cover range
   real(kind_phys), parameter :: phghmin =   5000._kind_phys ! Min pressure for hgh cloud cover range

contains

   !> \section arg_table_cloud_cover_diagnostics_init  Argument Table
   !! \htmlinclude cloud_cover_diagnostics_init.html
   subroutine cloud_cover_diagnostics_init(errmsg, errflg)
      use cam_history,         only: history_add_field
      use cam_history_support, only: horiz_only

      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      character(len=128) :: long_name
      character(len=*), parameter :: range_fmt = '(a, f7.0, a, f7.0, a)'

      errmsg = ''
      errflg = 0

      call history_add_field('CLOUD',  'Cloud fraction', 'lev', 'avg', 'fraction')
      call history_add_field('CLDTOT', 'Vertically-integrated total cloud', horiz_only, 'avg', 'fraction')

      ! Long names carry the pressure ranges from the module parameters so the history file stays self-describing.
      write(long_name, range_fmt) 'Vertically-integrated low cloud from ', plowmin, ' to ', plowmax, ' Pa'
      call history_add_field('CLDLOW', trim(long_name), horiz_only, 'avg', 'fraction')
      write(long_name, range_fmt) 'Vertically-integrated mid-level cloud from ', pmedmin, ' to ', pmedmax, ' Pa'
      call history_add_field('CLDMED', trim(long_name), horiz_only, 'avg', 'fraction')
      write(long_name, range_fmt) 'Vertically-integrated high cloud from ', phghmin, ' to ', phghmax, ' Pa'
      call history_add_field('CLDHGH', trim(long_name), horiz_only, 'avg', 'fraction')

   end subroutine cloud_cover_diagnostics_init

   !> \section arg_table_cloud_cover_diagnostics_run  Argument Table
   !! \htmlinclude cloud_cover_diagnostics_run.html
   subroutine cloud_cover_diagnostics_run(ncol, pver, pverp, pint, pmid, cld, errmsg, errflg)
      use cam_history,              only: history_out_field
      use cloud_optical_properties, only: cldovrlap

      integer,            intent(in)  :: ncol
      integer,            intent(in)  :: pver
      integer,            intent(in)  :: pverp
      real(kind_phys),    intent(in)  :: pint(:,:)  ! Air pressure at interfaces [Pa]
      real(kind_phys),    intent(in)  :: pmid(:,:)  ! Air pressure at midpoints [Pa]
      real(kind_phys),    intent(in)  :: cld(:,:)   ! Cloud fraction [fraction]
      character(len=*),   intent(out) :: errmsg
      integer,            intent(out) :: errflg

      ! Local variables
      integer         :: nmxrgn(ncol)        ! Number of maximally overlapped regions
      real(kind_phys) :: pmxrgn(ncol,pverp)  ! Maximum pressure of each maximally overlapped region
      real(kind_phys) :: cldtot(ncol)        ! Total random overlap cloud cover
      real(kind_phys) :: cldlow(ncol)        ! Low random overlap cloud cover
      real(kind_phys) :: cldmed(ncol)        ! Middle random overlap cloud cover
      real(kind_phys) :: cldhgh(ncol)        ! High random overlap cloud cover

      errmsg = ''
      errflg = 0

      ! Partition each column into maximally-overlapped regions of
      ! contiguous cloudy layers, then collapse to cover diagnostics.
      !$acc data present(pint, pmid, cld) &
      !$acc      create(nmxrgn, pmxrgn, cldtot, cldlow, cldmed, cldhgh)
      call cldovrlap(ncol, pver, pverp, pint, cld, nmxrgn, pmxrgn)
      call cldsav(ncol, pver, cld, pmid, cldtot, cldlow, cldmed, cldhgh, nmxrgn, pmxrgn)

      call history_out_field('CLOUD',  cld)
      call history_out_field('CLDTOT', cldtot)
      call history_out_field('CLDLOW', cldlow)
      call history_out_field('CLDMED', cldmed)
      call history_out_field('CLDHGH', cldhgh)
      !$acc end data

   end subroutine cloud_cover_diagnostics_run

   ! Compute total & 3 levels of cloud fraction assuming maximum-random overlap.
   ! Pressure ranges for the 3 cloud levels are specified.
   !
   ! Note: pmxrgn from cldovrlap uses the CAMRT-lineage factor of 10 (dyn cm-2)
   ! relative to pmid [Pa]; kept as is for bit-for-bit parity with CAM.
   !
   ! Runs on the device in OpenACC builds (inputs and outputs present).
   ! Each column is independent, so the per-column kernels below do the same
   ! operations in the same order as the original (ityp, k, i) loops.
   !
   ! Original author: W. Collins
   subroutine cldsav(ncol, pver, cld, pmid, cldtot, cldlow, cldmed, cldhgh, nmxrgn, pmxrgn)

      ! Input arguments
      integer,         intent(in) :: ncol           ! number of atmospheric columns
      integer,         intent(in) :: pver           ! number of vertical layers
      real(kind_phys), intent(in) :: cld(:,:)       ! Cloud fraction
      real(kind_phys), intent(in) :: pmid(:,:)      ! Level pressures
      real(kind_phys), intent(in) :: pmxrgn(:,:)    ! Maximum values of pressure for each
      !    maximally overlapped region.
      !    0->pmxrgn(i,1) is range of pressure for
      !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
      !    2nd region, etc
      integer,         intent(in) :: nmxrgn(:)      ! # of maximally overlapped regions

      ! Output arguments
      real(kind_phys), intent(out) :: cldtot(:)     ! Total random overlap cloud cover
      real(kind_phys), intent(out) :: cldlow(:)     ! Low random overlap cloud cover
      real(kind_phys), intent(out) :: cldmed(:)     ! Middle random overlap cloud cover
      real(kind_phys), intent(out) :: cldhgh(:)     ! High random overlap cloud cover

      ! Local workspace
      integer  :: i,k                    ! Longitude,level indices
      integer  :: irgn                   ! Max-overlap region index
      integer  :: max_nmxrgn             ! maximum value of nmxrgn over columns
      integer  :: ityp                   ! Type counter
      real(kind_phys) :: clrsky          ! Max-random clear sky fraction
      real(kind_phys) :: clrskymax       ! Maximum overlap clear sky fraction
      real(kind_phys) :: pmin, pmax      ! Pressure range of this cloud type
      real(kind_phys) :: ptypmin(4)
      real(kind_phys) :: ptypmax(4)

      ptypmin = [phghmin, plowmin, pmedmin, phghmin]
      ptypmax = [plowmax, plowmax, pmedmax, phghmax]

      !$acc data present(cld, pmid, pmxrgn, nmxrgn, cldtot, cldlow, cldmed, cldhgh)

      ! Initialize region number
      max_nmxrgn = -1
      !$acc parallel loop gang vector reduction(max:max_nmxrgn)
      do i=1,ncol
         max_nmxrgn = max(max_nmxrgn,nmxrgn(i))
      end do

      do ityp = 1, 4
         pmin = ptypmin(ityp)
         pmax = ptypmax(ityp)
         !$acc parallel loop gang vector private(irgn, clrsky, clrskymax)
         do i=1,ncol
            irgn = 1
            !$acc loop seq
            do k =1,max_nmxrgn-1
               if (pmxrgn(i,irgn) < pmin .and. irgn < nmxrgn(i)) then
                  irgn = irgn + 1
               end if
            end do

            ! Compute cloud amount by estimating clear-sky amounts
            clrsky    = 1.0_kind_phys
            clrskymax = 1.0_kind_phys
            !$acc loop seq
            do k = 1, pver
               if (pmid(i,k) >= pmin .and. pmid(i,k) <= pmax) then
                  if (pmxrgn(i,irgn) < pmid(i,k) .and. irgn < nmxrgn(i)) then
                     irgn = irgn + 1
                     clrsky = clrsky * clrskymax
                     clrskymax = 1.0_kind_phys
                  end if
                  clrskymax = min(clrskymax,1.0_kind_phys-cld(i,k))
               end if
            end do
            if (ityp == 1) cldtot(i) = 1.0_kind_phys - (clrsky * clrskymax)
            if (ityp == 2) cldlow(i) = 1.0_kind_phys - (clrsky * clrskymax)
            if (ityp == 3) cldmed(i) = 1.0_kind_phys - (clrsky * clrskymax)
            if (ityp == 4) cldhgh(i) = 1.0_kind_phys - (clrsky * clrskymax)
         end do
      end do

      !$acc end data

   end subroutine cldsav

end module cloud_cover_diagnostics
