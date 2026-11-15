module atmos_phys_interp_utils
   ! Linear interpolation utilities used by the physics schemes.
   !
   ! Taken from the host models' interpolate_data module (CAM, CAM-SIMA:
   ! "Gathered from a number of places and put into the current format by
   ! Jim Edwards"), only the routines the physics uses, so that schemes do
   ! not depend on a host module. The arithmetic is unchanged. Differences
   ! from interpolate_data:
   !   * errors return through errmsg/errflg instead of stopping the model
   !     (lininterp_init, the one-call lininterp and vertinterp);
   !   * vertinterp takes the gas constant and 1/g as arguments, needed only
   !     when extrapolating, instead of using the host's physconst;
   !   * kinds are kind_phys.
   ! A new module name, as the hosts keep their own interpolate_data.

   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: interp_type
   public :: lininterp_init
   public :: lininterp
   public :: lininterp_finish
   public :: vertinterp

   type interp_type
      real(kind_phys), pointer :: wgts(:)
      real(kind_phys), pointer :: wgtn(:)
      integer, pointer         :: jjm(:)
      integer, pointer         :: jjp(:)
   end type interp_type

   interface lininterp
      module procedure lininterp_full1d
      module procedure lininterp1d
      module procedure lininterp2d2d
      module procedure lininterp2d1d
      module procedure lininterp3d2d
   end interface lininterp

   integer, parameter, public :: extrap_method_zero = 0
   integer, parameter, public :: extrap_method_bndry = 1
   integer, parameter, public :: extrap_method_cycle = 2

contains

   subroutine lininterp_full1d(arrin, yin, nin, arrout, yout, nout, errmsg, errflg)
      ! Interpolate arrin, on mesh yin, to mesh yout in one call (weights
      ! computed, used and freed here; boundary values beyond yin's range).
      integer,          intent(in)  :: nin, nout
      real(kind_phys),  intent(in)  :: arrin(nin), yin(nin), yout(nout)
      real(kind_phys),  intent(out) :: arrout(nout)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      type(interp_type) :: interp_wgts

      call lininterp_init(yin, nin, yout, nout, extrap_method_bndry, interp_wgts, &
                          errmsg=errmsg, errflg=errflg)
      if (errflg /= 0) return
      call lininterp1d(arrin, nin, arrout, nout, interp_wgts)
      call lininterp_finish(interp_wgts)

   end subroutine lininterp_full1d

   subroutine lininterp_init(yin, nin, yout, nout, extrap_method, interp_wgts, &
                             cyclicmin, cyclicmax, errmsg, errflg)
      !
      ! Description:
      !   Initialize a variable of type(interp_type) with weights for linear interpolation.
      !       this variable can then be used in calls to lininterp1d and lininterp2d.
      !   yin is a 1d array of length nin of locations to interpolate from - this array must
      !       be monotonic but can be increasing or decreasing
      !   yout is a 1d array of length nout of locations to interpolate to, this array need
      !       not be ordered
      !   extrap_method determines how to handle yout points beyond the bounds of yin
      !       if 0 set values outside output grid to 0
      !       if 1 set to boundary value
      !       if 2 set to cyclic boundaries
      !         optional values cyclicmin and cyclicmax can be used to set the bounds of the
      !         cyclic mapping - these default to 0 and 360.
      !   On an error (errflg /= 0), interp_wgts is not usable.
      !

      integer,         intent(in) :: nin
      integer,         intent(in) :: nout
      real(kind_phys), intent(in) :: yin(:)           ! input mesh
      real(kind_phys), intent(in) :: yout(:)          ! output mesh
      integer,         intent(in) :: extrap_method    ! if 0 set values outside output grid to 0
                                                      ! if 1 set to boundary value
                                                      ! if 2 set to cyclic boundaries
      type(interp_type), intent(out) :: interp_wgts
      real(kind_phys), intent(in), optional :: cyclicmin, cyclicmax
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      real(kind_phys) :: cmin, cmax
      real(kind_phys) :: extrap
      real(kind_phys) :: dyinwrap
      real(kind_phys) :: ratio
      real(kind_phys) :: avgdyin
      integer :: j, icount
      integer :: jj
      real(kind_phys), pointer :: wgts(:)
      real(kind_phys), pointer :: wgtn(:)
      integer, pointer :: jjm(:)
      integer, pointer :: jjp(:)
      logical :: increasing

      errmsg = ''
      errflg = 0
      !
      ! Check validity of input coordinate arrays: must be monotonically increasing,
      ! and have a total of at least 2 elements
      !
      if (nin < 2) then
         errmsg = 'LININTERP: Must have at least 2 input points for interpolation'
         errflg = 1
         return
      end if
      if (present(cyclicmin)) then
         cmin = cyclicmin
      else
         cmin = 0._kind_phys
      end if
      if (present(cyclicmax)) then
         cmax = cyclicmax
      else
         cmax = 360._kind_phys
      end if
      if (cmax <= cmin) then
         errmsg = 'LININTERP: cyclic min value must be < max value'
         errflg = 1
         return
      end if
      increasing = .true.
      icount = 0
      do j = 1, nin-1
         if (yin(j) > yin(j+1)) icount = icount + 1
      end do
      if (icount == nin-1) then
         increasing = .false.
         icount = 0
      end if
      if (icount > 0) then
         errmsg = 'LININTERP: Non-monotonic input coordinate array found'
         errflg = 1
         return
      end if
      allocate(interp_wgts%jjm(nout), &
               interp_wgts%jjp(nout), &
               interp_wgts%wgts(nout), &
               interp_wgts%wgtn(nout))

      jjm => interp_wgts%jjm
      jjp => interp_wgts%jjp
      wgts => interp_wgts%wgts
      wgtn => interp_wgts%wgtn

      !
      ! Initialize index arrays for later checking
      !
      jjm = 0
      jjp = 0

      extrap = 0._kind_phys
      select case (extrap_method)
      case (extrap_method_zero)
         !
         ! For values which extend beyond boundaries, set weights
         ! such that values will be 0.
         !
         do j = 1, nout
            if (increasing) then
               if (yout(j) < yin(1)) then
                  jjm(j) = 1
                  jjp(j) = 1
                  wgts(j) = 0._kind_phys
                  wgtn(j) = 0._kind_phys
                  extrap = extrap + 1._kind_phys
               else if (yout(j) > yin(nin)) then
                  jjm(j) = nin
                  jjp(j) = nin
                  wgts(j) = 0._kind_phys
                  wgtn(j) = 0._kind_phys
                  extrap = extrap + 1._kind_phys
               end if
            else
               if (yout(j) > yin(1)) then
                  jjm(j) = 1
                  jjp(j) = 1
                  wgts(j) = 0._kind_phys
                  wgtn(j) = 0._kind_phys
                  extrap = extrap + 1._kind_phys
               else if (yout(j) < yin(nin)) then
                  jjm(j) = nin
                  jjp(j) = nin
                  wgts(j) = 0._kind_phys
                  wgtn(j) = 0._kind_phys
                  extrap = extrap + 1._kind_phys
               end if
            end if
         end do
      case (extrap_method_bndry)
         !
         ! For values which extend beyond boundaries, set weights
         ! such that values will just be copied.
         !
         do j = 1, nout
            if (increasing) then
               if (yout(j) <= yin(1)) then
                  jjm(j) = 1
                  jjp(j) = 1
                  wgts(j) = 1._kind_phys
                  wgtn(j) = 0._kind_phys
                  extrap = extrap + 1._kind_phys
               else if (yout(j) > yin(nin)) then
                  jjm(j) = nin
                  jjp(j) = nin
                  wgts(j) = 1._kind_phys
                  wgtn(j) = 0._kind_phys
                  extrap = extrap + 1._kind_phys
               end if
            else
               if (yout(j) > yin(1)) then
                  jjm(j) = 1
                  jjp(j) = 1
                  wgts(j) = 1._kind_phys
                  wgtn(j) = 0._kind_phys
                  extrap = extrap + 1._kind_phys
               else if (yout(j) <= yin(nin)) then
                  jjm(j) = nin
                  jjp(j) = nin
                  wgts(j) = 1._kind_phys
                  wgtn(j) = 0._kind_phys
                  extrap = extrap + 1._kind_phys
               end if
            end if
         end do
      case (extrap_method_cycle)
         !
         ! For values which extend beyond boundaries, set weights
         ! for circular boundaries
         !
         dyinwrap = yin(1) + (cmax-cmin) - yin(nin)
         avgdyin = abs(yin(nin)-yin(1))/(nin-1._kind_phys)
         ratio = dyinwrap/avgdyin
         if (ratio < 0.9_kind_phys .or. ratio > 1.1_kind_phys) then
            write(errmsg, '(a,es12.5,a,es12.5)') 'LININTERP: Bad dyinwrap value = ', &
                 dyinwrap, ' avg = ', avgdyin
            errflg = 1
            return
         end if

         do j = 1, nout
            if (increasing) then
               if (yout(j) <= yin(1)) then
                  jjm(j) = nin
                  jjp(j) = 1
                  wgts(j) = (yin(1)-yout(j))/dyinwrap
                  wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
               else if (yout(j) > yin(nin)) then
                  jjm(j) = nin
                  jjp(j) = 1
                  wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
                  wgtn(j) = (yout(j)-yin(nin))/dyinwrap
               end if
            else
               if (yout(j) > yin(1)) then
                  jjm(j) = nin
                  jjp(j) = 1
                  wgts(j) = (yin(1)-yout(j))/dyinwrap
                  wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
               else if (yout(j) <= yin(nin)) then
                  jjm(j) = nin
                  jjp(j) = 1
                  wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
                  wgtn(j) = (yout(j)+(cmax-cmin)-yin(nin))/dyinwrap
               end if

            end if
         end do
      end select

      !
      ! Loop though output indices finding input indices and weights
      !
      if (increasing) then
         do j = 1, nout
            do jj = 1, nin-1
               if (yout(j) > yin(jj) .and. yout(j) <= yin(jj+1)) then
                  jjm(j) = jj
                  jjp(j) = jj + 1
                  wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                  wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                  exit
               end if
            end do
         end do
      else
         do j = 1, nout
            do jj = 1, nin-1
               if (yout(j) <= yin(jj) .and. yout(j) > yin(jj+1)) then
                  jjm(j) = jj
                  jjp(j) = jj + 1
                  wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                  wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                  exit
               end if
            end do
         end do
      end if

      !
      ! Check that interp/extrap points have been found for all outputs
      !
      icount = 0
      do j = 1, nout
         if (jjm(j) == 0 .or. jjp(j) == 0) icount = icount + 1
         ratio = wgts(j) + wgtn(j)
         if ((ratio < 0.9_kind_phys .or. ratio > 1.1_kind_phys) .and. extrap_method /= 0) then
            write(errmsg, '(a,i0,a,2es12.5,a,2(1x,i0),a,l1,a,i0)') &
                 'LININTERP: Bad weight computed in lininterp_init: point ', j, &
                 ' weights', wgts(j), wgtn(j), ' indices', jjm(j), jjp(j), &
                 ' increasing ', increasing, ' extrap_method ', extrap_method
            errflg = 1
            return
         end if
      end do
      if (icount > 0) then
         errmsg = 'LININTERP: Point found without interp indices'
         errflg = 1
         return
      end if

   end subroutine lininterp_init

   subroutine lininterp1d(arrin, n1, arrout, m1, interp_wgts)
      !-----------------------------------------------------------------------
      !
      ! Purpose: Do a linear interpolation from input mesh to output
      !          mesh with weights as set in lininterp_init.
      !
      !
      ! Author: Jim Edwards
      !
      !-----------------------------------------------------------------------
      !
      ! Arguments
      !
      integer, intent(in) :: n1                 ! number of input latitudes
      integer, intent(in) :: m1                 ! number of output latitudes

      real(kind_phys),   intent(in)  :: arrin(n1)    ! input array of values to interpolate
      type(interp_type), intent(in)  :: interp_wgts
      real(kind_phys),   intent(out) :: arrout(m1)   ! interpolated array

      !
      ! Local workspace
      !
      integer :: j                ! latitude indices
      integer, pointer :: jjm(:)
      integer, pointer :: jjp(:)

      real(kind_phys), pointer :: wgts(:)
      real(kind_phys), pointer :: wgtn(:)

      jjm => interp_wgts%jjm
      jjp => interp_wgts%jjp
      wgts => interp_wgts%wgts
      wgtn => interp_wgts%wgtn

      !
      ! Do the interpolation
      !
      do j = 1, m1
         arrout(j) = arrin(jjm(j))*wgts(j) + arrin(jjp(j))*wgtn(j)
      end do

   end subroutine lininterp1d

   subroutine lininterp2d2d(arrin, n1, n2, arrout, m1, m2, wgt1, wgt2)
      !-----------------------------------------------------------------------
      !
      ! Arguments
      !
      integer,           intent(in)  :: n1, n2, m1, m2
      real(kind_phys),   intent(in)  :: arrin(n1,n2)    ! input array of values to interpolate
      type(interp_type), intent(in)  :: wgt1, wgt2
      real(kind_phys),   intent(out) :: arrout(m1,m2)   ! interpolated array
      !
      ! locals
      !
      integer :: i, j                ! indices
      integer, pointer :: iim(:), jjm(:)
      integer, pointer :: iip(:), jjp(:)

      real(kind_phys), pointer :: wgts1(:), wgts2(:)
      real(kind_phys), pointer :: wgtn1(:), wgtn2(:)

      real(kind_phys) :: arrtmp(n1,m2)

      jjm => wgt2%jjm
      jjp => wgt2%jjp
      wgts2 => wgt2%wgts
      wgtn2 => wgt2%wgtn

      iim => wgt1%jjm
      iip => wgt1%jjp
      wgts1 => wgt1%wgts
      wgtn1 => wgt1%wgtn

      do j = 1, m2
         do i = 1, n1
            arrtmp(i,j) = arrin(i,jjm(j))*wgts2(j) + arrin(i,jjp(j))*wgtn2(j)
         end do
      end do

      do j = 1, m2
         do i = 1, m1
            arrout(i,j) = arrtmp(iim(i),j)*wgts1(i) + arrtmp(iip(i),j)*wgtn1(i)
         end do
      end do

   end subroutine lininterp2d2d

   subroutine lininterp2d1d(arrin, n1, n2, arrout, m1, wgt1, wgt2, fldname)
      !-----------------------------------------------------------------------
      !
      ! Arguments
      !
      integer,           intent(in)  :: n1, n2, m1
      real(kind_phys),   intent(in)  :: arrin(n1,n2)    ! input array of values to interpolate
      type(interp_type), intent(in)  :: wgt1, wgt2
      real(kind_phys),   intent(out) :: arrout(m1)      ! interpolated array
      character(len=*),  intent(in), optional :: fldname(:)
      !
      ! locals
      !
      integer :: i                ! indices
      integer, pointer :: iim(:), jjm(:)
      integer, pointer :: iip(:), jjp(:)

      real(kind_phys), pointer :: wgts(:), wgte(:)
      real(kind_phys), pointer :: wgtn(:), wgtw(:)

      jjm => wgt2%jjm
      jjp => wgt2%jjp
      wgts => wgt2%wgts
      wgtn => wgt2%wgtn

      iim => wgt1%jjm
      iip => wgt1%jjp
      wgtw => wgt1%wgts
      wgte => wgt1%wgtn

      do i = 1, m1
         arrout(i) = arrin(iim(i),jjm(i))*wgtw(i)*wgts(i)+arrin(iip(i),jjm(i))*wgte(i)*wgts(i) + &
                     arrin(iim(i),jjp(i))*wgtw(i)*wgtn(i)+arrin(iip(i),jjp(i))*wgte(i)*wgtn(i)
      end do

   end subroutine lininterp2d1d

   subroutine lininterp3d2d(arrin, n1, n2, n3, arrout, m1, len1, wgt1, wgt2)
      !-----------------------------------------------------------------------
      !
      ! Arguments
      !
      integer,           intent(in)  :: n1, n2, n3, m1, len1   ! m1 is to len1 as ncols is to pcols
      real(kind_phys),   intent(in)  :: arrin(n1,n2,n3)        ! input array of values to interpolate
      type(interp_type), intent(in)  :: wgt1, wgt2
      real(kind_phys),   intent(out) :: arrout(len1, n3)       ! interpolated array

      !
      ! locals
      !
      integer :: i, k               ! indices
      integer, pointer :: iim(:), jjm(:)
      integer, pointer :: iip(:), jjp(:)

      real(kind_phys), pointer :: wgts(:), wgte(:)
      real(kind_phys), pointer :: wgtn(:), wgtw(:)

      jjm => wgt2%jjm
      jjp => wgt2%jjp
      wgts => wgt2%wgts
      wgtn => wgt2%wgtn

      iim => wgt1%jjm
      iip => wgt1%jjp
      wgtw => wgt1%wgts
      wgte => wgt1%wgtn

      do k = 1, n3
         do i = 1, m1
            arrout(i,k) = arrin(iim(i),jjm(i),k)*wgtw(i)*wgts(i)+arrin(iip(i),jjm(i),k)*wgte(i)*wgts(i) + &
                          arrin(iim(i),jjp(i),k)*wgtw(i)*wgtn(i)+arrin(iip(i),jjp(i),k)*wgte(i)*wgtn(i)
         end do
      end do

   end subroutine lininterp3d2d

   subroutine lininterp_finish(interp_wgts)
      type(interp_type), intent(inout) :: interp_wgts

      deallocate(interp_wgts%jjm, &
                 interp_wgts%jjp, &
                 interp_wgts%wgts, &
                 interp_wgts%wgtn)

      nullify(interp_wgts%jjm, &
              interp_wgts%jjp, &
              interp_wgts%wgts, &
              interp_wgts%wgtn)

   end subroutine lininterp_finish

   !=========================================================================================

   subroutine vertinterp(ncol, ncold, nlev, pin, pout, arrin, arrout, errmsg, errflg, &
                         extrapolate, ln_interp, ps, phis, tbot, rair, rga)

      ! Vertically interpolate input array to output pressure level.  The
      ! interpolation is linear in either pressure (the default) or ln pressure.
      !
      ! If above the top boundary then return top boundary value.
      !
      ! If below bottom boundary then use bottom boundary value, or optionally
      ! for T or Z use the extrapolation algorithm from ECMWF (which was taken
      ! from the NCL code base). Extrapolation also needs rair (dry air gas
      ! constant, J K-1 kg-1) and rga (reciprocal of gravity, s2 m-1).

      !------------------------------Arguments--------------------------------
      integer,          intent(in)  :: ncol              ! number active columns
      integer,          intent(in)  :: ncold             ! column dimension
      integer,          intent(in)  :: nlev              ! vertical dimension
      real(kind_phys),  intent(in)  :: pin(ncold,nlev)   ! input pressure levels
      real(kind_phys),  intent(in)  :: pout              ! output pressure level
      real(kind_phys),  intent(in)  :: arrin(ncold,nlev) ! input  array
      real(kind_phys),  intent(out) :: arrout(ncold)     ! output array (interpolated)
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      character(len=1), optional, intent(in) :: extrapolate ! either 'T' or 'Z'
      logical,          optional, intent(in) :: ln_interp   ! set true to interolate
                                                            ! in ln(pressure)
      real(kind_phys),  optional, intent(in) :: ps(ncold)   ! surface pressure
      real(kind_phys),  optional, intent(in) :: phis(ncold) ! surface geopotential
      real(kind_phys),  optional, intent(in) :: tbot(ncold) ! temperature at bottom level
      real(kind_phys),  optional, intent(in) :: rair        ! dry air gas constant
      real(kind_phys),  optional, intent(in) :: rga         ! reciprocal of gravity

      !---------------------------Local variables-----------------------------
      real(kind_phys) :: alpha
      logical  :: linear_interp
      logical  :: do_extrapolate_T, do_extrapolate_Z

      integer  :: i, k                    ! indices
      integer  :: kupper(ncold)           ! Level indices for interpolation
      real(kind_phys) :: dpu              ! upper level pressure difference
      real(kind_phys) :: dpl              ! lower level pressure difference
      logical  :: found(ncold)            ! true if input levels found
      logical  :: error                   ! error flag
      !----------------------------------------------------------------------------

      errmsg = ''
      errflg = 0

      do_extrapolate_T = .false.
      do_extrapolate_Z = .false.
      if (present(extrapolate)) then

         if (extrapolate == 'T') do_extrapolate_T = .true.
         if (extrapolate == 'Z') do_extrapolate_Z = .true.

         if (.not. do_extrapolate_T .and. .not. do_extrapolate_Z) then
            errmsg = 'VERTINTERP: extrapolate must be T or Z'
            errflg = 1
            return
         end if
         if (.not. present(ps)) then
            errmsg = 'VERTINTERP: ps required for extrapolation'
            errflg = 1
            return
         end if
         if (.not. present(phis)) then
            errmsg = 'VERTINTERP: phis required for extrapolation'
            errflg = 1
            return
         end if
         if (do_extrapolate_Z) then
            if (.not. present(tbot)) then
               errmsg = 'VERTINTERP: tbot required for extrapolating Z'
               errflg = 1
               return
            end if
         end if
         if (.not. (present(rair) .and. present(rga))) then
            errmsg = 'VERTINTERP: rair and rga required for extrapolation'
            errflg = 1
            return
         end if
         alpha = 0.0065_kind_phys*rair*rga
      end if

      linear_interp = .true.
      if (present(ln_interp)) then
         if (ln_interp) linear_interp = .false.
      end if

      do i = 1, ncol
         found(i)  = .false.
         kupper(i) = 1
      end do
      error = .false.

      ! Find indices of upper pressure bound
      do k = 1, nlev - 1
         do i = 1, ncol
            if ((.not. found(i)) .and. pin(i,k)<pout .and. pout<=pin(i,k+1)) then
               found(i)  = .true.
               kupper(i) = k
            end if
         end do
      end do

      do i = 1, ncol

         if (pout <= pin(i,1)) then

            ! if above top pressure level then use value at top (no interpolation)
            arrout(i) = arrin(i,1)

         else if (pout >= pin(i,nlev)) then

            if (do_extrapolate_T) then
               ! use ECMWF algorithm to extrapolate T
               arrout(i) = extrapolate_T()

            else if (do_extrapolate_Z) then
               ! use ECMWF algorithm to extrapolate Z
               arrout(i) = extrapolate_Z()

            else
               ! use bottom value
               arrout(i) = arrin(i,nlev)
            end if

         else if (found(i)) then
            if (linear_interp) then
               ! linear interpolation in p
               dpu = pout - pin(i,kupper(i))
               dpl = pin(i,kupper(i)+1) - pout
               arrout(i) = (arrin(i,kupper(i))*dpl + arrin(i,kupper(i)+1)*dpu)/(dpl + dpu)
            else
               ! linear interpolation in ln(p)
               arrout(i) = arrin(i,kupper(i)) + (arrin(i,kupper(i)+1) - arrin(i,kupper(i)))* &
                                                log(pout/pin(i,kupper(i))) /                 &
                                                log(pin(i,kupper(i)+1)/pin(i,kupper(i)))
            end if
         else
            error = .true.
         end if
      end do

      ! Error check
      if (error) then
         errmsg = 'VERTINTERP: output pressure not bracketed by the input levels'
         errflg = 1
      end if

   contains

      !======================================================================================

      real(kind_phys) function extrapolate_T()

         ! local variables
         real(kind_phys) :: sixth
         real(kind_phys) :: tstar
         real(kind_phys) :: hgt
         real(kind_phys) :: alnp
         real(kind_phys) :: t0
         real(kind_phys) :: tplat
         real(kind_phys) :: tprime0
         !-------------------------------------------------------------------------

         sixth  = 1._kind_phys / 6._kind_phys
         tstar = arrin(i,nlev)*(1._kind_phys + alpha*(ps(i)/pin(i,nlev) - 1._kind_phys))
         hgt   = phis(i)*rga

         if (hgt < 2000._kind_phys) then
            alnp = alpha*log(pout/ps(i))
         else
            t0 = tstar + 0.0065_kind_phys*hgt
            tplat = min(t0, 298._kind_phys)

            if (hgt <= 2500._kind_phys) then
               tprime0 = 0.002_kind_phys*((2500._kind_phys - hgt)*t0 + &
                                          (hgt - 2000._kind_phys)*tplat)
            else
               tprime0 = tplat
            end if

            if (tprime0 < tstar) then
               alnp = 0._kind_phys
            else
               alnp = rair*(tprime0 - tstar)/phis(i) * log(pout/ps(i))
            end if

         end if

         extrapolate_T = tstar*(1._kind_phys + alnp + 0.5_kind_phys*alnp**2 + sixth*alnp**3)

      end function extrapolate_T

      !======================================================================================

      real(kind_phys) function extrapolate_Z()

         ! local variables
         real(kind_phys) :: sixth
         real(kind_phys) :: hgt
         real(kind_phys) :: tstar
         real(kind_phys) :: t0
         real(kind_phys) :: alph
         real(kind_phys) :: alnp
         !-------------------------------------------------------------------------

         sixth  = 1._kind_phys / 6._kind_phys
         hgt   = phis(i)*rga
         tstar = tbot(i)*(1._kind_phys + alpha*(ps(i)/pin(i,nlev) - 1._kind_phys))
         t0 = tstar + 0.0065_kind_phys*hgt

         if (tstar <= 290.5_kind_phys .and. t0 > 290.5_kind_phys) then
            alph = rair/phis(i) * (290.5_kind_phys - tstar)

         else if (tstar > 290.5_kind_phys .and. t0 > 290.5_kind_phys) then
            alph  = 0._kind_phys
            tstar = 0.5_kind_phys*(290.5_kind_phys + tstar)

         else
            alph = alpha

         end if

         if (tstar < 255._kind_phys) then
            tstar = 0.5_kind_phys*(tstar + 255._kind_phys)
         end if
         alnp = alph*log(pout/ps(i))

         extrapolate_Z = hgt - rair*tstar*rga*log(pout/ps(i))* &
                         (1._kind_phys + 0.5_kind_phys*alnp + sixth*alnp**2)

      end function extrapolate_Z

   end subroutine vertinterp

end module atmos_phys_interp_utils
