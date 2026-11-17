module atmos_phys_interp_point
   ! Linear interpolation to a single point, callable from OpenACC kernels.
   !
   ! These routines give the same results as atmos_phys_interp_utils'
   ! lininterp_init (extrap_method_bndry) followed by lininterp for one
   ! output point (nout = m1 = 1), without its allocated interp_type: the
   ! weights come back as scalars. The search, boundary handling and the
   ! order of the arithmetic are those of lininterp_init, lininterp1d and
   ! lininterp2d1d, so results are bit-for-bit.
   !
   ! Errors are returned as an integer code (0 = none); use
   ! interp_point_errmsg on the host to turn a code into a message.

   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: interp_point_wgts        ! weights on a 1-d mesh
   public :: interp_point_wgts_blend  ! weights on a mesh blended from two table rows
   public :: interp_point_1d          ! interpolate a 1-d array
   public :: interp_point_2d          ! interpolate a 2-d array (bilinear)
   public :: interp_point_errmsg

   integer, parameter, public :: interp_point_err_none        = 0
   integer, parameter, public :: interp_point_err_short_mesh  = 1
   integer, parameter, public :: interp_point_err_nonmonotone = 2
   integer, parameter, public :: interp_point_err_bad_weight  = 3
   integer, parameter, public :: interp_point_err_no_bracket  = 4

contains

   subroutine interp_point_wgts(yin, nin, yout, jjm, jjp, wgts, wgtn, err)
      ! Weights to interpolate from mesh yin (monotonic, increasing or
      ! decreasing) to the point yout; beyond the mesh, the boundary value.
      !$acc routine seq
      integer,         intent(in)  :: nin
      real(kind_phys), intent(in)  :: yin(:)
      real(kind_phys), intent(in)  :: yout
      integer,         intent(out) :: jjm, jjp
      real(kind_phys), intent(out) :: wgts, wgtn
      integer,         intent(out) :: err

      integer :: j, jj, icount
      logical :: increasing
      real(kind_phys) :: ratio

      jjm = 0
      jjp = 0
      wgts = 0._kind_phys
      wgtn = 0._kind_phys
      err = interp_point_err_none
      if (nin < 2) then
         err = interp_point_err_short_mesh
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
         err = interp_point_err_nonmonotone
         return
      end if

      ! Beyond the mesh: copy the boundary value
      if (increasing) then
         if (yout <= yin(1)) then
            jjm = 1
            jjp = 1
            wgts = 1._kind_phys
            wgtn = 0._kind_phys
         else if (yout > yin(nin)) then
            jjm = nin
            jjp = nin
            wgts = 1._kind_phys
            wgtn = 0._kind_phys
         end if
      else
         if (yout > yin(1)) then
            jjm = 1
            jjp = 1
            wgts = 1._kind_phys
            wgtn = 0._kind_phys
         else if (yout <= yin(nin)) then
            jjm = nin
            jjp = nin
            wgts = 1._kind_phys
            wgtn = 0._kind_phys
         end if
      end if

      ! Inside the mesh
      if (increasing) then
         do jj = 1, nin-1
            if (yout > yin(jj) .and. yout <= yin(jj+1)) then
               jjm = jj
               jjp = jj + 1
               wgts = (yin(jj+1)-yout)/(yin(jj+1)-yin(jj))
               wgtn = (yout-yin(jj))/(yin(jj+1)-yin(jj))
               exit
            end if
         end do
      else
         do jj = 1, nin-1
            if (yout <= yin(jj) .and. yout > yin(jj+1)) then
               jjm = jj
               jjp = jj + 1
               wgts = (yin(jj+1)-yout)/(yin(jj+1)-yin(jj))
               wgtn = (yout-yin(jj))/(yin(jj+1)-yin(jj))
               exit
            end if
         end do
      end if

      if (jjm == 0 .or. jjp == 0) then
         ! yout is NaN
         err = interp_point_err_no_bracket
         return
      end if
      ratio = wgts + wgtn
      if (ratio < 0.9_kind_phys .or. ratio > 1.1_kind_phys) then
         err = interp_point_err_bad_weight
      end if

   end subroutine interp_point_wgts

   subroutine interp_point_wgts_blend(tab, nin, im, ip, bs, bn, yout, &
                                      jjm, jjp, wgts, wgtn, err)
      ! As interp_point_wgts, on the mesh yin(j) = tab(im,j)*bs + tab(ip,j)*bn, j = 1..nin
      ! (rows im and ip of tab blended with the weights of an earlier
      ! interpolation). Each yin(j) is computed as lininterp1d computes it,
      ! so the result equals interpolating that row with lininterp first.
      !$acc routine seq
      integer,         intent(in)  :: nin
      real(kind_phys), intent(in)  :: tab(:,:)
      integer,         intent(in)  :: im, ip
      real(kind_phys), intent(in)  :: bs, bn
      real(kind_phys), intent(in)  :: yout
      integer,         intent(out) :: jjm, jjp
      real(kind_phys), intent(out) :: wgts, wgtn
      integer,         intent(out) :: err

      integer :: j, jj, icount
      logical :: increasing
      real(kind_phys) :: ratio
      real(kind_phys) :: y1, y2, ylo, yhi

      jjm = 0
      jjp = 0
      wgts = 0._kind_phys
      wgtn = 0._kind_phys
      err = interp_point_err_none
      if (nin < 2) then
         err = interp_point_err_short_mesh
         return
      end if

      increasing = .true.
      icount = 0
      do j = 1, nin-1
         y1 = tab(im,j)*bs + tab(ip,j)*bn
         y2 = tab(im,j+1)*bs + tab(ip,j+1)*bn
         if (y1 > y2) icount = icount + 1
      end do
      if (icount == nin-1) then
         increasing = .false.
         icount = 0
      end if
      if (icount > 0) then
         err = interp_point_err_nonmonotone
         return
      end if

      ! Beyond the mesh: copy the boundary value
      y1 = tab(im,1)*bs + tab(ip,1)*bn
      y2 = tab(im,nin)*bs + tab(ip,nin)*bn
      if (increasing) then
         if (yout <= y1) then
            jjm = 1
            jjp = 1
            wgts = 1._kind_phys
            wgtn = 0._kind_phys
         else if (yout > y2) then
            jjm = nin
            jjp = nin
            wgts = 1._kind_phys
            wgtn = 0._kind_phys
         end if
      else
         if (yout > y1) then
            jjm = 1
            jjp = 1
            wgts = 1._kind_phys
            wgtn = 0._kind_phys
         else if (yout <= y2) then
            jjm = nin
            jjp = nin
            wgts = 1._kind_phys
            wgtn = 0._kind_phys
         end if
      end if

      ! Inside the mesh
      do jj = 1, nin-1
         ylo = tab(im,jj)*bs + tab(ip,jj)*bn
         yhi = tab(im,jj+1)*bs + tab(ip,jj+1)*bn
         if (increasing) then
            if (yout > ylo .and. yout <= yhi) then
               jjm = jj
               jjp = jj + 1
               wgts = (yhi-yout)/(yhi-ylo)
               wgtn = (yout-ylo)/(yhi-ylo)
               exit
            end if
         else
            if (yout <= ylo .and. yout > yhi) then
               jjm = jj
               jjp = jj + 1
               wgts = (yhi-yout)/(yhi-ylo)
               wgtn = (yout-ylo)/(yhi-ylo)
               exit
            end if
         end if
      end do

      if (jjm == 0 .or. jjp == 0) then
         ! yout is NaN
         err = interp_point_err_no_bracket
         return
      end if
      ratio = wgts + wgtn
      if (ratio < 0.9_kind_phys .or. ratio > 1.1_kind_phys) then
         err = interp_point_err_bad_weight
      end if

   end subroutine interp_point_wgts_blend

   function interp_point_1d(arrin, jjm, jjp, wgts, wgtn) result(arrout)
      ! lininterp1d for one output point
      !$acc routine seq
      real(kind_phys), intent(in) :: arrin(:)
      integer,         intent(in) :: jjm, jjp
      real(kind_phys), intent(in) :: wgts, wgtn
      real(kind_phys) :: arrout

      arrout = arrin(jjm)*wgts + arrin(jjp)*wgtn

   end function interp_point_1d

   function interp_point_2d(arrin, iim, iip, wgtw, wgte, jjm, jjp, wgts, wgtn) &
        result(arrout)
      ! lininterp2d1d for one output point: (iim, iip, wgtw, wgte) are the
      ! weights in the first dimension, (jjm, jjp, wgts, wgtn) in the second
      !$acc routine seq
      real(kind_phys), intent(in) :: arrin(:,:)
      integer,         intent(in) :: iim, iip, jjm, jjp
      real(kind_phys), intent(in) :: wgtw, wgte, wgts, wgtn
      real(kind_phys) :: arrout

      arrout = arrin(iim,jjm)*wgtw*wgts+arrin(iip,jjm)*wgte*wgts + &
               arrin(iim,jjp)*wgtw*wgtn+arrin(iip,jjp)*wgte*wgtn

   end function interp_point_2d

   function interp_point_errmsg(err) result(msg)
      ! Host-side message for an error code from interp_point_wgts*
      integer, intent(in) :: err
      character(len=128)  :: msg

      select case (err)
      case (interp_point_err_none)
         msg = ''
      case (interp_point_err_short_mesh)
         msg = 'LININTERP: Must have at least 2 input points for interpolation'
      case (interp_point_err_nonmonotone)
         msg = 'LININTERP: Non-monotonic input coordinate array found'
      case (interp_point_err_bad_weight)
         msg = 'LININTERP: Bad weight computed in lininterp_init'
      case (interp_point_err_no_bracket)
         msg = 'LININTERP: Point found without interp indices'
      case default
         msg = 'LININTERP: unknown error'
      end select

   end function interp_point_errmsg

end module atmos_phys_interp_point
