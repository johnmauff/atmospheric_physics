! Per-column implicit finite-volume vertical diffusion on plain arrays.
!
! The arithmetic of coords_1d's Coords1D(ifc), linear_1d_operators'
! diffusion_operator, TriDiagDecomp and left_div, and of
! vertical_diffusion_solver's fin_vol_solve / vdiff_lu_solver's
! fin_vol_lu_decomp, done operation by operation for one column, so that
! callers get bit-for-bit the same results without the derived types (which
! cannot be used in GPU code). Implemented: a diffusion coefficient
! kv*dpidz_sq at interfaces, an optional diagonal term (coef_q), the default
! (fixed-flux) boundaries or a fixed layer above the top interface, and an
! optional fixed value at the top (BoundaryData). Not implemented: advection
! terms, coef_q_weight, other boundary types (use linear_1d_operators).
!
! Arrays are passed whole, (ncol, ...), with the column index icol; the
! routines are OpenACC `routine seq`, called once per column from a column
! kernel (or from host code on the CPU).
module atmos_phys_fin_vol_solver

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: fin_vol_coords_column   ! Coords1D(pint) for one column
  public :: fin_vol_decomp_column   ! LU decomposition of (I - dt*L) for one column
  public :: fin_vol_left_div_column ! Solve with a decomposition for one column
  public :: fin_vol_left_div_column_3d ! The same for q(:,:,m) of a 3-D array

contains

  ! Pressure coordinates of column icol from its interface pressures, as
  ! Coords1D(ifc): mid = 0.5*(ifc(k)+ifc(k+1)), del = ifc(k+1)-ifc(k),
  ! rdel = 1/del, rdst = 1/(mid(k+1)-mid(k)).
  subroutine fin_vol_coords_column(icol, nlev, pint, mid, del, rdel, rdst)
    !$acc routine seq
    integer,         intent(in)    :: icol
    integer,         intent(in)    :: nlev        ! Number of layers
    real(kind_phys), intent(in)    :: pint(:,:)   ! Interface pressures (ncol, nlev+1) [Pa]
    real(kind_phys), intent(inout) :: mid(:,:)    ! Midpoints (ncol, nlev) [Pa]
    real(kind_phys), intent(inout) :: del(:,:)    ! Layer thickness (ncol, nlev) [Pa]
    real(kind_phys), intent(inout) :: rdel(:,:)   ! 1/del (ncol, nlev) [Pa-1]
    real(kind_phys), intent(inout) :: rdst(:,:)   ! 1/(midpoint distance) (ncol, nlev-1) [Pa-1]

    integer :: k

    do k = 1, nlev
      mid(icol,k) = 0.5_kind_phys * (pint(icol,k) + pint(icol,k+1))
      del(icol,k) = pint(icol,k+1) - pint(icol,k)
      rdel(icol,k) = 1._kind_phys/del(icol,k)
    end do
    do k = 1, nlev-1
      rdst(icol,k) = 1._kind_phys/(mid(icol,k+1) - mid(icol,k))
    end do

  end subroutine fin_vol_coords_column

  ! LU decomposition of column icol for the implicit step
  ! (I - dt*(L + diag(coef_q))) q_new = q, where L is the finite-volume
  ! diffusion operator with coefficient kv*dpidz_sq at interfaces
  ! (fin_vol_lu_decomp / fin_vol_solve with coef_q_diff=kv*dpidz_sq and, if
  ! present, coef_q). Top boundary: fixed flux (default) or, if top_width is
  ! present, a fixed layer of that width (BoundaryFixedLayer); bottom: fixed
  ! flux. ca, ze and dnom are TriDiagDecomp's arrays.
  subroutine fin_vol_decomp_column(icol, nlev, dt, del, rdel, rdst, kv, dpidz_sq, &
                                   ca, ze, dnom, coef_q, top_width)
    !$acc routine seq
    integer,         intent(in)    :: icol
    integer,         intent(in)    :: nlev          ! Number of layers
    real(kind_phys), intent(in)    :: dt            ! Time step [s]
    real(kind_phys), intent(in)    :: del(:,:)      ! Layer thickness (ncol, nlev) [Pa]
    real(kind_phys), intent(in)    :: rdel(:,:)     ! 1/del (ncol, nlev) [Pa-1]
    real(kind_phys), intent(in)    :: rdst(:,:)     ! 1/(midpoint distance) (ncol, nlev-1) [Pa-1]
    real(kind_phys), intent(in)    :: kv(:,:)       ! Diffusivity at interfaces (ncol, nlev+1)
    real(kind_phys), intent(in)    :: dpidz_sq(:,:) ! (dp/dz)**2 at interfaces (ncol, nlev+1)
    real(kind_phys), intent(inout) :: ca(:,:)       ! (ncol, nlev)
    real(kind_phys), intent(inout) :: ze(:,:)       ! (ncol, nlev)
    real(kind_phys), intent(inout) :: dnom(:,:)     ! (ncol, nlev)
    real(kind_phys), intent(in), optional :: coef_q(:,:) ! Diagonal term (ncol, nlev) [s-1]
    real(kind_phys), intent(in), optional :: top_width   ! Fixed-layer width above the top interface [Pa]

    integer         :: k
    real(kind_phys) :: c                 ! Scaling of the operator, -dt
    real(kind_phys) :: d_top             ! Diffusion coefficient at the top interface
    real(kind_phys) :: d_kp1, d_k        ! Diffusion coefficient at interfaces k+1, k
    real(kind_phys) :: left_bound        ! Operator's left (top) boundary term, unscaled
    real(kind_phys) :: right_bound       ! Operator's right (bottom) boundary term, unscaled
    real(kind_phys) :: spr, sub          ! Super- and subdiagonal of row k, unscaled
    real(kind_phys) :: sub_km1           ! Subdiagonal of row k-1, unscaled
    real(kind_phys) :: diag              ! Diagonal of row k
    real(kind_phys) :: spr_s, sub_s      ! Scaled super- and subdiagonal of row k

    c = -dt

    ! diffusion_operator: boundary terms
    if (present(top_width)) then
      d_top = kv(icol,1)*dpidz_sq(icol,1)
      left_bound = 2._kind_phys*d_top*rdel(icol,1) / (top_width + del(icol,1))
    else
      left_bound = 0._kind_phys
    end if
    right_bound = 0._kind_phys

    ! Bottom row (k = nlev): diag = -right_bound - sub(nlev-1) (deriv_diag),
    ! + coef_q, scaled, + 1 (add_to_diag); new_TriDiagDecomp's last row.
    d_k = kv(icol,nlev)*dpidz_sq(icol,nlev)
    sub_km1 = d_k*rdst(icol,nlev-1)*rdel(icol,nlev)
    diag = - right_bound
    diag = diag - sub_km1
    if (present(coef_q)) diag = diag + coef_q(icol,nlev)
    diag = diag*c
    diag = diag + 1._kind_phys
    ca(icol,nlev) = - (right_bound*c)
    dnom(icol,nlev) = 1._kind_phys / diag

    ! Rows nlev-1 .. 1 (new_TriDiagDecomp's backward loop)
    do k = nlev - 1, 1, -1
      d_kp1 = kv(icol,k+1)*dpidz_sq(icol,k+1)
      spr = d_kp1*rdst(icol,k)*rdel(icol,k)
      sub = d_kp1*rdst(icol,k)*rdel(icol,k+1)
      diag = - spr
      if (k == 1) then
        diag = diag - left_bound
      else
        d_k = kv(icol,k)*dpidz_sq(icol,k)
        sub_km1 = d_k*rdst(icol,k-1)*rdel(icol,k)
        diag = diag - sub_km1
      end if
      if (present(coef_q)) diag = diag + coef_q(icol,k)
      spr_s = spr*c
      sub_s = sub*c
      diag = diag*c
      diag = diag + 1._kind_phys

      ca(icol,k) = - spr_s
      ze(icol,k+1) = - sub_s * dnom(icol,k+1)
      dnom(icol,k) = 1._kind_phys / (diag - ca(icol,k)*ze(icol,k+1))
    end do
    ze(icol,1) = - (left_bound*c)

  end subroutine fin_vol_decomp_column

  ! Solve column icol of q in place with a decomposition from
  ! fin_vol_decomp_column (TriDiagDecomp's left_div). If top_value is present,
  ! the top boundary holds that value (BoundaryData on the left).
  subroutine fin_vol_left_div_column(icol, nlev, ca, ze, dnom, q, top_value)
    !$acc routine seq
    integer,         intent(in)    :: icol
    integer,         intent(in)    :: nlev        ! Number of layers
    real(kind_phys), intent(in)    :: ca(:,:)     ! (ncol, nlev)
    real(kind_phys), intent(in)    :: ze(:,:)     ! (ncol, nlev)
    real(kind_phys), intent(in)    :: dnom(:,:)   ! (ncol, nlev)
    real(kind_phys), intent(inout) :: q(:,:)      ! (ncol, nlev)
    real(kind_phys), intent(in), optional :: top_value

    integer :: k

    if (present(top_value)) then
      q(icol,1) = q(icol,1) + ze(icol,1)*top_value
    end if

    ! Backward sweep, keeping zf in q
    q(icol,nlev) = q(icol,nlev) * dnom(icol,nlev)
    do k = nlev - 1, 1, -1
      q(icol,k) = (q(icol,k) + ca(icol,k)*q(icol,k+1)) * dnom(icol,k)
    end do
    ! Forward sweep (q(1) = zf(1) already)
    do k = 2, nlev
      q(icol,k) = q(icol,k) + ze(icol,k)*q(icol,k-1)
    end do

  end subroutine fin_vol_left_div_column

  ! fin_vol_left_div_column for q(:,:,m) of a 3-D array (no top value), so
  ! device callers do not pass an array section.
  subroutine fin_vol_left_div_column_3d(icol, nlev, m, ca, ze, dnom, q)
    !$acc routine seq
    integer,         intent(in)    :: icol
    integer,         intent(in)    :: nlev        ! Number of layers
    integer,         intent(in)    :: m           ! Index in the third dimension of q
    real(kind_phys), intent(in)    :: ca(:,:)     ! (ncol, nlev)
    real(kind_phys), intent(in)    :: ze(:,:)     ! (ncol, nlev)
    real(kind_phys), intent(in)    :: dnom(:,:)   ! (ncol, nlev)
    real(kind_phys), intent(inout) :: q(:,:,:)    ! (ncol, nlev, :)

    integer :: k

    q(icol,nlev,m) = q(icol,nlev,m) * dnom(icol,nlev)
    do k = nlev - 1, 1, -1
      q(icol,k,m) = (q(icol,k,m) + ca(icol,k)*q(icol,k+1,m)) * dnom(icol,k)
    end do
    do k = 2, nlev
      q(icol,k,m) = q(icol,k,m) + ze(icol,k)*q(icol,k-1,m)
    end do

  end subroutine fin_vol_left_div_column_3d

end module atmos_phys_fin_vol_solver
