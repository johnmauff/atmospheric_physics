! This module contains code computing the effective diffusion of
! constituents and dry static energy due to gravity wave breaking.
!
! The routines work on one column (column index icol, arrays passed whole)
! and are OpenACC `routine seq`, so gw_drag_prof_column can call them from a
! column kernel or from host code. The implicit solve uses the plain-array
! column solver in atmos_phys_fin_vol_solver (the arithmetic of
! vdiff_lu_solver's fin_vol_lu_decomp and TriDiagDecomp's left_div, operation
! by operation) instead of the derived types, which cannot be used in GPU
! code. The solver arrays index levels from 1, so ktop must be 1 (gw_common's
! ktop parameter).
module gw_diffusion

  use ccpp_kinds, only: kind_phys

  implicit none
  private

  public :: gw_ediff_column
  public :: gw_diff_tend_column
  public :: gw_diff_tend_column_3d

contains

  ! Calculate effective diffusivity associated with GW forcing, and the LU
  ! decomposition (ca, ze, dnom) for the diffusion solve, for column icol.
  ! Author: F. Sassi, Jan 31, 2001
  subroutine gw_ediff_column(icol, ngwv, kbot, ktop, tend_level, &
                             gwut, ubm, nm, rho, dt, prndl, gravit, &
                             del, rdel, rdst, c, vramp, &
                             egwdffi, egwdffm, dpidz_sq, ca, ze, dnom, ro_adjust)
    !$acc routine seq

    use atmos_phys_fin_vol_solver, only: fin_vol_decomp_column

    ! Column index and gravity wave spectrum dimension.
    integer, intent(in) :: icol, ngwv
    ! Bottom and top levels to operate on.
    integer, intent(in) :: kbot, ktop
    ! Per-column bottom index where tendencies are applied.
    integer, intent(in) :: tend_level(:)
    ! GW zonal wind tendencies at midpoint.
    real(kind_phys), intent(in) :: gwut(:, :, -ngwv:)
    ! Projection of wind at midpoints.
    real(kind_phys), intent(in) :: ubm(:, :)
    ! Brunt-Vaisalla frequency.
    real(kind_phys), intent(in) :: nm(:, :)
    ! Density at interfaces.
    real(kind_phys), intent(in) :: rho(:, :)
    ! Time step.
    real(kind_phys), intent(in) :: dt
    ! Inverse Prandtl number.
    real(kind_phys), intent(in) :: prndl
    ! Acceleration due to gravity.
    real(kind_phys), intent(in) :: gravit
    ! Pressure coordinates: layer thickness, its inverse, and the inverse
    ! distance between midpoints (Coords1D's del, rdel, rdst).
    real(kind_phys), intent(in) :: del(:, :), rdel(:, :), rdst(:, :)
    ! Wave phase speeds for each column.
    real(kind_phys), intent(in) :: c(:, -ngwv:)
    ! Coefficient to ramp down diffusion coeff.
    real(kind_phys), intent(in) :: vramp(:)
    ! Effective gw diffusivity at interfaces.
    real(kind_phys), intent(inout) :: egwdffi(:, :)
    ! Work arrays (ncol, pver), (ncol, pverp): effective gw diffusivity at
    ! midpoints, and (dp/dz)^2 == (gravit*rho)^2.
    real(kind_phys), intent(inout) :: egwdffm(:, :), dpidz_sq(:, :)
    ! LU decomposition (ncol, pver).
    real(kind_phys), intent(inout) :: ca(:, :), ze(:, :), dnom(:, :)
    ! Adjustment parameter for IGWs.
    real(kind_phys), intent(in), optional :: ro_adjust(:, -ngwv:, :)

    ! Temporary used to hold gw_diffusivity for one level and wavenumber.
    real(kind_phys) :: egwdff_lev
    ! Effective gw diffusivity at one midpoint, summed over waves.
    real(kind_phys) :: egwdff_sum
    ! dp/dz == gravit*rho at one interface.
    real(kind_phys) :: dpidz
    ! Level and wave indices.
    integer :: k, l

    do k = 1, size(egwdffi, 2)
      egwdffi(icol, k) = 0._kind_phys
    end do
    do k = 1, size(egwdffm, 2)
      egwdffm(icol, k) = 0._kind_phys
    end do

    ! Calculate effective diffusivity at midpoints (summed over waves in
    ! wave order, as before).
    do k = ktop, kbot
      egwdff_sum = egwdffm(icol, k)
      do l = -ngwv, ngwv
        egwdff_lev = &
          prndl*0.5_kind_phys*gwut(icol, k, l)*(c(icol, l) - ubm(icol, k))/nm(icol, k)**2
        ! IGWs need ro_adjust factor.
        if (present(ro_adjust)) then
          egwdff_lev = egwdff_lev*ro_adjust(icol, l, k)**2
        end if
        egwdff_sum = egwdff_sum + egwdff_lev
      end do
      egwdffm(icol, k) = egwdff_sum*vramp(k)
    end do

    ! Interpolate effective diffusivity to interfaces.
    ! Assume zero at top and bottom interfaces.
    do k = ktop + 1, kbot
      egwdffi(icol, k) = 0.5_kind_phys*(egwdffm(icol, k - 1) + egwdffm(icol, k))
    end do

    ! Do not calculate diffusivities below level where tendencies are
    ! actually allowed.
    do k = ktop + 1, kbot
      if (k > tend_level(icol)) egwdffi(icol, k) = 0.0_kind_phys
    end do

    ! Calculate (dp/dz)^2 (in a scalar, stored once).
    do k = 1, size(dpidz_sq, 2)
      dpidz = rho(icol, k)*gravit
      dpidz_sq(icol, k) = dpidz*dpidz
    end do

    ! Decompose the diffusion matrix: levels ktop..kbot, coefficient
    ! egwdffi*dpidz_sq at their interfaces, default boundaries.
    call fin_vol_decomp_column(icol, kbot - ktop + 1, dt, del, rdel, rdst, &
                               egwdffi, dpidz_sq, ca, ze, dnom)

  end subroutine gw_ediff_column

  !==========================================================================

  ! Calculates tendencies from effective diffusion due to gravity wave
  ! breaking.
  !
  ! Method:
  ! A constituent flux on interfaces is given by:
  !
  !              rho * (w'q') = rho * Deff qz
  !
  ! where (all evaluated on interfaces):
  !
  !        rho   = density
  !        qz    = constituent vertical gradient
  !        Deff  = effective diffusivity
  !
  ! An effective diffusivity is calculated by adding up the diffusivities
  ! from all waves (see gw_ediff_column). The tendency is calculated by
  ! invoking LU decomposition and solving as for a regular diffusion equation.
  !
  ! Author: Sassi - Jan 2001
  subroutine gw_diff_tend_column(icol, kbot, ktop, q, dt, ca, ze, dnom, dq, qnew)
    !$acc routine seq

    use atmos_phys_fin_vol_solver, only: fin_vol_left_div_column

    ! Column index.
    integer, intent(in) :: icol
    ! Bottom and top levels to operate on.
    integer, intent(in) :: kbot, ktop
    ! Constituent to diffuse (ncol, pver).
    real(kind_phys), intent(in) :: q(:, :)
    ! Time step.
    real(kind_phys), intent(in) :: dt
    ! LU decomposition.
    real(kind_phys), intent(in) :: ca(:, :), ze(:, :), dnom(:, :)
    ! Constituent tendencies.
    real(kind_phys), intent(inout) :: dq(:, :)
    ! Work array (ncol, pver): temporary storage for constituent.
    real(kind_phys), intent(inout) :: qnew(:, :)

    integer :: k

    do k = 1, size(q, 2)
      qnew(icol, k) = q(icol, k)
    end do
    call fin_vol_left_div_column(icol, kbot - ktop + 1, ca, ze, dnom, qnew)

    ! Evaluate tendency to be reported back.
    do k = 1, size(q, 2)
      dq(icol, k) = (qnew(icol, k) - q(icol, k))/dt
    end do

  end subroutine gw_diff_tend_column

  ! gw_diff_tend_column for constituent m of q(:,:,:), tendency into
  ! dq(:,:,m), so device callers do not pass array sections.
  subroutine gw_diff_tend_column_3d(icol, kbot, ktop, m, q, dt, ca, ze, dnom, dq, qnew)
    !$acc routine seq

    use atmos_phys_fin_vol_solver, only: fin_vol_left_div_column

    integer, intent(in) :: icol
    integer, intent(in) :: kbot, ktop
    ! Constituent index.
    integer, intent(in) :: m
    real(kind_phys), intent(in) :: q(:, :, :)
    real(kind_phys), intent(in) :: dt
    real(kind_phys), intent(in) :: ca(:, :), ze(:, :), dnom(:, :)
    real(kind_phys), intent(inout) :: dq(:, :, :)
    real(kind_phys), intent(inout) :: qnew(:, :)

    integer :: k

    do k = 1, size(q, 2)
      qnew(icol, k) = q(icol, k, m)
    end do
    call fin_vol_left_div_column(icol, kbot - ktop + 1, ca, ze, dnom, qnew)

    do k = 1, size(q, 2)
      dq(icol, k, m) = (qnew(icol, k) - q(icol, k, m))/dt
    end do

  end subroutine gw_diff_tend_column_3d

end module gw_diffusion
