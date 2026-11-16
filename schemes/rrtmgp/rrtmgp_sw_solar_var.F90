!-------------------------------------------------------------------------------
! This module uses the solar irradiance data 
! to provide a spectral scaling factor
! to approximate the spectral distribution of irradiance
! when the radiation scheme might use a different solar source function
!-------------------------------------------------------------------------------
module rrtmgp_sw_solar_var

  implicit none

  private
  public :: rrtmgp_sw_solar_var_run

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> \section arg_table_rrtmgp_sw_solar_var_run Argument Table
!! \htmlinclude rrtmgp_sw_solar_var_run.html
!!
  subroutine rrtmgp_sw_solar_var_run(toa_flux, constant_dim_two, band2gpt_sw, nswbands, sol_irrad, wave_end, nbins, sol_tsi, &
                                     nday, dosw, do_spectral_scaling, sfac, eccf, errmsg, errflg)
     use rrtmgp_sw_solar_var_setup, only: irrad, radbinmax, radbinmin
     use ccpp_kinds,                only : kind_phys

     ! Arguments 
     real(kind_phys),    intent(inout) :: toa_flux(:,:)         ! top-of-atmosphere flux to be scaled (columns,gpts)
     real(kind_phys),    intent(in)    :: sol_tsi               ! total solar irradiance
     real(kind_phys),    intent(in)    :: sol_irrad(:)          ! solar irradiance
     real(kind_phys),    intent(in)    :: wave_end(:)           ! wavelength endpoints
     integer,            intent(in)    :: nday                  ! number of daytime points
     integer,            intent(in)    :: nbins                 ! number of bins
     integer,            intent(in)    :: constant_dim_two      ! dimension for band2gpt_sw
     integer,            intent(in)    :: band2gpt_sw(:,:)      ! array for converting shortwave band limits to g-points
     integer,            intent(in)    :: nswbands              ! number of shortwave bands
     logical,            intent(in)    :: do_spectral_scaling   ! flag to do spectral scaling
     logical,            intent(in)    :: dosw                  ! flag to do shortwave radiation
     real(kind_phys),    intent(in)    :: eccf                  ! Earth-Sun distance factor
     real(kind_phys),    intent(out)   :: sfac(:,:)             ! scaling factors (columns,gpts)
     character(len=*),   intent(out)   :: errmsg
     integer,            intent(out)   :: errflg

     ! Local variables 
     integer :: i, j, ib
     real(kind_phys) :: s                                  ! g-point sum of toa_flux
     real(kind_phys) :: band_irrad(nswbands)               ! irrad, as data the kernel reads
     integer         :: band_gpts(2, nswbands)             ! band2gpt_sw, likewise

     ! Initialize error variables
     errflg = 0
     errmsg = ''

     if (.not. dosw .or. nday == 0) then
        return
     end if

     if (do_spectral_scaling) then
        ! Determine target irradiance for each band (on the host)
        call integrate_spectrum(nbins, nswbands, wave_end, radbinmin, radbinmax, sol_irrad, irrad)
        band_irrad(:) = irrad(:nswbands)
     else
        band_irrad(:) = 0.0_kind_phys
     end if
     band_gpts(:, :) = band2gpt_sw(:, :nswbands)

     ! One column per thread; the g-point sums add left to right from zero,
     ! as SUM does on the host: the inner loops are sequential (loop seq), or
     ! the compiler may vectorize them into reductions of another order.
     !$acc parallel loop private(s) copyin(band_irrad, band_gpts) copy(toa_flux) copyout(sfac)
     do i = 1, size(toa_flux, 1)
        if (do_spectral_scaling) then
           !$acc loop seq
           do ib = 1, nswbands
              s = 0.0_kind_phys
              !$acc loop seq
              do j = band_gpts(1, ib), band_gpts(2, ib)
                 s = s + toa_flux(i, j)
              end do
              !$acc loop seq
              do j = band_gpts(1, ib), band_gpts(2, ib)
                 sfac(i, j) = band_irrad(ib) / s
              end do
           end do
        else
           s = 0.0_kind_phys
           !$acc loop seq
           do j = 1, size(toa_flux, 2)
              s = s + toa_flux(i, j)
           end do
           !$acc loop seq
           do j = 1, size(toa_flux, 2)
              sfac(i, j) = sol_tsi / s
           end do
        end if
        !$acc loop seq
        do j = 1, size(toa_flux, 2)
           toa_flux(i, j) = toa_flux(i, j) * sfac(i, j) * eccf
        end do
     end do

  end subroutine rrtmgp_sw_solar_var_run


!-------------------------------------------------------------------------------
! private method.........
!-------------------------------------------------------------------------------

  subroutine integrate_spectrum( nsrc, ntrg, src_x, min_trg, max_trg, src, trg )

    use ccpp_tuvx_utils, only : rebin
    use ccpp_kinds,      only : kind_phys

    implicit none

    !---------------------------------------------------------------
    !	... dummy arguments
    !---------------------------------------------------------------
    integer,  intent(in)  :: nsrc                  ! dimension source array
    integer,  intent(in)  :: ntrg                  ! dimension target array
    real(kind_phys), intent(in)  :: src_x(nsrc+1)         ! source coordinates
    real(kind_phys), intent(in)  :: max_trg(ntrg)         ! target coordinates
    real(kind_phys), intent(in)  :: min_trg(ntrg)         ! target coordinates
    real(kind_phys), intent(in)  :: src(nsrc)             ! source array
    real(kind_phys), intent(out) :: trg(ntrg)             ! target array
 
    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    real(kind_phys) :: trg_x(2), targ(1)         ! target coordinates
    integer  :: i

    do i = 1, ntrg

       trg_x(1) = min_trg(i)
       trg_x(2) = max_trg(i)

       call rebin( nsrc, 1, src_x, trg_x, src(1:nsrc), targ(:) )
       ! W m-2 nm-1 --> W m-2
       trg( i ) = targ(1)*(trg_x(2)-trg_x(1))

    end do

  end subroutine integrate_spectrum

end module rrtmgp_sw_solar_var
