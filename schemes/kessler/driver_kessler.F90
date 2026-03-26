program test_kessler_driver

  use kessler
  implicit none

  integer, parameter :: kind_phys = selected_real_kind(12)

  ! Dimensions
  integer :: ncol, nz
  integer :: lyr_surf, lyr_toa
  real(kind_phys) :: dt

  ! Init constants
  real(kind_phys) :: lv_in, pref_in, rhoqr_in

  ! Arrays
  real(kind_phys), allocatable :: cpair(:,:), rair(:,:), rho(:,:)
  real(kind_phys), allocatable :: z(:,:), pk(:,:)
  real(kind_phys), allocatable :: theta(:,:), qv(:,:), qc(:,:), qr(:,:)
  real(kind_phys), allocatable :: precl(:)
  real(kind_phys), allocatable :: relhum(:,:)
  real(kind_phys), allocatable :: arr(:)
  real(kind_phys) :: u1, u2, ztmp

  character(len=64)  :: scheme_name
  character(len=512) :: errmsg
  integer            :: errflg
  integer, allocatable :: seed_values(:)
  integer :: seed_size


  integer :: i, k

  !------------------------------------------------------
  ! Set grid size
  !------------------------------------------------------
  ncol = 128
  nz   = 56
  dt   = 60.0_kind_phys

  lyr_surf = 1
  lyr_toa  = nz

  !------------------------------------------------------
  ! Initialize Kessler constants
  !------------------------------------------------------
  lv_in    = 2.5e6_kind_phys     ! J/kg
  pref_in  = 100000.0_kind_phys  ! Pa (1000 hPa)
  rhoqr_in = 1000.0_kind_phys    ! kg/m^3

  call kessler_init(lv_in, pref_in, rhoqr_in, errmsg, errflg)

  if (errflg /= 0) then
     print *, 'Initialization error: ', trim(errmsg)
     stop
  end if

  !------------------------------------------------------
  ! Allocate arrays
  !------------------------------------------------------
  allocate(cpair(ncol,nz), rair(ncol,nz), rho(ncol,nz))
  allocate(z(ncol,nz), pk(ncol,nz))
  allocate(theta(ncol,nz), qv(ncol,nz), qc(ncol,nz), qr(ncol,nz))
  allocate(precl(ncol))
  allocate(relhum(ncol,nz))
  allocate(arr(ncol))

  ! Query the size of the RNG state
  call random_seed(size=seed_size)
  allocate(seed_values(seed_size))

  ! Fill with fixed values for deterministic sequence
  seed_values = [(i, i=1, seed_size)]  ! or any fixed sequence

  ! Set the seed
  call random_seed(put=seed_values)

  do i = 1, ncol
     ! Box-Muller transform to generate standard normal
     call random_number(u1)
     call random_number(u2)
     ztmp = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265 * u2)

     ! Scale to mean=1, stddev=0.1
     arr(i) = 1.0 + 0.1 * ztmp
  end do
  !------------------------------------------------------
  ! Simple initialization
  !------------------------------------------------------
  do i = 1, ncol
     do k = 1, nz

        cpair(i,k) = 1004.0_kind_phys
        rair(i,k)  = 287.0_kind_phys

        z(i,k)   = arr(i) * (100.0_kind_phys * real(k-1, kind_phys))
        rho(i,k) = arr(i) * (1.2_kind_phys * exp(-z(i,k)/8000.0_kind_phys))

        pk(i,k)    = arr(i) * (1.0_kind_phys)
        theta(i,k) = arr(i) * (300.0_kind_phys - 0.006_kind_phys*z(i,k))

        qv(i,k) = arr(i) * (0.010_kind_phys)
        qc(i,k) = arr(i) * (0.01_kind_phys)
        qr(i,k) = arr(i) * (0.01_kind_phys)

     end do
  end do

  !------------------------------------------------------
  ! Run microphysics
  !------------------------------------------------------
  call kessler_run(ncol, nz, dt, lyr_surf, lyr_toa, &
                   cpair, rair, rho, z, pk, &
                   theta, qv, qc, qr, &
                   precl, relhum, scheme_name, errmsg, errflg)

  if (errflg /= 0) then
     print *, 'Run error: ', trim(errmsg)
  else
     print *, 'Scheme name: ', trim(scheme_name)
     print *, 'theta: ', SUM(theta)
     print *, 'qv: ', SUM(qv)
     print *, 'qc: ', SUM(qc)
     print *, 'qr: ', SUM(qr)
     print *, 'Precip (m/s): ', SUM(precl)
     print *, 'relnum: ', SUM(relhum)
  end if

end program test_kessler_driver
