module dadadj
  !=======================================================================
  ! GFDL style dry adiabatic adjustment
  !
  ! Method:
  ! if stratification is unstable, adjustment to the dry adiabatic lapse
  ! rate is forced subject to the condition that enthalpy is conserved.
  !=======================================================================

  use ccpp_kinds, only:  kind_phys

  implicit none
  private
  save

  public :: dadadj_init ! init routine
  public :: dadadj_run  ! main routine

  integer  :: nlvdry  ! number of layers from top of model to apply the adjustment
  integer  :: niter   ! number of iterations for convergence

CONTAINS

  !> \section arg_table_dadadj_init Argument Table
  !! \htmlinclude dadadj_init.html
  subroutine dadadj_init(dadadj_nlvdry, dadadj_niter, nz, errmsg, errflg)
    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------
    integer,            intent(in)  :: dadadj_nlvdry
    integer,            intent(in)  :: dadadj_niter
    integer,            intent(in)  :: nz
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    errmsg = ''
    errflg = 0

    if (dadadj_nlvdry >= nz .or. dadadj_nlvdry < 0) then
       errflg = 1
       write(errmsg,*) 'dadadj_init: dadadj_nlvdry=',dadadj_nlvdry,' but must be less than the number of vertical levels ',&
       '(',nz,'), and must be a positive integer.'
    end if

    nlvdry = dadadj_nlvdry
    niter = dadadj_niter

  end subroutine dadadj_init

  !> \section arg_table_dadadj_run Argument Table
  !! \htmlinclude dadadj_run.html
  subroutine dadadj_run( &
       ncol, nz, dt, pmid, pint, pdel, state_t, state_q, cappa, cpair, s_tend, &
       q_tend, dadpdf, scheme_name, errmsg, errflg)

    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------
    integer,  intent(in)        :: ncol  ! number of atmospheric columns
    integer,  intent(in)        :: nz    ! number of atmospheric levels
    real(kind_phys), intent(in) :: dt          ! physics timestep
    real(kind_phys), intent(in) :: pmid(:,:)   ! pressure at model levels
    real(kind_phys), intent(in) :: pint(:,:)   ! pressure at model interfaces
    real(kind_phys), intent(in) :: pdel(:,:)   ! vertical delta-p
    real(kind_phys), intent(in) :: cappa(:,:)  ! variable Kappa
    real(kind_phys), intent(in) :: cpair(:,:)  ! heat capacity of air
    real(kind_phys), intent(in) :: state_t(:,:)   ! temperature (K)
    real(kind_phys), intent(in) :: state_q(:,:)   ! specific humidity
    real(kind_phys), intent(out) :: s_tend(:,:)   ! dry air enthalpy tendency
    real(kind_phys), intent(out) :: q_tend(:,:)   ! specific humidity tendency
    real(kind_phys), intent(out) :: dadpdf(:,:)  ! PDF of where adjustments happened

    character(len=64),  intent(out) :: scheme_name
    character(len=512), intent(out) :: errmsg
    integer,            intent(out) :: errflg

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------

    integer                      :: i,k      ! longitude, level indices
    integer                      :: jiter    ! iteration index
    integer                      :: nlv      ! nlvdry, for the device
    integer                      :: nit      ! niter, for the device
    ! Intermediate constants, per column (each column is adjusted on its own).
    real(kind_phys), allocatable :: c1dad(:,:)
    real(kind_phys), allocatable :: c2dad(:,:)
    real(kind_phys), allocatable :: c3dad(:,:)
    real(kind_phys), allocatable :: c4dad(:,:)
    real(kind_phys) :: gammad    ! dry adiabatic lapse rate (deg/Pa)
    real(kind_phys) :: zeps      ! convergence criterion (deg/Pa)
    real(kind_phys) :: rdenom    ! reciprocal of denominator of expression
    real(kind_phys) :: dtdp      ! delta-t/delta-p
    real(kind_phys) :: zepsdp    ! zeps*delta-p
    real(kind_phys) :: zgamma    ! intermediate constant
    real(kind_phys) :: qave      ! mean q between levels
    real(kind_phys) :: cappaint  ! Kappa at level intefaces
    real(kind_phys) :: t(ncol,nz)
    real(kind_phys) :: q(ncol,nz)

    logical :: ilconv          ! .TRUE. ==> convergence was attained
    logical :: done            ! .TRUE. ==> the column is finished
    logical :: dodad           ! .TRUE. ==> do dry adjustment
    logical :: failed(ncol)    ! .TRUE. ==> no convergence in the column

    !-----------------------------------------------------------------------

    errmsg = ''
    errflg = 0
    scheme_name = 'DADADJ'
    nlv = nlvdry
    nit = niter

    allocate(c1dad(ncol,nlvdry), stat=errflg)
    if (errflg /= 0) then
       errmsg = trim(scheme_name)//': Allocate of c1dad(ncol,nlvdry) failed'
       return
    end if
    allocate(c2dad(ncol,nlvdry), stat=errflg)
    if (errflg /= 0) then
       errmsg = trim(scheme_name)//': Allocate of c2dad(ncol,nlvdry) failed'
       return
    end if
    allocate(c3dad(ncol,nlvdry), stat=errflg)
    if (errflg /= 0) then
       errmsg = trim(scheme_name)//': Allocate of c3dad(ncol,nlvdry) failed'
       return
    end if
    allocate(c4dad(ncol,nlvdry), stat=errflg)
    if (errflg /= 0) then
       errmsg = trim(scheme_name)//': Allocate of c4dad(ncol,nlvdry) failed'
       return
    end if

    ! Columns are independent: one pass per column finds unstable
    ! stratification, makes the dry adiabatic adjustment and computes the
    ! tendencies. Note: nlvdry ****MUST**** be < pver
    !$acc parallel loop gang vector &
    !$acc   copyin(pmid, pint, pdel, cappa, cpair, state_t, state_q) &
    !$acc   copyout(s_tend, q_tend, dadpdf, failed) &
    !$acc   create(t, q, c1dad, c2dad, c3dad, c4dad)
    COL: do i = 1, ncol

       do k = 1, nz
          t(i,k) = state_t(i,k)
          q(i,k) = state_q(i,k)
          dadpdf(i,k) = 0._kind_phys
       end do
       failed(i) = .false.

       ! Find gridpoints with unstable stratification
       zeps = 2.0e-5_kind_phys           ! set convergence criteria
       cappaint = 0.5_kind_phys*(cappa(i,2) + cappa(i,1))
       gammad = cappaint*0.5_kind_phys*(t(i,2) + t(i,1))/pint(i,2)
       dtdp = (t(i,2) - t(i,1))/(pmid(i,2) - pmid(i,1))
       dodad = (dtdp + zeps) > gammad
       do k = 2, nlv
          cappaint = 0.5_kind_phys*(cappa(i,k+1) + cappa(i,k))
          gammad = cappaint*0.5_kind_phys*(t(i,k+1) + t(i,k))/pint(i,k+1)
          dtdp = (t(i,k+1) - t(i,k))/(pmid(i,k+1) - pmid(i,k))
          dodad = dodad .or. (dtdp + zeps) > gammad
          if ((dtdp + zeps) > gammad) then
             dadpdf(i,k) = 1._kind_phys
          end if
       end do

       ! Make a dry adiabatic adjustment
       if (dodad) then

          do k = 1, nlv
             cappaint = 0.5_kind_phys*(cappa(i,k+1) + cappa(i,k))
             c1dad(i,k) = cappaint*0.5_kind_phys*(pmid(i,k+1)-pmid(i,k))/pint(i,k+1)
             c2dad(i,k) = (1._kind_phys - c1dad(i,k))/(1._kind_phys + c1dad(i,k))
             rdenom = 1._kind_phys/(pdel(i,k)*c2dad(i,k) + pdel(i,k+1))
             c3dad(i,k) = rdenom*pdel(i,k)
             c4dad(i,k) = rdenom*pdel(i,k+1)
          end do

          ! Iterate to convergence; if niter iterations do not converge,
          ! double zeps and start again, until zeps exceeds 1.e-4.
          done = .false.
          DBLZEP: do while (.not. done)

             ilconv = .false.
             do jiter = 1, nit
                ilconv = .true.

                do k = 1, nlv
                   zepsdp = zeps*(pmid(i,k+1) - pmid(i,k))
                   zgamma = c1dad(i,k)*(t(i,k) + t(i,k+1))

                   if ((t(i,k+1)-t(i,k)) >= (zgamma+zepsdp)) then
                      ilconv = .false.
                      t(i,k+1) = t(i,k)*c3dad(i,k) + t(i,k+1)*c4dad(i,k)
                      t(i,k) = c2dad(i,k)*t(i,k+1)
                      qave = (pdel(i,k+1)*q(i,k+1) + pdel(i,k)*q(i,k))/(pdel(i,k+1)+ pdel(i,k))
                      q(i,k+1) = qave
                      q(i,k) = qave
                   end if

                end do

                if (ilconv) exit ! convergence => this column is done

             end do

             if (ilconv) then
                done = .true.
             else
                zeps = zeps + zeps
                if (zeps > 1.e-4_kind_phys) then
                   failed(i) = .true.
                   done = .true.
                end if
             end if
          end do DBLZEP

       end if

       do k = 1, nz
          s_tend(i,k) = (t(i,k) - state_t(i,k))/dt*cpair(i,k)
          q_tend(i,k) = (q(i,k) - state_q(i,k))/dt
       end do

    end do COL

    deallocate(c1dad, c2dad, c3dad, c4dad)

    ! Error return at the first column that did not converge.
    do i = 1, ncol
       if (failed(i)) then
          errflg = i
          write(errmsg,*) trim(scheme_name)//': Convergence failure at column ',i,' zeps > 1.e-4 '// &
               '(errflg set to failing column index)'
          return
       end if
    end do

 end subroutine dadadj_run

end module dadadj
