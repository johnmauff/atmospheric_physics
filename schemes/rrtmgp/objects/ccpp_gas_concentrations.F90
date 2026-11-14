module ccpp_gas_concentrations
   ! CCPP wrapper for ty_gas_concs DDT from RRTMGP
   !
   ! RRTMGP records the number of columns and layers of the gas
   ! concentrations privately (from the first array given to set_vmr) and
   ! offers no way to read them back. The wrapper keeps a copy, updated by
   ! its own init/set_vmr/reset, which forward to RRTMGP, so callers can
   ! size a get_vmr buffer from the object itself (get_ncol, get_nlay).
   ! The copy is right only while every init, set_vmr and reset goes
   ! through the wrapper; read-only uses (gas optics) pass gas_concs as
   ! before. 0 means no gas has that dimension (set as a scalar or a
   ! profile), as in RRTMGP.
   use mo_rte_kind,           only: wp
   use mo_gas_concentrations, only: ty_gas_concs

   implicit none
   private

   !> \section arg_table_ty_gas_concs_ccpp Argument Table
   !! \htmlinclude ty_gas_concs_ccpp.html
   type, public :: ty_gas_concs_ccpp
      type(ty_gas_concs) :: gas_concs
      integer, private   :: ncol = 0
      integer, private   :: nlay = 0
   contains
      procedure, public  :: init
      procedure, public  :: reset
      procedure, private :: set_vmr_scalar
      procedure, private :: set_vmr_1d
      procedure, private :: set_vmr_2d
      generic,   public  :: set_vmr => set_vmr_scalar, set_vmr_1d, set_vmr_2d
      procedure, public  :: get_ncol
      procedure, public  :: get_nlay
   end type

contains

   function init(this, gas_names) result(error_msg)
      class(ty_gas_concs_ccpp),       intent(inout) :: this
      character(len=*), dimension(:), intent(in)    :: gas_names
      character(len=128)                            :: error_msg

      ! RRTMGP's init resets the object, its shape included, unless it
      ! rejects the names first
      error_msg = this%gas_concs%init(gas_names)
      if (len_trim(error_msg) == 0) then
         this%ncol = 0
         this%nlay = 0
      end if
   end function init

   subroutine reset(this)
      class(ty_gas_concs_ccpp), intent(inout) :: this

      call this%gas_concs%reset()
      this%ncol = 0
      this%nlay = 0
   end subroutine reset

   function set_vmr_scalar(this, gas, w) result(error_msg)
      class(ty_gas_concs_ccpp), intent(inout) :: this
      character(len=*),         intent(in)    :: gas
      real(wp),                 intent(in)    :: w
      character(len=128)                      :: error_msg

      ! A scalar sets no dimension
      error_msg = this%gas_concs%set_vmr(gas, w)
   end function set_vmr_scalar

   function set_vmr_1d(this, gas, w) result(error_msg)
      class(ty_gas_concs_ccpp), intent(inout) :: this
      character(len=*),         intent(in)    :: gas
      real(wp), dimension(:),   intent(in)    :: w
      character(len=128)                      :: error_msg

      ! RRTMGP rejects a profile whose size differs from nlay, once set
      error_msg = this%gas_concs%set_vmr(gas, w)
      if (len_trim(error_msg) == 0) then
         this%nlay = size(w)
      end if
   end function set_vmr_1d

   function set_vmr_2d(this, gas, w) result(error_msg)
      class(ty_gas_concs_ccpp), intent(inout) :: this
      character(len=*),         intent(in)    :: gas
      real(wp), dimension(:,:), intent(in)    :: w
      character(len=128)                      :: error_msg

      ! RRTMGP rejects an array whose shape differs from ncol/nlay, once set
      error_msg = this%gas_concs%set_vmr(gas, w)
      if (len_trim(error_msg) == 0) then
         this%ncol = size(w, 1)
         this%nlay = size(w, 2)
      end if
   end function set_vmr_2d

   pure integer function get_ncol(this)
      class(ty_gas_concs_ccpp), intent(in) :: this

      get_ncol = this%ncol
   end function get_ncol

   pure integer function get_nlay(this)
      class(ty_gas_concs_ccpp), intent(in) :: this

      get_nlay = this%nlay
   end function get_nlay

end module ccpp_gas_concentrations
