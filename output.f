ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Ouput module
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     TODO: At the moment, all classes are publicly visible. Scoping
c           should really be private. However, this is a major job
c           to change, and requires a possible rethinking of the
c           structure. Deferred to later dates
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      module output
      implicit none
      private

      type polytype
        private
        character (len=999) name
        character (len=999) :: type  !Type of data stored: real, int, log, str
        real:: r
        integer:: i
        !logical:: l
        character (len=999) :: s
      end type polytype
      public ::polytype

      type variable
        private
        character (len=999) name
        character (len=999) desc
        character (len=999) units
        character (len=999) fmt !Output format
        character (len=999) type !int, real, log, str, NF90_BYTE etc 
        real :: rescale(2)   !Scalar, bias
      end type variable
      public :: variable

      interface construct
        module procedure polytype_construct_real
        module procedure polytype_construct_integer
        module procedure polytype_construct_logical
        module procedure polytype_construct_string
        module procedure var_construct
      end interface
      public :: construct
      
      interface get_name
        module procedure var_get_name
        module procedure polytype_get_name
      end interface
      public :: get_name
      public :: get_desc
      public :: get_units
      public :: get_fmt
      interface get_type
        module procedure var_get_type
        module procedure polytype_get_type
      end interface
      public :: get_type
      public :: get_scalar
      public :: get_offset

      interface set_fmt
        module procedure set_fmt_single
        module procedure set_fmt_vec
      end interface
      public :: set_fmt
      interface set_type
        module procedure set_type_single
        module procedure set_type_vec
      end interface
      public :: set_type
      interface set_scalar
        module procedure set_scalar_single
        module procedure set_scalar_vec
      end interface
      public :: set_scalar
      interface set_offset
        module procedure set_offset_single
        module procedure set_offset_vec
      end interface
      public :: set_offset
       
      public :: get_data_string
      public :: get_data_int
      public :: get_data_real

c     ============================================================
      contains
c     ============================================================

c     ============================================================
c     Polytype class
c     ============================================================
      subroutine polytype_construct_real(this, name, rin)
        type(polytype),intent(inout) :: this
        character (len=*) :: name
        real,intent(in) :: rin
        !Assign arguments to real slot
        this%name = name
        this%type ="real"        
        this%r=rin
      end subroutine polytype_construct_real

      subroutine polytype_construct_integer(this, name, iin)
        type(polytype),intent(inout) :: this
        character (len=*) :: name
        integer,intent(in) :: iin
        !Assign arguments to integer slot
        this%name = name
        this%type ="int"
        this%i=iin
      end subroutine polytype_construct_integer

      subroutine polytype_construct_logical(this, name, lin)
        type(polytype),intent(inout) :: this
        character (len=*) :: name
        logical,intent(in) :: lin
        !Assign arguments to boolean slot
        this%name = name
        this%type ="int"
        if(lin) then
          this%i=1
        else 
          this%i =0
        endif
      end subroutine polytype_construct_logical

      subroutine polytype_construct_string(this, name, str)
        type(polytype),intent(inout) :: this
        character (len=*) :: name
        character (len=*),intent(in) :: str
        !Assign arguments to boolean slot
        this%name = name
        this%type ="str"
        this%s(1:len(str))=str
      end subroutine polytype_construct_string

      character (len=999) function polytype_get_name(this)
        type(polytype),intent(in) :: this
        polytype_get_name = this%name
      end function

      character (len=999) function polytype_get_type(this)
        type(polytype),intent(in) :: this
        polytype_get_type = this%type
      end function

      integer function get_data_int(this)
       type(polytype),intent(in) :: this
       get_data_int = this%i
      end function

      real function get_data_real(this)
       type(polytype),intent(in) :: this
       get_data_real = this%r
      end function

      subroutine  get_data_string(this,fmt,strout)
        type(polytype), intent(in) :: this
        character (*), intent(in) :: fmt
        character (*) :: strout
        !Format the string according to the data type
        select case(this%type)
        case ("real")
            write(strout,fmt=fmt) this%r
        case ("int")
            write(strout,fmt=trim(fmt)) this%i
        !case ("log")
        !    write(strout,fmt=trim(fmt)) this%l
        case ("str")
            write(strout,fmt=trim(fmt)) this%s
        case default
          write(*,*) "ERROR: Variable '"// trim(this%name) // 
     +       "' has an unknown type: "//trim(this%type)
          stop
        end select
      end subroutine 

c     ============================================================
c     Variable class
c     ============================================================
      subroutine var_construct(this,name,desc,units,fmt,type)
        type(variable),intent(out) :: this
        character(len=*) name, desc,units,fmt,type
        type(variable) :: tmp
        !Setup defaults
        this%name =name
        this%desc =desc
        this%units=units
        this%type =type
        this%fmt = fmt
        this%rescale = (/1.0, 0.0/)
      end subroutine

      character (len=999) function var_get_name(this)
        type(variable),intent(in) :: this
        var_get_name = this%name
      end function

      character (len=999) function get_desc(this)
        type(variable),intent(in) :: this
        get_desc = this%desc
      end function

      character (len=999) function get_fmt(this)
        type(variable),intent(in) :: this
        get_fmt = this%fmt
      end function

      subroutine set_fmt_single(this,fmt)
       type(variable), intent(inout) :: this
       character(*), intent(in) :: fmt
       this%fmt = fmt
      end subroutine

      subroutine set_fmt_vec(this,fmt)
       type(variable), intent(inout) :: this(:)
       character(*), intent(in) :: fmt
       integer i
       do i=1,size(this)
         call set_fmt_single(this(i), fmt)
       enddo
      end subroutine

      character (len=999) function get_units(this)
        type(variable),intent(in) :: this
        get_units = this%units
      end function

      character (len=999) function var_get_type(this)
        type(variable),intent(in) :: this
        var_get_type = this%type
      end function

      subroutine set_type_single(this, new_type)
        type(variable),intent(inout) :: this
        character(*), intent(in) :: new_type
        this%type = new_type
      end subroutine

      subroutine set_type_vec(this, new_type)
        type(variable),intent(inout) :: this(:)
        character(*), intent(in) :: new_type
        integer i
        do i=1,size(this)
          call set_type_single(this(i),new_type)
        enddo
      end subroutine

      real function get_scalar(this)
        type(variable),intent(in) :: this
        get_scalar = this%rescale(1)
      end function

      subroutine set_scalar_single(this, scalar)
        type(variable),intent(inout) :: this
        real, intent(in) :: scalar
        this%rescale(1) = scalar 
      end subroutine

      subroutine set_scalar_vec(this, scalar)
        type(variable),intent(inout) :: this(:)
        real, intent(in) :: scalar
        integer i
        do i=1,size(this)
          call set_scalar_single(this(i),scalar)
        enddo
      end subroutine

      real function get_offset(this)
        type(variable),intent(in) :: this
        get_offset = this%rescale(2)
      end function

      subroutine set_offset_single(this, offset)
        type(variable),intent(inout) :: this
        real, intent(in) :: offset
        this%rescale(2) = offset
      end subroutine

      subroutine set_offset_vec(this, offset)
        type(variable),intent(inout) :: this(:)
        real, intent(in) :: offset
        integer i
        do i=1,size(this)
           call set_offset_single(this(i),offset)
        enddo
      end subroutine

      end module output

