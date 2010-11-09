c     TODO: At the moment, all classes are publicly visible. Scoping
c           should really be private. However, this is a major job
c           to change, and requires a possible rethinking of the
c           structure. Deferred to later dates

      module output
      implicit none
      private

      type metadata
!        private      
        character (len=999) name
        character (len=999) desc
        character (len=999) units
        character (len=10)  default_fmt        !Suggested format specifier
        character (len=4) type    !Native representation of data
      end type metadata
      public:: metadata

      type polytype
!        private
        character (len=999) name
        character (len=4) :: type  !Type of data stored: real, int, log, str
        real:: r
        integer:: i
        logical:: l
        character (len=999) :: s
      end type polytype
      public ::polytype

      type output_var
!        private
        character (len=999)  name
        type(metadata) metadata
        character (len=20) fmt
        integer var_type    !eg short, double, integer, logical etc
        integer varid
        real ::  scalar  !No scaling by default
        real ::  offset
        type(polytype):: dat
      end type output_var
      public :: output_var

      interface construct
        module procedure polytype_construct_real
        module procedure polytype_construct_integer
        module procedure polytype_construct_logical
        module procedure polytype_construct_string
        module procedure output_var_construct
        module procedure metadata_construct
      end interface
      public :: construct
      
      interface get_name
        module procedure output_var_get_name
        module procedure polytype_get_name
        module procedure metadata_get_name
      end interface
      public :: get_name

      interface store
        module procedure polytype_store
        module procedure metadata_store
      end interface
      public :: store
c     ============================================================
      contains
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
        this%type ="log"
        this%l=lin
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


      subroutine output_var_construct(this,name,fmt)
        type(output_var)    :: this
        character(*) :: name
        character (*), optional:: fmt
        !Setup defaults
        this%name =name
        this%fmt=""
        this%scalar=1.0
        this%offset=0.0
        !Put all optional arguments into slots
        if(present(fmt)) this%fmt=fmt
      end subroutine

      subroutine metadata_construct(this,name,desc,units,fmt,type)
        type(metadata),intent(inout) :: this
        character(len=*) name, desc,units,fmt,type
        this=metadata(name,desc,units,fmt,type)
      end subroutine metadata_construct
      
      
      character (len=999) function output_var_get_name(this)
        type(output_var),intent(in) :: this
        output_var_get_name = this%name
      end function

      character (len=999) function polytype_get_name(this)
        type(polytype),intent(in) :: this
        polytype_get_name = this%name
      end function

      character (len=999) function metadata_get_name(this)
        type(metadata),intent(in) :: this
        metadata_get_name = this%name
      end function

      subroutine polytype_store(this,bucket)
        type(output_var), intent(inout) :: this
        type(polytype), intent(in) :: bucket
        this%dat = bucket
      end subroutine

      subroutine metadata_store(this,meta)
        type(output_var), intent(inout) :: this
        type(metadata), intent(in) :: meta
        this%metadata = meta
      end subroutine

      end module output
