      program calc_tmat
      implicit none
c---------------------------------------------------------
c     compute tmat(i,j) = transport probability from i to j
c 
c     compile: ifort calc_tmat.f -o calc_tmat
c     usage  : calc_tmat <inputfile>
c---------------------------------------------------------
      integer, external :: iargc
      character*999     :: inputfile,program_name
      integer           :: nc,ic,id,nread,from,to,nproj
      real              :: finpos(2)
      real,allocatable    :: tmat(:,:), sw(:,:), ne(:,:),cc(:,:)
      real,allocatable    :: surv(:), sink(:), retent(:)
      integer,allocatable :: nemit(:)
      logical             :: in_any
c---------------------------------------------------------

c.....get program name and start time
      if (iargc() /= 1) then   ! no input file
         call getarg(0, program_name)
         write(*,*) "usage: ", trim(adjustl(program_name)),
     +              " <inputfile>"
         stop "syntax error"
      endif
      call getarg(1,inputfile)
      write(*,*) "reading data from file ", trim(adjustl(inputfile))
      open(22, file=inputfile)
      read(22,*) nc
      write(*,*) "contains ",nc,"cells"
      allocate(tmat(nc,nc))
      tmat=0.0
      allocate(sw(nc,2))
      allocate(ne(nc,2))
      allocate(cc(nc,2))
      allocate(nemit(nc))
      allocate(surv(nc))
      allocate(sink(nc))
      allocate(retent(nc))
c
c     read emission boxes
c
      nemit = 0
      do ic=1,nc    
         read(22,*) id, sw(id,:), ne(id,:), nemit(id)
         cc(id,:) = 0.5*(sw(id,:) + ne(id,:)) ! compute cell center
      enddo
      nread=0
      nproj=0
      write(*,*) "read ",nc,"cells"
c     
c     read particle positions
c
      do 
         read(22,*,end=777,err=779) from, finpos
         nread = nread+1
         call project_position(finpos,in_any,to)
         if (in_any) then
            tmat(from,to) = tmat(from,to) + 1
            nproj         = nproj + 1
         endif
      enddo
 779  stop "error reading particles"
 777  write(*,*) "read ", nread, "particles"
      write(*,*) nproj, "particles were inside cells"
      write(*,*) "unweighted avg survival = ",1.0*nproj/nread
c     
c     normalize transport matrix by emissions
c      
      do ic=1,nc    
         tmat(ic,:) = tmat(ic,:)/max(1,nemit(ic))
      enddo
      write(*,*) "normalized tmat"
      surv = sum(tmat,dim=2)
      sink = sum(tmat,dim=1)
      do ic=1,nc    
         retent(ic) = tmat(ic,ic)
      enddo
      write(*,*) minval(surv),  " < surv   < ", maxval(surv)
      write(*,*) minval(sink)  ," < sink   < ", maxval(sink)
      write(*,*) minval(retent)," < retent < ", maxval(retent)
      
c     -----------------------------------
                contains
c     -----------------------------------
      subroutine project_position(xy,inside,which)
      real,intent(in)     :: xy(:)
      logical,intent(out) :: inside
      integer,intent(out) :: which
      integer             :: i
      do i=1,nc
         if ((sw(i,1)<xy(1)) .and.
     +       (xy(1)<ne(i,1)) .and.
     +       (sw(i,2)<xy(2)) .and.
     +       (xy(2)<ne(i,2))) exit       
      enddo
      if (i == (nc+1)) then
         inside = .false.
      else
         inside = .true.
         which  = i
      endif 
      end subroutine project_position
      
      
      end program
      
      
