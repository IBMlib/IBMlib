ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Generate topography information from the BGM file
c     
c     Make scope public in module atlantis_grid before compiling
c     ---------------------------------------------------

c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c
c     Input parameters:
c          bgmfile               : (preprocessed) BGM file name defining Atlantis geometry
c          layer_widths          : [m] integer vector defining widths of vertical layers
c          samplings_per_layer   : how many times a layer i subdivided at vertical integral sampling
c          horiz_sampling_lenght : [m] lenght scale for generating horizontal integration meshes
c          start_time            : YYYY MM DD SSS  - start time of time frame genaration
c          end_time              : YYYY MM DD SSS  - end time of time frame genaration
c          dt_sampling           : [sec] minor time step to generate each time frame (dt_sampling < dt_frames)
c          fnametemp             : template to create file name; 'TAG' is replaced with <propertyname>
c
c     make ibmrun
c     ibmrun task_providers/Atlantis/simpar
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use atlantis_grid   ! NB: make default scope public

c     ------------ declarations ------------      
      implicit none
      type(clock),target  :: start_time
      character(len=999)  :: bgmfile, lwc
      integer, parameter  :: lwmax = 99  ! number of layers
      real                :: lwbuf(lwmax), drhoz
      integer             :: nsamppts, nlayers, idum4(4)
      integer             :: dt_sampling
      integer             :: ibox,ilay,iunits,ixy,atlay
      integer             :: nboxes,istat
      real                :: layer_width,cwd,cwdavg,area,a
      real                :: mwd, full_area
      character(len=3),parameter :: sep = "   " 
      type(AtlantisBox), pointer  :: box
      real, parameter     :: sediment_width = 0.5 ! meters

c     ------------   show time starts  ------------
      call init_run_context()
c.....    
      call read_control_data(simulation_file, "bgmfile", bgmfile)
      call read_control_data(simulation_file, "layer_widths", lwc)
      lwbuf = -1.0
      read(lwc,*,end=455) lwbuf
c.....detect the actual provided number of layers: nlayers
 455  do nlayers = 1, lwmax
         if (lwbuf(nlayers) < 0) exit
      enddo
      nlayers = nlayers - 1

      call read_control_data(simulation_file, "samplings_per_layer", 
     +                       nsamppts)
      call read_control_data(simulation_file, "horiz_sampling_lenght", 
     +                       drhoz)
      
c.....set clocks + load topography
      
      call read_control_data(simulation_file,"start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call init_physical_fields(start_time)
      call update_physical_fields()
      call initialize_atlantis_grid(bgmfile, lwbuf(1:nlayers), nsamppts, 
     +                              drhoz, .true.)

c     Volume per box, per depth layer
c     Please assume a sediment depth of 0.5 meters.
c     Depth layer thickness ("nominal dz")  - just for those boxes with depth >200, I know rest is fixed
c      # of wet layers per box

      call find_free_IO_unit(iunits) 
      open(iunits, file="box_info.dat")

      write(iunits, 946) "ibox", sep, "layer", sep, 
     +                   "volume", sep,
     +                   "nwet", sep, 
     +                   "lay_width", sep, 
     +                   "avg_depth"  

      nboxes  = size(Boxes)
c
c     ------- loop over all boxes -------
c
      write(*,533) "ibox", "maxwd", "avgwd", "fill"
      do ibox = 0, nboxes-1
         box    => Boxes(ibox)
         if (box%nwet == 0) cycle  ! no do process designated dry boxes
         area   = 0.0
         cwdavg = 0.0
         mwd    = 0.0
         do ixy = 1, size(box%xy, 2)  ! horizontal mesh loop
            call interpolate_wdepth(box%xy(:,ixy), cwd, istat)
            if (istat == 0) then !   assume same, no padding if /= 0 
               a      = box%da(ixy) ! area associated with this sampling, ignore slope
               area   = area + a    ! wet area
               cwdavg = cwdavg + cwd*a
               mwd    = max(mwd, cwd) ! max water depth in this box
            endif
         enddo
         cwdavg = cwdavg/area
         full_area = sum(box%da) ! of box, including dry points
         write(*,534) ibox, mwd, cwdavg, area/full_area

 533     format(a12,3x,3(a12,1x))
 534     format(i12,3x,3(f12.7,1x))

c        ------- loop over all layers for this box (wet+sediment)-------
         do atlay = 0, nlayers
            if ((atlay >= box%nwet).and.(atlay < nlayers)) cycle ! drop non-wet/non-sediment
            ilay = box%nwet - atlay ! corresponding internal layer index
            if (atlay == 0) then
               layer_width = cwdavg - zgrid_bounds(ilay)
            elseif (atlay == nlayers) then
               layer_width = sediment_width 
            else
               layer_width = zgrid_bounds(ilay+1)-zgrid_bounds(ilay)
            endif
            write(iunits, 947) ibox, sep, atlay, sep, 
     +                         full_area*layer_width, sep, 
     +                         box%nwet, sep,        
     +                         layer_width, sep, cwdavg
         enddo ! atlay
      enddo    ! ibox
      call close_atlantis_grid()
      close(iunits)

      write(*,*) "normal end of simulation"

 946  format(1x, a6, a, a6, a, a15,   a, a6, a, a8, a, a8)
 947  format(1x, i6, a, i6, a, e15.7, a, i6, a, f8.3, a, f8.3)    
      
      end program
