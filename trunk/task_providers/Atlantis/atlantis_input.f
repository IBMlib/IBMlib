ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Generate hydrodynamic input for Atlantis
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
c          dt_frames             : [sec] period between subsequent time frames
c          dt_sampling           : [sec] minor time step to generate each time frame (dt_sampling < dt_frames)
c          hydro_fname           : name of hydrodynamic fluxes 
c          temp_fname            : name of hydrodynamic fluxes 
c          salt_fname            : name of hydrodynamic fluxes 
c
c     make ibmrun
c     ibmrun task_providers/Atlantis/simpar
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use atlantis_grid
c     ------------ declarations ------------      
      implicit none
      type(clock),target  :: start_time, end_time
      character(len=999) :: bgmfile, lwc
      integer, parameter :: lwmax = 99  ! number of layers
      real               :: lwbuf(lwmax), drhoz
      integer            :: nsamppts, nlayers, idum4(4)
      integer            :: dt_frames, dt_sampling
      character(len=999) :: hydro_fname, temp_fname, salt_fname
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
      
      
c.....set clocks 
      
      call read_control_data(simulation_file,"start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file,"end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file,"dt_frames", dt_frames) 
      call read_control_data(simulation_file,"dt_sampling",dt_sampling)
      call read_control_data(simulation_file,"hydro_fname", hydro_fname)
      call read_control_data(simulation_file,"temp_fname", temp_fname)
      call read_control_data(simulation_file,"salt_fname", salt_fname)
      call init_physical_fields(start_time)
      call update_physical_fields()
      call initialize_atlantis_grid(bgmfile, lwbuf(1:nlayers), nsamppts, 
     +                              drhoz)

c.....currently set time offset to start time
      call make_input_files(start_time, end_time, 
     +                      dt_frames, dt_sampling,
     +                      start_time, hydro_fname,
     +                      temp_fname, salt_fname)

      call close_atlantis_grid()

      write(*,*) "normal end of simulation"
      
      end program
