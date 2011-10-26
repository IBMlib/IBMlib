ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Test a single state object in forward time simulation, 
c     at a fixed spatial position (no space update)
c     ---------------------------------------------------
c     $Rev: 147 $
c     $LastChangedDate: 2010-11-19 00:33:57 +0100 (Fri, 19 Nov 2010) $
c     $LastChangedBy: mpay $ 
c
c     make ibmrun
c     ibmrun task_providers/test/testpar
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program ibmrun
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particle_tracking
      use particle_state
      use output

      implicit none
c     ------------ declarations ------------      
      integer                  :: idum4(4),istat,istep, iunit
      type(clock),target       :: start_time, end_time
      type(clock),pointer      :: current_time
      type(state_attributes)   :: state
      type(spatial_attributes) :: space
      real                     :: r3(3),elapsed_days,time_step
      real                     :: lenght,weight,surv
      type(variable)           :: ask4len,ask4dw,ask4surv
      type(polytype)           :: bucket
c     ------------   show time starts  ------------
      call init_run_context()

c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "particle_time_step", 
     +                       time_step)

      call init_physical_fields(start_time)
      call init_particle_tracking()
      call init_particle_state()
      call update_physical_fields()
      call read_control_data(simulation_file, "position", r3)
      call set_tracer_position(space,r3) ! do not set remaining attributes, because it is static
      call init_state_attributes(state, space, 1.0, " ", 1)   ! invoke default init state 
 
      call construct(ask4len,  "std_length", "","","","")
      call construct(ask4dw,   "dry_weight", "","","","")
      call construct(ask4surv, "survival"  , "","","","") 

      call find_free_IO_unit(iunit)
      open(iunit,file="growth.trace")
c     =====================  main time loop =====================
      current_time => get_master_clock()
      istep   = 0
      do while (compare_clocks(current_time, end_time) <= 0)
         write(*,372) istep
         call update_particle_state(state, space, time_step)
c        -------- access/write tracer state -------- 
         call get_property(state,ask4len,  bucket,istat)
         if (istat /= 0) stop "state does not provide std_length"
         lenght = get_data_real(bucket)
         call get_property(state,ask4dw,   bucket,istat)
         if (istat /= 0) stop "state does not provide dry_weight"
         weight = get_data_real(bucket)
         call get_property(state,ask4surv, bucket,istat)
         if (istat /= 0) stop "state does not provide survival"
         surv   = get_data_real(bucket)
         elapsed_days = istep*time_step/86400.
         write(iunit,*) elapsed_days, lenght, weight, surv
c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         call update_physical_fields()
         istep = istep + 1
      enddo     
       
 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))
c     ----------------- close down ---------------------------
      call close_physical_fields()
      call close_particle_tracking()
      call close_particle_state()
      close(iunit)

      end program
