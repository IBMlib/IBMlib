ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Task driver to generate connectivity matrix connecting 
c     release boxes to settlement habitats
c     Settlement habitats are loaded and managed by particle_state
c     ---------------------------------------------------
c     $Rev: 147 $
c     $LastChangedDate: 2010-11-19 00:33:57 +0100 (Fri, 19 Nov 2010) $
c     $LastChangedBy: mpay $ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program tracker
      use input_parser
      use time_tools
      use run_context
      use physical_fields
      use particles
      use particle_tracking
      use particle_state     !  get_settlement_habitats
      use connectivity

      implicit none
c     ------------ declarations ------------      
      integer      :: idum4(4)
      type(clock),target  :: start_time, end_time
      type(clock),pointer :: current_time
      type(particle_ensemble)     :: par_ens
      type(emission_box), pointer :: emitboxes(:)
      type(state_attributes), pointer :: state_stack(:)
      type(tmatnc_file_handler)       :: ncfh 
      integer :: i, last, final_year
      integer :: idum(4), istep, nemax
      integer :: nsrc, ndest
      logical :: is_settled
      integer :: frombox, tobox
      real    :: time_step, xyz(3), survival
      real, allocatable    :: tmat(:,:), tmat_prior(:)
      integer, allocatable :: num_emit(:)
      character*256        :: filename

c     ------------   show time starts  ------------
      write(*,*) " >>>>>  connectivity run  <<<<<"

      call init_run_context()

c.....set clocks 
      call read_control_data(simulation_file, "start_time", idum4)
      call set_clock(start_time, idum4(1), idum4(2), idum4(3), idum4(4))
      call read_control_data(simulation_file, "end_time", idum4)
      call set_clock(end_time, idum4(1), idum4(2), idum4(3), idum4(4))
      final_year = idum4(1) ! store to data annotation
      call read_control_data(simulation_file, "particle_time_step", 
     +                       time_step)
     
      call init_physical_fields(start_time)
      call init_particles()
      call update_physical_fields()  ! need topology to verify emission boxes
      call setup_ensemble_from_file(par_ens, "emitbox", emitboxes)
      
      current_time => get_master_clock()
c     =====================  main time loop =====================
      istep   = 0
      last    = 0

      do while (compare_clocks(current_time, end_time) <= 0)
         write(*,372) istep
         call update_physical_fields()
         
c        -------- propagate tracers  --------
         
         call generate_particles(par_ens, emitboxes, time_step)
         call update_particles(par_ens, time_step)

c        -------- write tracer state -------- 
c        -------- loop control       --------
         call add_seconds_to_clock(current_time, nint(time_step))
         istep = istep + 1

         call get_last_particle_number(par_ens, last) 
         write(*,*) "number of free particles =", last

      enddo

 372  format(25("*"), " main time loop step ",i5.5, " ", 25("*"))

c
c     -------- generate connectivity from ensemble --------
c
      write(*,*) "generating connectivity matrix:"
      nsrc  = size(emitboxes)
      ndest = size(get_settlement_habitats())  ! particular particle_state method
      allocate( tmat(nsrc, ndest) )
      allocate( num_emit(nsrc)    )
      allocate( tmat_prior(nsrc)  )
      do i = 1, nsrc
         call get_current_emissions(emitboxes(i), num_emit(i), nemax)
         tmat_prior(i) = 1.0/ndest/max(1, num_emit(i))
      enddo
      state_stack => get_active_particle_states(par_ens)
      call get_connectivity_matrix(state_stack, num_emit,
     +                                   tmat_prior, tmat)


c$$$c     --- misc dumps: begin
c$$$      open(56, file="finalpos.dat")
c$$$      do i=1,last
c$$$c         write(*,'(72("-"))') 
c$$$c          call write_state_attributes(state_stack(i))
c$$$         call get_particle_position(get_particle(par_ens,i), xyz)
c$$$c          write(*,*) "position=", xyz
c$$$         call inquire_settling_success(state_stack(i), is_settled, 
c$$$     +                                 frombox, tobox, survival)  
c$$$         write(56,*) xyz, is_settled, frombox
c$$$      enddo
c$$$      close(56)
c$$$c     --- misc dumps:end
      

c
c     -------- dump connectivity to file as netCDF --------
c
      
      call read_control_data(simulation_file, 
     +                       "connectivity_matrix_file", filename)
      write(*,*) "writing connectivity matrix to file ", 
     +           trim(adjustl(filename))
      call tmatnc_start(filename, nsrc, ndest, ncfh)   ! initialize netCDF write of connectivity matrix
      call tmatnc_add_tmat(ncfh, tmat, 1.0*final_year)
      call tmatnc_close(ncfh)
c
c     ----------------- close down ---------------------------
c     
      call close_physical_fields()
      call close_particles()

      deallocate( tmat       )
      deallocate( num_emit   )
      deallocate( tmat_prior )

      write(*,*) "normal termination of simulation"

      end program
