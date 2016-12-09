      subroutine resolve_corresp_dataset(aclock, tag)
c     ------------------------------------------
c     Resolve tag of data set corresponding to aclock
c
c     tag in BIMS-ECO PBI for daily averaged data sets is DDDD,
c     where DDDD is days since 1990-01-01
c     In test mode, a fixed frame can be recycled by setting
c     fixed_hydro_frame = DDDD in simulation file
c     ------------------------------------------
      type(clock), intent(in)                :: aclock
      character(len=tag_lenght), intent(out) :: tag
      integer                                :: nhours, ndays
c     ------------------------------------------
      if (fixed_tag == not_set_c) then
         call get_period_length_hour(time_offset, aclock, nhours) ! round to nearest hour
         ndays = int(nhours/24.0) 
         write(tag, 427) ndays 
 427     format(i4.4)           !  DDDD
      else
         tag = fixed_tag
         write(*,*) "resolve_corresp_dataset: applying fixed_tag:", tag
      endif

      end subroutine resolve_corresp_dataset
