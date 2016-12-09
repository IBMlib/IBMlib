      subroutine resolve_corresp_dataset(aclock, tag)
c     ------------------------------------------
c     Resolve tag of data set corresponding to aclock
c
c     tag in ANONYMOUS for daily averaged data sets is YYYY_MM_DD
c     In test mode, a fixed frame can be recycled by setting
c     fixed_hydro_frame = YYYY_MM_DD in simulation file
c     ------------------------------------------
      type(clock), intent(in)                :: aclock
      character(len=tag_lenght), intent(out) :: tag
      integer                                :: year, month, day
c     ------------------------------------------
      if (fixed_tag == not_set_c) then
         call get_date_from_clock(aclock, year, month, day)
         write(tag, 427) year, month, day
 427     format(i4.4,"_",i2.2,"_",i2.2)           !  YYYY_MM_DD
      else
         tag = fixed_tag
         write(*,*) "resolve_corresp_dataset: applying fixed_tag:", tag
      endif

      end subroutine resolve_corresp_dataset
