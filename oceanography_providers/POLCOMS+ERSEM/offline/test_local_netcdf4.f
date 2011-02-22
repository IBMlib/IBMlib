      program jj
      use netcdf
      implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ifort -I/usr/local/include  test_local_netcdf4.f -L/usr/local/lib -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -lz
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer :: ncid, varid    ! This will be the netCDF ID for the file and data variable.
      real    :: lat(224)
      character*200::fname="data/Monthly.PolcomsErsem.2000.01.HDF.gz.nc"
c     Open the file. NF90_NOWRITE tells netCDF we want read-only access to
c     the file.
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid) )

c     Get the varid of the data variable, based on its name.
      call NetCDFcheck( nf90_inq_varid(ncid, "lat", varid) )

c     Read the data.
      call NetCDFcheck( nf90_get_var(ncid, varid, lat) )

      write(*,*) lat
c     Close the file, freeing all resources.
      call NetCDFcheck( nf90_close(ncid) )

      end program 


      subroutine NetCDFcheck(status)
c     ------------------------------------------
c     This subroutine supports the recommended reading 
c     style of NetCDF4:
c       call NetCDFcheck( nf90_get_var(ncid, varid, data_in) )
c     Private to this module
c     ------------------------------------------
      use netcdf
      integer, intent ( in) :: status
      
      if(status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         stop "NetCDFcheck:Stopped"
      end if
      end subroutine NetCDFcheck 
