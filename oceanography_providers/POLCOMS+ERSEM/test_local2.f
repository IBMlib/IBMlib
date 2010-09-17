      program jj
      use netcdf
      implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ifort -I/usr/local/include  test_local2.f -L/usr/local/lib -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -lz
c     
c     NS test point : (ilon,ilat) == (150,130)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer, parameter ::  z = 40 ;
      integer, parameter ::  n_zfaces = 2
      integer, parameter ::  time = 1
      integer, parameter ::  lon = 198 
      integer, parameter ::  lat = 224 
      character*200::fname="data/Monthly.PolcomsErsem.2000.01.HDF.gz.nc"

      real :: ucurM(lon, lat, z,time)
      real :: vcurM(lon, lat, z,time)
      real :: wcurM(lon, lat, z,time)
      real :: ETWM(lon, lat, z,time)

      real :: zbnd(n_zfaces,lon, lat, z,time)
      real :: bathymetry(lon,lat)
      real :: depth(lon, lat, z,time)
      real :: pdepthM(lon, lat, z,time)

      
      integer :: ncid, varid    ! This will be the netCDF ID for the file and data variable.
      integer :: iz,ilon,ilat,iface
      
c
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid) )
c
      call NetCDFcheck( nf90_inq_varid(ncid, "ucurM", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, ucurM) )
c
      call NetCDFcheck( nf90_inq_varid(ncid, "vcurM", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, vcurM) )
c
      call NetCDFcheck( nf90_inq_varid(ncid, "wcurM", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, wcurM) )
c
      call NetCDFcheck( nf90_inq_varid(ncid, "ETWM", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, ETWM) )
c
      call NetCDFcheck( nf90_inq_varid(ncid, "zbnd", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, zbnd) )
c
      call NetCDFcheck( nf90_inq_varid(ncid, "bathymetry", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, bathymetry) )
c
      call NetCDFcheck( nf90_inq_varid(ncid, "depth", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, depth) )
c
      call NetCDFcheck( nf90_inq_varid(ncid, "pdepthM", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, pdepthM) )

      call NetCDFcheck( nf90_close(ncid) )
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c$$$c     ----- geo map -----
c$$$      do ilon=1,lon
c$$$         do ilat=1,lat
c$$$            if (abs(bathymetry(ilon,ilat)<0.99*nf90_fill_float)) then
c$$$               write(*,*) ilon,ilat
c$$$            endif
c$$$         enddo
c$$$      enddo

      do ilon=1,lon
         do ilat=1,lat
            write(*,*) zbnd(2,ilon,ilat,z,1)
         enddo
      enddo

c$$$c     NS point == (150,130)
c$$$c      ilon = 150; ilat = 130
c$$$      ilon = 1; ilat = 1
c$$$      write(*,*) "bathymetry=", bathymetry(ilon,ilat)
c$$$      write(*,*) "  iz   zbnd(1:2)    depth   pdepthM"
c$$$      do iz=1,z
c$$$         write(*,*) iz, zbnd(1:2,ilon,ilat,iz,1), depth(ilon,ilat,iz,1),
c$$$     +                  pdepthM(ilon,ilat,iz,1)
c$$$      enddo

      

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
