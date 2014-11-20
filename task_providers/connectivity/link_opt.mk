#
#  --- exported link options to final build ---
#
LINKFLAGS_TASK = -I/usr/local/include
#LINKLIBS_TASK  = -L/usr/lib  -L/usr/local/lib -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -lz
# ASCH 02Oct2014: avoid /usr/lib first, prefer netcdf in /usr/local/lib 
LINKLIBS_TASK  =   -L/usr/local/lib -L/usr/lib -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -lz
