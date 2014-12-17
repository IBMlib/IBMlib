#
#  --- exported link options to final build ---
#
LINKFLAGS_PHYSICAL = -I/usr/local/include
# LINKLIBS_PHYSICAL  = -L/usr/lib  -L/usr/local/lib -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -lz
# ASCH 25Nov2014: avoid /usr/lib first, prefer netcdf in /usr/local/lib 
LINKLIBS_TASK  =   -L/usr/local/lib -L/usr/lib -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -lz
