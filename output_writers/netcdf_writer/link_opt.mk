#Netcdf_writer link options

OUTPUT_MODS    += netcdf_writer.mod  

LINKLIBS_OUTPUT += -L/usr/lib  -L/usr/local/lib -lnetcdff -lnetcdf  -lhdf5_hl -lhdf5 -lz
