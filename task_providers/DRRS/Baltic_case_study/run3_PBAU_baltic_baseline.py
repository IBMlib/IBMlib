#!/usr/bin/env python
###############################################################################################
#  Baseline 3 PBAU run (smagorinsky diffusivity)
#  
###############################################################################################
#
#  nohup nice -15 run_baseline_D1dot5.py  &
#
#  include stoke drift + wind from ECMWF
#
#  sudo mount -t nfs 10.128.22.171:/volume1/OPEC_20yr_hindcast  /mnt
#  sudo umount /mnt
#
#  grep reshuffled *log|awk '{print $2/1440.,$10}'|xmgrp
#  make clean; tar --exclude='./runs_*' -zcvf ../tarball.tgz *
#  make ibmrun; link
#  -------------------- HPC voodoo (LSF) --------------------
#  ssh asch@login.gbar.dtu.dk
#  scp ../tarball.tgz   asch@login.gbar.dtu.dk:CLAIM/WP1/IBMlib
#  scp run2b_PBAU_baltic_baseline.py asch@login.gbar.dtu.dk:CLAIM/WP1/IBMlib
#  scp asch@login.gbar.dtu.dk:CLAIM/WP1/IBMlib/*.nc  runs_15Jan2020
#  module load intel/2017.4.196
#  module load netcdf-fortran/4.4.4-intel-2017update4
#  FCFLAGS += -I/appl/netcdf-fortran/4.4.4/intel-2017update4/include 
#  LINKLIBS  += -L/appl/netcdf-fortran/4.4.4/intel-2017update4/lib
#
#  LSF voodoo: https://www.hpc.dtu.dk/?page_id=3356
#  bjobs        # check  Job status     
#  bkill JOB_ID # Delete a job
#  -------------------- old sparky -------------------- 
#  old sparky: 10.128.22.166
#  ssh asbjorn@10.128.22.166
#  scp ../tarball.tgz   asbjorn@10.128.22.166:DTU/CLAIM/Deliverables/D4.2_plastic_distribution_maps/IBMlib
#  scp  asbjorn@10.128.22.166:DTU/CLAIM/Deliverables/D4.2_plastic_distribution_maps/IBMlib/plast_distrib* .
#####################################################################################

from datetime import datetime, timedelta  # avoid importing time class, con
from numpy    import *   
import os
import sys
import string

(scriptstem, ext) = os.path.splitext(sys.argv[0])


simintval    = (2008,2012) # incl end year
#hydroDBroot  = "/mnt"
hydroDBroot  = "/work3/asch/HBM_physics"
npart        = 100000  # per source

# ===========================================================

batch_template = """#!/usr/bin/env bash
cd %s
module load intel/2017.4.196
module load netcdf-fortran/4.4.4-intel-2017update4
echo 'begin run'
./ibmrun %s 
echo 'end run'
"""

# ===========================================================

sim_template = """
! --------------------------------------------------------------
!   HBM hydrography
! --------------------------------------------------------------

hydroDBpath            = %(hydroDBroot)s   ! refers to a specific PBI (OPEC)
cmod_data_set          = coarse  ! value of X above
include_biogeochem     = F       ! T/F, optional, default T

vertical_diffusivity_scheme    = constant
vertical_diffusivity_constant  = 0.1  ! m2/sec
horizontal_diffusivity_scheme  = smagorinsky
smagorinsky_constant           = 0.2
smagorinsky_schmidtnum         = 1.0
smagorinsky_eddyvisc_max       = 200

wind_wave_height   = ../wave_BS_2008-2012_1h.nc  shww       ! file varname               
stokes_drift       = ../wave_BS_2008-2012_1h.nc  ust  vst   ! file uvarname vvarname
wind               = ../wind_BS_2008-2012_1h.nc  u10  v10   ! file uvarname vvarname
stokes_decay_depth = 4                                ! meters

! --------------------------------------------------------------
!   biology
! --------------------------------------------------------------
sink_rate  = %(sink_rate)f       ! per day
k_transfer = %(k_transfer)f      ! coalessed drag factor; unitless ; v = uv10*k_transfer

! -----------------------------------------------------------------------------------------------------------------
!   main
! -----------------------------------------------------------------------------------------------------------------
start_time     = %(year_start)d 01 01 30   ! format = YYYY MM DD S  where 0 < S < 86400 is second in the day
end_time       = %(year_end)d   12 31 0    ! format = YYYY MM DD S  where 0 < S < 86400 is second in the day

advec_intg_method  = euler    !  euler/rk2
particle_time_step = 1800     ! seconds

emitbox  = %(year_start)d 01 01 10  %(year_start)d 01 01 20     13.0000 53.8000 -0.00001   30.0000 65.74 -0.00001  %(npart)d

reshuf_memory_time = 30 ! unit = days
wet_point_file     = wet_13+.lonlat  ! allowed source points
source_point_file  = %(sources)s
open_boundary_line = 13 55 -1 0   ! x0,y0,vx0,vy0 : (x0,y0) = basept (vx0,vy0) = forbidden dir

outputfile           = %(outputfile)s
output_grid_dims     = 102 120
output_grid_SWcorner = 13.08333 53.75
output_grid_spacings = 0.166667 0.1
output_frequency     = 1   ! unit days

"""


for (srate, wpi, ktrans) in ((0.003,  10,  0.07),):
    # 
    jobtag      = "wpi_%dg_%d-%d_srate_%5.3f_hdiff=smag_ktrans_%4.2f" % (wpi, simintval[0], simintval[1], srate, ktrans)
    simfilename = "simpar_" + jobtag
    srcfilename = "all_input_%dg_13+.lonlat" % wpi
    job_params  = {"hydroDBroot" : hydroDBroot,
                   "sink_rate"   : srate,
                   "k_transfer"  : ktrans,
                   "year_start"  : simintval[0],
                   "year_end"    : simintval[1],
                   "npart"       : npart,
                   "outputfile"  : "PBAU_" + jobtag + ".nc",
                   "sources"     : srcfilename}
    #
    #  -------------------- generate simulation file: simfile--------------------
    #
    #
    f = open(simfilename, "w")
    f.write(sim_template % job_params)
    f.close()
    #
    # -------- run job --------
    #
    # os.system("ibmrun %s 2>&1" % simfilename)  
    #
    # -------- submit batch job --------
    jobfile_name = "run_%s" % jobtag
    jobfile      = open(jobfile_name, "w")
    jobfile.write(batch_template % (os.getcwd(), simfilename))
    jobfile.close()
    os.system("chmod +x %s"    % jobfile_name)
    os.system('bsub -n 1 -W 24:00 -R "rusage[mem=16GB]" -J PBAU -o %s.o < %s' % (jobfile_name, jobfile_name))  # system specific command
    #

