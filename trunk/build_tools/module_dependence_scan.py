#!/usr/bin/env python
# -----------------------------------------------------------------
#  Module dependencies scan 
# ---------------------------------------------
# $Rev: 147 $
# $LastChangedDate: 2010-11-19 00:33:57 +0100 (Fri, 19 Nov 2010) $
# $LastChangedBy: mpay $ 
#
# Quick'n'dirty scanning of module "use" associations of Fortran90
#
# Usage:
#   module_dependence_scan.py [options] [f90filename*]
#
# Scan for module "use" associations of Fortran90 in fortran source files provided
# as command arguments. Assume a module modulename is found in file modulename.f
#
# Scan for lines conforming to
#    use modulename1[[ ], ...]
#    use modulename2[[ ], ...]  etc
# Write a list:
#    f90file.mod: modulename1.mod modulename2.mod
#    f90file.o:   modulename1.mod modulename2.mod
#    .... 
# etc. to stdout, i.e. in format usable for make dependences
# No output generated, if no module dependences are detected
# Match is case insensitive (as Fortran)
# Corresponding filename is generated in lowercase 
# (NB: there can be a compiler/platform issue here)
#                                   
# Options: 
#    -x modulename : exclude modulename from dependency list 
#                    (typically contains system resources)
#                    NB: modulename should be without extension !
#
# Extended version 2 / Jun 16, 2010
# -----------------------------------------------------------------
import sys
import os
import string 
import re
import getopt
#
# modext is the module code file extension used by the actual compiler
# without leading/trailing spaces
#
modext = ".mod" # ifort 
objext = ".o"   # all ? 
#
use_assoc = re.compile("^[ ]*use[ ]*(?P<modname>[\w]+)[,\s]", re.IGNORECASE|re.MULTILINE)
#
# ------------ split command ------------
#
try:
    opts, args = getopt.getopt(sys.argv[1:], "x:")
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    sys.exit()
#
# ------------ parse options ------------
#
exclude_modules = []
for o, a in opts:
    if o == "-x": 
        exclude_modules.append(a)
    else:
        assert False, ("unhandled option"+o)
#
# ------------ parse arguments (assumed F90 files)  ------------
#
for filename in args:
    (root, ext) = os.path.splitext(filename) # do not check extension
    dependencies = ""
    # --- read and parse F90 file line-by-line ---
    for line in open(filename).readlines():
        m = use_assoc.match(line)
        if m is not None:
            modname = string.lower(m.group('modname'))
            if not (modname in exclude_modules):
                dependencies = dependencies + modname + modext + " "
    # --- flush result for this F90 file, if any dependencies
    #     write nothing if no dependencies (to avoid a phony target)
    if len(dependencies)>0:
        sys.stdout.write(root + modext + " : " + dependencies + "\n")
        sys.stdout.write(root + objext + " : " + dependencies + "\n") 
#



