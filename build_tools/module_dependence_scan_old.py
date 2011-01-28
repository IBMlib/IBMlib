#!/usr/bin/env python
# -----------------------------------------------------------------
# Quick'n'dirty scanning of module "use" associations of Fortran90
#
# Usage:
#   module_dependence_scan.py [f90file.f [ ...]]
#
# Scan for module "use" associations of Fortran90 in fortran source files provided
# as command arguments. Assume a module modulename is found in file modulename.f
#
# Scan for lines conforming to
#    use modulename1[[ ], ...]
#    use modulename2[[ ], ...]  etc
# Write a list:
#    f90file.o: modulename1.o modulename2.o
# etc. to stdout, i.e. in format usable for make dependences
# No output generated, if no module dependences are detected
# Match is case insensitive (as Fortran)
# Corresponding filename is generated in lowercase 
# (NB: there can be a compiler/platform issue here)
#                                   
# -----------------------------------------------------------------
import sys
import os
import string 
import re
#
# modext is the module code file extension used by the actual compiler
# NB: remember a trailing space after extension
#
modext = ".mod " # ifort 
#
use_assoc = re.compile("^[ ]*use[ ]*(?P<modname>[\w]+)[,\s]", re.IGNORECASE|re.MULTILINE)
for filename in sys.argv[1:]:
    (root, ext) = os.path.splitext(filename) # do not check extension
    sys.stdout.write(root + modext + ": ") 
    for line in open(filename).readlines():
        m = use_assoc.match(line)
        if m is not None:
            modname = string.lower(m.group('modname'))
            sys.stdout.write(modname+modext)
    sys.stdout.write("\n") 
#

