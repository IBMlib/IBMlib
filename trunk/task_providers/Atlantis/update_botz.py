#!/usr/bin/env python
# -*- coding: utf-8 -*-  
##########################################################################################
#
#  Update botz entries in BGM file
#  usage: cat BGMfile|update_botz.py new_botz.dat > updated_BGMfile
##########################################################################################

from string import *
from sys    import *

botzdata = open(argv[1], "r").readlines()
replace_list = {}
for line in botzdata:
    items = split(line)
    replace_list["box%s.botz" % items[0]] = "-"+items[1]
maxwcbotz = min(map(float, replace_list.values())) # maxwcbotz < 0
replace_list["maxwcbotz"] = "%f" % maxwcbotz
#
for line in stdin.readlines():
    items = split(line)
    pline = line
    if len(items)>0:
        for tag in replace_list.keys():
            if items[0] == tag: pline = tag + " " + replace_list[tag]
    print pline[:-1] # skip trailing EOL
            
