#!/usr/bin/env python
# -----------------------------------------------------------------
# ar extension with archiving merging capability
# -----------------------------------------------------------------
# $Rev:  $
# $LastChangedDate:  $
# $LastChangedBy: $ 
#
# Usage:
#    arm.py <archive_file>  {<archive_object>}+
#
# Some linkers (e.g. the ifort linker) is not able to extract 
# objects from an archive with multiple types, e.g. archives
# nested into an archives. This is a well-complianed issue
# in many fora without an elegant solution apart from explicit flattening.
#
# This archive wrapper will flatten a set of objects {<archive_object>}+
# by (recursively) inflating file objects with .a extension until
# all objects has a .o extension and then archive all objects into
# <archive_file>. Only last suffix is considered. It is assumed that
# all .a objects can be inflated by ar. Vertical name conflicts are
# not resolved (only last inflated kept), horizontal conflicts by
# archiving order.
#  
# Archive order: objects will be added in the order of the list
# {<archive_object>}+ Nested objects will be added in the order 
# of the nested archive(s) at the point of the parent archive, i.e
# archiving order is vertical first then horizontal.
#       
# Requires Python >= 2.3
#                      
# Options: currently no options
#    
# First edition: asc@11Mar2011
# -----------------------------------------------------------------
import sys
import os
import string 
import getopt
import tempfile 


def recursively_inflate(x):
    # -----------------------------------------
    # Function to recursively inflate content of 
    # archive x into tmpdir/objects
    # tmpdir is created locally to allow arbitrary recursion
    #
    # x is a path to an archive like a/b/c.a
    #
    # return a list of inflated objects (otherwise
    # archiving order can not be guarentied for the
    # final target archive) and the created tmpdir
    # -----------------------------------------
    (head, tail) = os.path.split(x)
    wkdir   = tempfile.mkdtemp(dir=head)
    content = os.popen("ar -t " + x).readlines()
    # -- prepend wkdir and strip trailing EOL
    for i in range(len(content)):
        content[i] = os.path.join(wkdir, content[i])[:-1]
    # -- inflate x
    inflate = "cd %s; ar -x ../%s" % (wkdir,tail)
    os.system(inflate)
    # -- recursively check inflated objects
    #    do not store lower temporary directories
    for i in range(len(content)):
        (root, ext) = os.path.splitext(content[i])
        if ext==".a": # replace archive in list by sectioning
            content[i:i+1],dummy = recursively_inflate(content[i])
    return content, wkdir


#
usage = "Usage: ", sys.argv[0] + " <archive_file>  {<archive_object>}+"
#
# ------------ split command ------------
#
try:
    opts, args = getopt.getopt(sys.argv[1:], "")
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    sys.exit()
#
# ------------ parse options ------------
#
#exclude_modules = []
#for o, a in opts:
#    if o == "-x": 
#        exclude_modules.append(a)
#    else:
#        assert False, ("unhandled option"+o)
#
# ------------ parse arguments   ------------
#
if len(args)<2: 
    print "no objects for archiving provided"
    print usage
    sys.exit()


#
# ------------ process sources  ------------
#
target         = args[0]
source_objects = args[1:]
archlist       = [] # ordered list of .o objects to archive
tmpdirs        = [] # only create it if it is needed 
for obj in source_objects:
    (root, ext) = os.path.splitext(obj) 
    if   ext==".o":
        archlist.append(obj) # need no further processing
    elif ext==".a":
        extracted_objects, dir = recursively_inflate(obj)
        archlist = archlist + extracted_objects
        tmpdirs.append(dir)
    else:
        print "do not know how to flatten object", obj, "with extension", ext
        sys.exit()
    
#
# --------- create flat archive  ------------
#
# remember to delete it first
#

archlist_as_string = " ".join(archlist)
os.system("rm -f %s; ar rcs %s %s" % (target, target, archlist_as_string))

#
# --- clean up temporary directories ---
#
for dir in tmpdirs:
    os.system("rm -rf %s" % dir)  
       
