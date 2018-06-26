#!/usr/bin/python

import argparse
import re
#import os, sys
import glob
import numpy

#setup arguments to use in argparse (argument parsing)

parser = argparse.ArgumentParser(description='Average seperate files columnwize')
parser.add_argument('input', help="input file patern")
parser.add_argument('outfile', help="output file")
parser.add_argument("-s",help="Only sum, instead of average",action="store_true")
args = parser.parse_args()

print args

nfiles = 0
print "Search dir for pattern",args.input
for file in glob.glob(args.input):
    print "Processing",file
    data = numpy.loadtxt(file)
    if (nfiles) :
        avgdata += data # add to array
    else :
        avgdata = data # create array
    nfiles += 1

if (args.s):
    print "Saving sum"
else:
    print "Averaging"
    avgdata /= nfiles
    
print "Done, read in",nfiles,"files, writing outfile",args.outfile

numpy.savetxt(args.outfile,avgdata,fmt="%f")

