#!/usr/bin/python

import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
import math

#setup arguments to use in argparse (argument parsing)

parser = argparse.ArgumentParser(description='Get seperate files and total dipoles for each mol from dipole.out output (time mux muy muz tdp')
parser.add_argument('input', help="input file")
parser.add_argument('-ts',type=float,default=-1,help="time step i.e dt")
parser.add_argument("-td",help="Do total dipole only",action="store_false")
args = parser.parse_args()

print args

fp = open(args.input)

# read header
while(1):
    line = fp.readline()
    if not line:
        break
    if not (line[0] == '#'):
        break

nmol = int(line.split()[1])
step = float(line.split()[0])
print "read in header, looping over configs, with",nmol,"molecules"
print "Opening files"
fo = []
for i in range(nmol):
    name = args.input + "." + str(i+1)
    fo.append(open(name,"w"))

name = args.input + ".t"
ft = open(name,"w")

nconf = 0;
while (line):
    dtot = [0,0,0]
    for i in range(nmol):
        line = fp.readline()
        if not line:
            print "Error reading in, config # ",nconf
            exit(1)
        dat = line.split()
        print dat
        fo[i].write(str(step) + " " + dat[1] + " " + dat[2] + " " + dat[3] + "\n")
        dtot[0] += float(dat[1])
        dtot[1] += float(dat[2])
        dtot[2] += float(dat[3])
    nconf += 1

    # write total file
    ft.write(str(step) + " " + str(dtot[0]) + " " + str(dtot[1]) + " " + str(dtot[2]) + "\n")
        
    line = fp.readline() # read header of next config
    if line: # error check header of next config
        if not (int(line.split()[1]) == nmol):
            print "Error number of atoms changed! at config",nconf
            exit(1)
        step = float(line.split()[0])

print dat
print "Done, read in",nconf,"configs"

