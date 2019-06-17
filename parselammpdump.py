#!/usr/bin/python

import argparse

#setup arguments to use in argparse (argument parsing)

parser = argparse.ArgumentParser(description='Parse lammps custum dump id mol type x y z vx vy vz')
parser.add_argument('input', help="input file")
args = parser.parse_args()

print args

fp = open(args.input)

# read header 

line = fp.readline() # ITEM: TIMESTEP
line = fp.readline() # ts
line = fp.readline() # ITEM: NUMBER OF ATOMS
line = fp.readline() # natoms
natoms = int(line.split()[0])
print "natoms =",natoms
line = fp.readline() # ITEM: BOX BOUNDS pp pp pp
line = fp.readline() # xlo xhi
line = fp.readline() # ylo yhi
line = fp.readline() # zlo zhi
line = fp.readline() # ITEM: ATOMS id mol type x y z vx vy vz

atoms = {} # use dictonaries, with atom id as keys
    
for i in range(natoms):
    line = fp.readline() # ITEM: ATOMS id mol type x y z vx vy vz
    id = int(line.split()[0])
    atoms[id] = line

for i in range(1,natoms+1):
    if not (i in atoms):
        print "key: ",i, "not found"
        exit(1)
exit(1)

print "Read in header, looping over configs, with",nmol,"molecules"
print "Opening files"

if(args.td):
    fo = []
    for i in range(1,nmol):
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
        if(args.td):
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

