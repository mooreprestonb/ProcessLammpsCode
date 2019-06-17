#!/usr/bin/python
# read lammps data and reorder...

import sys
import re
import argparse
import numpy as np

def indexline(al,value):
    n = 0
    for x in al:
        n += 1
        #print x.rstrip()
        w = x.rstrip().split()
        if w:  # not empty )
            if w[0] == value: 
                break
    if n == len(al):
        print "Warning, can't fine", value
#        exit(1)    
    return n
#----------------------------

parser = argparse.ArgumentParser(description="Reorder lammps init files")
parser.add_argument("infile",help="Lammps file name to process")
parser.add_argument("outfile",help="Lammps output file to create")
args = parser.parse_args()

#lammpsfile = "lammps.test.init"
lammpsfile = args.infile
ofile = args.outfile

print "Reading in Lammps file",lammpsfile

f = open(lammpsfile,"r")
lines = list(f)
f.close()

f = open(ofile,"w")

l=0
#Find natoms
for x in lines:
    l += 1
    #print x.rstrip()
    w = x.rstrip().split()
    if len(w) > 1:  # list has at least 2 elements (and not empty ;-) )
        if w[1] == "atoms": # found keyword "atoms"
            natoms = int(w[0])
            if (l != 3):
                print "Atoms not on line 3! ",l
                exit(1)
            break
print "natoms in header = ",natoms

#Find atoms types
for x in lines:
    l += 1
    #print x.rstrip()
    if(re.search("atom types",x)):
        w = x.rstrip().split()
        atmtypes = int(w[0])
        break
print "Atoms types = ",atmtypes

# find box dimensions
l=0
lbox = -1
box = []
for x in lines:
    l += 1
    #print x.rstrip()
    w = x.rstrip().split()
    if (lbox > -1):  # box line
        lbox += 1
        box.extend([float(w[0]),float(w[1])])
        if (lbox > 2): # read in three lines!
            break
    if (not w and l>3):  # first empty line after line 3
        lbox += 1
print "box read in", box

#find masses
nmass = -1
mass = []
for x in lines:
    if (nmass > -1): # reading in masses
        w = x.rstrip().split()
        if w: # not a blank line
            # print w,natm
            if not w[0].isdigit():
                # print w[0],"not digit!"
                break
            nmass += 1
            mass.append([int(w[0]),float(w[1])])
        
    if(re.search("Masses",x)):
        print "found Masses"
        nmass = 0
    
print "Masses found = ",nmass
#print mass, mass[0][1]
if (nmass != atmtypes):  # error check
    print "Something wrong, atom types do not match", atmtypes, " masses found !=", nmass
    exit()

#read natoms
natm = -1
atms = []
indx = []
gindx = []
hl = 0
nl = 0
for x in lines:
    w = x.rstrip().split()
    nl += 1
    if w: # not a blank line
        # print w,natm
        if (natm > -1):
            if not w[0].isdigit():
                # print w[0],"not digit!"
                break
            natm += 1
            atms.append(w) # atom line
            indx.append(int(w[0]))  # atom index
            gindx.append(int(w[1])) # group index
        if w[0] == "Atoms": # start counting Atoms!
            natm = 0
            hl = nl # set header line count!
            
print "natoms found = ",natm

if (natm != natoms):  # error check
    print "Something wrong, natoms do not match", natm, " found !=", natoms
    exit()

#print atms

indxs = np.argsort(indx) # index to sort array
gindxs = np.argsort(gindx) # group index to sort array
sgindx = sorted(set(gindx)) # sorted group index

#print indx, indxs
#print gindx, sgindx

#print header
for i in range(hl):
    f.write(lines[i])
f.write("\n")
print "Wrote header"

# print renumbered atoms and groups
ni = 0
aridx = []  # get indexed array
for i in indxs:
    aridx.append(i);
    words = atms[i];
    l = str(aridx.index(i)+1) + " "
    igrp = sgindx.index(int(atms[i][1]))+1
    l += str(igrp) + " "
    for w in words[2:]:
        l += w + " "
    f.write(l)
    f.write("\n")
    ni += 1
f.write("\n")
    
print "Wrote atoms with new index and group id's"

# Now for other data
nl = indexline(lines,"Atoms")
vl = indexline(lines,"Velocities")
bl = indexline(lines,"Bonds")
al = indexline(lines,"Angles")
dl = indexline(lines,"Dihedrals")
il = indexline(lines,"Impropers")
#dl = len(lines)+1

#print "Velocities"
#print
#for i in range(vl+1,bl-1):
#    w = lines[i].rstrip().split()
#    # print w
#    if w: # skip over blank lines
#        ni = int(w[0])
#        nn = aridx.index(indx.index(ni))+1
#        l = str(nn) + " "
#        for wds in w[1:]:
#            l += wds + " "
#        print l 
#print

f.write("Bonds\n\n")
for i in range(bl+1,al-1):
    w = lines[i].rstrip().split()
    # print w
    if w: # skip over blank lines
        l = w[0] + " " + w[1] + " "
        ni = int(w[2])
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + " "
        ni = int(w[3])
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + "\n"
        f.write(l) 
f.write("\n")
print "Wrote Bonds"

#Find nangles 
for i in range(hl):
    #print x.rstrip()
    if(re.search("angles",lines[i])):
        w = lines[i].rstrip().split()
        nangles = int(w[0])
        break

print "Nangles = ",nangles
if not (nangles == (dl-al-3)):
    print "Warning nangles do not match lines to dihedrals",nangles,dl-al-3

f.write("Angles\n\n")
for i in range(nangles):
    w = lines[i+al+1].rstrip().split()
    if w: # skip over blank lines
        l = w[0] + " " + w[1] + " "
        ni = int(w[2])
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + " "

        ni = int(w[3]) 
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + " "
        
        ni = int(w[4])
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + "\n"
        
        f.write(l) 
f.write("\n")
print "Wrote Angles",nangles

#Find ndihedrals
ndihedrals = 0
for i in range(hl):
    #print x.rstrip()
    if(re.search("dihedrals",lines[i])):
        w = lines[i].rstrip().split()
        ndihedrals = int(w[0])
        break

print "ndihedrals = ",ndihedrals
if not (ndihedrals == (il-dl-3)):
    print "Warning nangles do not match lines to impropers",ndihedrals,il-dl-3

f.write("Dihedrals\n\n")
for i in range(ndihedrals):
    w = lines[i+dl+1].rstrip().split()
    # print w
    if w: # skip over blank lines
        l = w[0] + " " + w[1] + " "  # dihedral # and type

        ni = int(w[2])
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + " "

        ni = int(w[3]) 
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + " "
        
        ni = int(w[4])
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + " "

        ni = int(w[5])
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + "\n"

        f.write(l)
f.write("\n")
print "Wrote Dihedrals",ndihedrals

#Find impropers
nimpropers = 0
for i in range(hl):
    #print x.rstrip()
    if(re.search("impropers",lines[i])):
        w = lines[i].rstrip().split()
        nimpropers = int(w[0])
        break

print "nimpropers = ",nimpropers
if not (nimpropers == (il-len(lines))):
    print "Warning nangles do not match lines to impropers",nimpropers,len(lines)


if(nimpropers):
    f.write("Impropers\n\n")
for i in range(nimpropers):
    w = lines[i+nimpropers+1].rstrip().split()
    # print w
    if w: # skip over blank lines
        l = w[0] + " " + w[1] + " "  # dihedral # and type

        ni = int(w[2])
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + " "

        ni = int(w[3]) 
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + " "
        
        ni = int(w[4])
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + " "

        ni = int(w[5])
        n1 = aridx.index(indx.index(ni))+1
        l += str(n1) + "\n"

        f.write(l)
if(nimpropers):
    f.write("\n")

print "Wrote Impropers",nimpropers
