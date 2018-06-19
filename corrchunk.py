#!/usr/bin/python

import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
import math
from scipy.fftpack import dct, idct # direct and inverse cosine transform

def output (fo, hdr, ts, time, x, sd, md):
    mnum = x.shape[0]
    print ts,"/",time, mnum

    sdc = np.array(sd) # total system dipole corr
    mdc = np.array(md) # molecular dipole corr

    if (ts < time-mnum):
        sdc /= time
        mdc /= time
    else:
        for i in range(mnum): # uneven averaging
            sdc[i] = sdc[i]/((time-i)
            if(args.td):
                mdc[i] = mdc[i]/(time-i)

    outa = np.column_stack((x,sdc,mdc))  # need to pass one argument to get x y
    np.savetxt(args.output+str(i),outa,fmt="%f",header=hdr) # output with x y z format

                
#setup arguments to use in argparse (argument parsing)

parser = argparse.ArgumentParser(description='Autocorrelate dipole of lammps chunk output (time mux muy muz tdp')
parser.add_argument('input', help="input file")
parser.add_argument('output', help="output file")
parser.add_argument('-n',type=int,default=-1,help="#'s to correllate")
parser.add_argument('-ts',type=float,default=-1,help="time step i.e dt")
parser.add_argument('-rp',type=int,default=10,help="Report progress every \"rp\" steps")
parser.add_argument("-av",help="do NOT subtract average",action="store_true")
parser.add_argument("-s",help="Scale C[0] to 1",action="store_true")
parser.add_argument("-gr",help="Graph resutls",action="store_true")
parser.add_argument("-td",help="Do total dipole only",action="store_false")
args = parser.parse_args()

print args
#data = np.loadtxt(args.input) # read in special
with open(args.input) as f:
    lines = f.readlines()

print len(lines),"lines read in, parsing"

xt = []
data = []
dpt = []
ts = []
nmol = -1
dt = -1
time = 0
tol = 1e-8

for l in lines:
    w = l.split()
    if (not w[0] == '#'): # skip comments (first 3 lines)
        if (len(w) == 2):
            time += 1
            ts.append(float(w[0]))
            if (len(ts)==2):
                dt = ts[1]-ts[0]

            if (dt != -1): # check step
                if (abs(dt-(ts[-1]-ts[-2]))>tol):
                    print "Time step not equil?",time,dt,ts[-2],ts[-1]
                
            if (nmol == -1):
                nmol = int(w[1])
            
            if (nmol != int(w[1])):
                print "Error nmol has changed"
                print nmol, w[1]
                exit()
            
        else:
            xt.append(float(w[0]))
            data.append([float(w[1]),float(w[2]),float(w[3])])
            # dpt.append(float(w[4])) # don't need total dipole

x = np.array(xt)
dy = np.array(data)
tdp = np.array(dpt)

print time,nmol,3,dy.shape
dy = np.reshape(dy,(time,nmol,3)) # create array with ntime x nmol x 3

if (args.n == -1): # set numbers to correlate
    mnum = time
else:
    mnum = min(args.n,time)
    
corr = np.zeros(mnum)
corrtd = np.zeros(mnum)
nd = 3 # 3 dimensions

print time, nd, mnum, args.rp

avg = np.zeros((nmol,nd))

for j in range(nmol):
    for i in range(nd):
        avg[j][i] = np.average(dy[:,j,i])

print "Averages: ",avg
print "Variance: ",np.var(dy)

if(not args.av):
    print "Subtracting off average"
    dy = dy - avg # subtract off average

print "Calculating total dipole for each step"
tdp = np.empty((time,3))
for i in range(time):
    tdp[i] = np.sum(dy[i],axis=0) # total dipole

outp = max(time/args.rp,1)
hdr = " ".join(sys.argv) # create header from sys.argv
hdr += "\n time TotalDipole MolDipole" # create header from sys.argv

outa = np.column_stack((np.array(ts),tdp))  # need to pass one argument to get x y
np.savetxt(args.input+".td",outa,fmt="%f",header=hdr) # output with x y z format

xt = np.arange(mnum)*dt # create xt [0, dt, 2dt, 3dt, .... (mnum-1)dt]

print "Correlating: Reporting progress every",outp,"steps, with a window of",mnum,"steps and a total time of ",time

for i in range(time):
    if (not (i%outp)):
        #print i
        output (args.output, hdr, i, time, xt, corrtd, corr)

    # slice out what we want i to i+mnum (time window)
    for j in range(nmol):
#        print i,j,min(time,i+mnum)
        if(args.td):  # correlate individual molecule dipoles
            y = dy[i:min(time,i+mnum),j,:]
            y0 = y[0] # slice out step i
            corrd = np.dot(y,y0)  # correlate step i using numpy dot product
            corr[:corrd.shape[0]] += corrd # add to corr with correct shape

    # total dipole section
    y = tdp[i:min(time,i+mnum),:]
    y0 = y[0] # slice out step i
    corrd = np.dot(y,y0)  # correlate step i using numpy dot product
    corrtd[:corrd.shape[0]] += corrd # add to corr with correct shape
    
#print "Corr:",corr

print "Done Correlating: writing out"
    
for i in range(mnum):
    corrtd[i] = corrtd[i]/(time-i)
    if(args.td):
        corr[i] = corr[i]/(time-i)

# print corr,corrtd
print ", Now C[0] = ",corr[0],corrtd[0]

if(args.s): # scale to 1
    print "Scaling corr[0] to 1, C[0] = ",corr[0]
    corrtd = corrtd/corrtd[0]
    if(args.td):
        corr = corr/corr[0]

outa = np.column_stack((xt,corrtd,corr))  # need to pass one argument to get x y
np.savetxt(args.output,outa,fmt="%f",header=hdr) # output with x y z format
    
if(args.gr): # use graphs
    print "Creating graphs"
    plt.plot(xt,corr,label="Correlated dipole")
    plt.plot(xt,corrtd,label="Correlated total dipole")
    plt.legend()
    plt.show()

print "Done"

