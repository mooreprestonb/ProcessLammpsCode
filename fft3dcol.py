#!/usr/bin/python

import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
import math
from scipy.fftpack import dct, idct # direct and inverst cosine transform

#setup arguments to use in argparse (argument parsing)

parser = argparse.ArgumentParser(description='Fourier Transform of column data')
parser.add_argument('input', help="input file")
parser.add_argument('output', nargs='?', help="output file")
parser.add_argument('-n',type=int,default=-1,help="column to transform")
parser.add_argument("-ct",help="Perform cosine transform",action="store_true")
parser.add_argument("-w",help="Scale from 1/fs to 1/cm",action="store_true")
parser.add_argument("-gr",help="Graph results",action="store_true")
args = parser.parse_args()
print args

data = np.loadtxt(args.input)

num = data.shape[0]
x = data[:,0]

print "Read in ",sys.argv[1]," with ",num
hdr = " ".join(sys.argv) # create header from sys.argv

xt = np.zeros(num)
xf = np.zeros(num)
dt = x[1]-x[0]
tol = 1e-8
fact = 1.0
if (args.w):
    fact = 1./2.99792458e-5

for i in range(num):
    xt[i] = i*dt
    xf[i] = fact*i/(dt*num)
    if (i<(num-1)):
        if (((x[i]-x[i+1])-dt) > tol) :
            print "Error: Steps not equil at step i=",i,x[i]-x[i+1],dt

print "Doing rfft", xt.shape,"Column:",args.n
if (num%2):
    oa = num/2+1
else:
    oa = num/2

outa = np.array(xf[:oa],ndmin=2)  # 
if(args.n == -1):
    for i in range(1,data.shape[1]):
        y = data[:,i]
        ftcorr = np.fft.rfft(y)
        tarray = np.array([np.real(ftcorr),np.imag(ftcorr)])
        #print outa.shape, ftcorr.shape, tarray.shape
        outa = np.append(outa,tarray,axis=0)
else:
    y = data[:,int(args.n)]
    ftcorr = np.fft.rfft(y)
    tarray = np.array([np.real(ftcorr),np.imag(ftcorr)])
    outa = np.append(outa,tarray,axis=0)  # 

if(args.output==None):
    output = args.input + ".fft"
else:
    output = args.output

np.savetxt(output,outa.T,fmt="%f",header=hdr + " fft") # output with x y format

if(args.ct):
    dctcorr = dct(y,type=1)
    output += ".cos"
    outa = np.column_stack((xf/2.0,dctcorr))  # 
    np.savetxt(output,outa,fmt="%f",header=(hdr + " cos")) # output with x y format
    
print "Done, creating graphs"

if(args.gr): # use graphs
    plt.plot(xt,y,label="Correlation")
    plt.legend()
    plt.show()
    plt.plot(xf[:ftcorr.shape[0]],np.real(ftcorr),label="Real")
    plt.plot(xf[:ftcorr.shape[0]],np.imag(ftcorr),label="Imag")
    if(args.ct):
        plt.plot(xf/2.0,dctcorr/2.0,label="Cos") # cosine has twice as many points
    plt.legend()
    plt.show()

