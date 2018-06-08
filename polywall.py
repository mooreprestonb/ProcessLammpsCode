#!/usr/bin/python

import math

numa = 5
numb = 10

# 1/2 unit cells
xdis = .5
ydis = .5*math.sqrt(3)
p1 = [0,0,0]
p2 = [xdis/2,ydis/2,0]

xunits = 100
yunits = 50
zunits = 50 
box  = [xunits,yunits*math.sqrt(3),zunits]

# make it go from -L/2 to L/2
xunits = 2*int(xunits/2)
yunits = 2*int(yunits/2)
natomswall = 2*xunits*yunits

natoms = natomswall+numa+numb

print "Lammps init wall file ",xunits,yunits,natomswall,numa,numb
print " "
print natoms, "atoms"
print "1 atom types"
print " "
print -box[0]/2,box[0]/2,"xlo xhi"
print -box[1]/2,box[1]/2,"ylo yhi"
print 0,box[2],"zlo zhi"
print " "
print "Masses"
print " "
print "1 1.0"
print " " 
print "Atoms # Full"
print " " 

ip = 0 # ith particle
gp = 1 # group id
tp = 1 # type
ch = 0 # charge
for i in range(xunits):
    for j in range(yunits):
        ip += 1
        xp = i*xdis+p1[0]-box[0]/4
        yp = j*ydis+p1[1]-box[1]/4
        zp = p1[2]
        print ip,gp,tp,ch,xp,yp,zp
        
        ip += 1
        xp = i*xdis+p2[0]-box[0]/4
        yp = j*ydis+p2[1]-box[1]/4
        zp = p2[2]
        print ip,gp,tp,ch,xp,yp,zp

# done with wall

inuma = ip+1
gp = 2
tp = 2
ch = 0
zbond = 1.5
zstart = zbond
for i in range(numa):
    ip += 1
    xp = 0
    yp = 0
    zp = i*zbond+zstart
    print ip,gp,tp,ch,xp,yp,zp

gp = 2
tp = 3
ch = 0
zbond = 1.5
zstart = zbond*(numa+1)
for i in range(numb):
    ip += 1
    xp = 0
    yp = 0
    zp = i*zbond+zstart
    print ip,gp,tp,ch,xp,yp,zp

#bonds

nbonds = numa+numb

print ""
print "Bonds"
print ""

polybond = xunits*yunits+xunits/2+1
print 1, 2, polybond, inuma
for i in range(1,nbonds):  # i type a b
    print i+1, 1, i+inuma-1, i+inuma
    
