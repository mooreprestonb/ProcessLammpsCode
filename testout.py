#!/usr/bin/python
import math
# create test data to correllate: sum of three cosines

W1 = 10
W2 = W1 * math.sqrt(2.0)
W3 = 3 * W1
A1 = 1
A2 = 2*A1
A3 = .5*A1
dt = .5

for i in range(2000):
    x = i*dt
    print x, A1*math.cos(2.0*math.pi*x/W1), A1*math.cos(2.0*math.pi*x/W1) + A2*math.cos(2.0*math.pi*x/W2) + A3*math.cos(2*math.pi*x/W3)
