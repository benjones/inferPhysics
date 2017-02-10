# -*- coding: utf-8 -*-
"""
Created on Mon Oct 03 09:39:43 2016

@author: Frank
"""

import matplotlib.pyplot as plt
import numpy as np
import array

f = open('predictedPath', 'rb')
f1 = open('actualPath','rb')

reconstructedPath = np.fromfile(f,float, -1, '')

# Added in actual points to see the difference between generated data and data using Mguess
# Will add in code to add in initial points from generated data so they don't get hard coded in and 
# be read from the file as well.
snapshots = np.fromfile(f1, float, -1,'')

ind = array.array('f')
for i in range(0, len(reconstructedPath)):
    ind.append(i/10.0)

index = np.array(ind)

#first index is for projectile motion, second is spring force.
snapInd = np.array([0,2,4,6,8,10,12,14])
#snapInd = np.array([0,2,4,5,7,10,12,14])
#snapInd = np.array([0,2,6,12,19])

# Black * is data generated from our Mguess, Red line is actual path.
plt.plot(snapInd, snapshots, 'k*')
plt.plot(index, reconstructedPath, 'r.')
plt.axis([-1, 16, -10, 75])
#plt.axis([-2, 11, -2, 75])
