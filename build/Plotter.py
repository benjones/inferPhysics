# -*- coding: utf-8 -*-
"""
Created on Mon Oct 03 09:39:43 2016

@author: Frank
"""

import matplotlib.pyplot as plt
import numpy as np

f = open('predictedPath.txt', 'rb')
f1 = open('actualPath.txt','rb')

xs = np.fromfile(f,float, -1, '')

# Added in actual points to see the difference between generated data and data using Mguess
# Will add in code to add in initial points from generated data so they don't get hard coded in and 
# be read from the file as well.
actual = np.fromfile(f1, float, -1,'')

# Black * is data generated from our Mguess, Red line is actual path.
plt.plot(xs, 'k*')
plt.plot(actual, 'r')

