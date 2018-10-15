# -*- coding: utf-8 -*-
"""Test Program for spktype01.


@author: Shushi Uetsuki (whiskie14142)
"""

from spktype21 import SPKType21
import csv
import numpy as np

bspfile = 'Ryugu-21.bsp'
csvfile = 'RyuguEPH.csv'
spkid = 2162173

dposlim = 1.0       # position difference limit = 1.0 kilometer
dvellim = 1e-7      # velocity difference limit = 0.0001 m/s

kernel = SPKType21.open(bspfile)
print()
print(kernel)
csvfile = open(csvfile)

print()
print('No.  JD          p_dist(km)              v_dist(km/s)\n')
count = 0
for row in csv.reader(csvfile):
    count += 1
    refjd = float(row[0])
    refpos = np.array([float(row[2]), float(row[3]), float(row[4])])
    refvel = np.array([float(row[5]), float(row[6]), float(row[7])])
    
    spkpos, spkvel = kernel.compute_type21(0, spkid, refjd)
    
    dpos = refpos - spkpos
    dvel = refvel - spkvel
    p_dist = np.sqrt(np.dot(dpos, dpos))
    v_dist = np.sqrt(np.dot(dvel, dvel))
    print(count, ' ', refjd, ' ', p_dist, ' ', v_dist)
csvfile.close()
kernel.close()
