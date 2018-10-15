# -*- coding: utf-8 -*-
"""Multi-object SPK file test for spktype21.
CPJV.bsp contains spk data for four asteroids.

@author: Shushi Uetsuki (whiskie14142)
"""

from spktype21 import SPKType21
import numpy as np

oname = ['', 'Ceres', 'Pallas', 'Juno','Vesta']
spkid = [0, 2000001, 2000002, 2000003, 2000004]

path21 = 'CPJV.bsp'
k21 = SPKType21.open(path21)

print()
print(k21)
print()

while(True):
    asteroid = int(input('Enter asteroid No.(1 - 4) or zero for exit ? '))
    print()
    if (asteroid < 1) or (asteroid > 4): break
    while(True):
        jd = float(input('Enter JD or 999 for break ? '))
        if jd <1000.0 : break
        pos21, vel21 = k21.compute_type21(0, spkid[asteroid], jd)
        print('POS of ' + oname[asteroid] + ': ', pos21)
        print('VEL of ' + oname[asteroid] + ': ', vel21)
        print()

k21.close()
