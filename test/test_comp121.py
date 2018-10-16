# -*- coding: utf-8 -*-
"""Compare results of spktype01 and spktype21.
@author: Shushi Uetsuki (whiskie14142)
"""

from spktype01 import SPKType01
from spktype21 import SPKType21
import numpy as np

oname = 'Ryugu'
spkid = 2162173

path01 = oname + '-01.bsp'
path21 = oname + '-21.bsp'

k01 = SPKType01.open(path01)
k21 = SPKType21.open(path21)

print('\nspktype01')
print(k01)
print('\nspktype21')
print(k21)
print()

while(True):
    jd = float(input('Enter JD, or 999 to exit >> '))
    if jd == 999.0:
        break
    else:
        print('\nspktype01')
        pos01, vel01 = k01.compute_type01(0, spkid, jd)
        print(pos01)
        print(vel01)
        print('\nspktype21')
        pos21, vel21 = k21.compute_type21(0, spkid, jd)
        print(pos21)
        print(vel21)
        print('\ndistances')
        dpos = pos01 - pos21
        dvel = vel01 - vel21
        print(np.sqrt(np.dot(dpos, dpos)), np.sqrt(np.dot(dvel, dvel)))
        print()

k01.close()
k21.close()
