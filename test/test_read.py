# -*- coding: utf-8 -*-
"""Read SPK file by spktype21.
The SPK file should be an 'One object' file.

@author: Shushi Uetsuki (whiskie14142)
"""

from spktype21 import SPKType21
import numpy as np

print()
path21 = input('Enter SPK file name (without extention) >> ') + '.bsp'
k21 = SPKType21.open(path21)
spkid = k21.segments[0].target

print()
print(k21)
print()

while(True):
    jd = float(input('Enter JD, or 999 to exit >> '))
    if jd <1000.0 : break
    pos21, vel21 = k21.compute_type21(0, spkid, jd)
    print('POS: ', pos21)
    print('VEL: ', vel21)
    print()

k21.close()
