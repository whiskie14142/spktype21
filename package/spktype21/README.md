# spktype21
A supporting module for [jplephem](https://pypi.org/project/jplephem/) to handle data type 21 (Version 0.1.0)

This module computes positions and velocities of a celestial small body, from a NASA SPICE SPK ephemeris kernel file of data type 21 (Extended Modified Difference Arrays).  See http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

You can get SPK files for many solar system small bodies from HORIZONS system of NASA/JPL.  See https://ssd.jpl.nasa.gov/?horizons

This module reads SPK files of data type 21, one of the types of binary SPK file. At the point of Oct. 2018, HORIZONS system provides files of type 21 as binary SPK files by default.  You can get type 21 binary SPK file for celestial small bodies through TELNET interface by answering back 'B' for 'SPK file format'. Also you can get type 21 binary SPK file from:
https://ssd.jpl.nasa.gov/x/spk.html

### Modules required
* jplephem (version 2.6 or later)
* numpy

### Usage
    from spktype21 import SPKType21
    kernel = SPKType21.open('path')
    position, velocity = kernel.compute_type21(center, target, jd)
    print(kernel)     ---- this line prints information of all segments
    kernel.close()
    
    where:
        path - path to the SPK file
        center - SPKID of central body (0 for SSB, 10 for Sun, etc.)
        target - SPKID of target body
        jd - time for computation (Julian date)
        position - a numpy array (x, y, z)
        velocity - a numpy array (xd, yd, zd)

### Modification Log
##### 0.1.0 October 15, 2018
* Beta Release
