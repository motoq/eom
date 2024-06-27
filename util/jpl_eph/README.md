JPL Sun, Moon, and Planetary Ephemerides
========================================

The utilities here are used to generate eom specific binary (unformatted)
ephemeris files.  These files are optional.  The use of the precision
sun and lunar ephemerides compared to the built in Meeus algorithms for
orbit propagation, access analysis, etd., is negligible.  However, this
is currently the only source of planetary ephemerides for use with eom.

The eom compatible sun and moon files (sun.emb, moon.emb) are earth
centered inertial allowing for direct interpolation.  The planetary
files (e.g., pluto.emb) are sun centered.  Interpolation of both sun
position and planet position is required, with the earth relative
position of the planet being the sum.

The JPL asc2eph application is used to generate a SPICE compatible
binary (unformatted) ephemerides file from one or more text (formatted)
ascpYYYY.4xx files.  A modified version of the JPL testeph program,
gen_eom_eph, is run to compute individual sun, moon, and planetary eom
compatible .emb binary files.  These files must be accessible from the
directory in which eomx is run (adding a directory to search is on the
to-do list).

More information about the JPL planetary ephemerides can be found here:

<https://ssd.jpl.nasa.gov/planets/eph_export.html>

Jacob Williams describes obtaining and the JPLEPH ephemeris file:

<https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit#ephemeris-files>

Obtain code and ASCII text ephemerides:

wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/*
wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/*

On a PC with Linux, the software can be built via:

  #edit asc2eph.f file to set NRECL = 4
$ gfortran asc2eph.f -o asc2eph

  # Convert ephemerides for 2020-2040 to binary JPLEPH file
$ cat header.405 ascp2020.405 | ./asc2eph

To use JPL's testeph program on a similar architecture, update
multiple instances of NRECL to equal 4, as was done with asc2eph.f.  In
addition, set KSIZE to be equal to 2036 for use with enabling the
FSIZER3 option with DE404 or DE440 ephemerides.  The original JPL code
is well documented.  The code here has been modified to work on a PC:

  # Test converted ephemerides
Build testeph after updating multiple NRECL=4 and uncommenting FSIZER3
to use KSIZE = 2036

$ gfortran testeph.f -o testeph

Note IEEE underflow errors occuring when running testeph are considered
normal...

$ cat testpo.405 | ./testeph

  # Generate eom compatible files
To generate the eom compatible sun, moon, and planetary .emb files:

$ gfortran gen_eom_eph.f -o gen_eom_eph
$ ./gen_eom_eph

