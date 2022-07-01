#eom

##Equations of Motion

The Equations of Motion project provides an application with built in
astrodynamics related functionality along with a library to aid in the
creation of custom tools.

###eomx

The **eomx** application is being designed to focus on orbit
propagation, orbit determination, relative orbital dynamics, and
attitude dynamics.  It is a command line program driven by a simple text
based modeling language.  The library allows rapid development of custom
applications.

Current functionality includes the Vinti 2-body Kepler1 and J2/J3 Vinti6
orbit propagators.  The modern IAU GCRF and ITRF reference frames are
supported.  Parsing of IERS EOP data still needs to be implemented,
however the The IAU 2000A and IAU 2006 precession-nutation theories are
fully functional along with internal support for a true equator, mean
equinox (GMST 1980) ECI reference frame to support legacy general
perturbations propagators (such as the Vinti, traditional secular J2,
and SGP based propagators).  Propagated ephemeris can be saved in STK
compatible **.e** file formats.  Other outputs are written as Matlab
(Octave compatible) functions.  When run, these **.m** files will plot
the data with appropriate formatting while optionally returning the
handle to the figure along with the raw data.

**eom** utilizes two mature libraries upon which built-in models will
rely and external libraries may leverage.  The first is the
International Astronomical Union Standards of Fundamental Astronomy (IAU
SOFA) C library <http://www.iausofa.org/>.  Installation of SOFA-Issue
2015-02-09 or newer is required to build VMSAT.  The second library is
the Eigen 3 C++ template library for linear algebra
<http://eigen.tuxfamily.org>.  In the future, Kitware’s VTK OpenGL
library will be used for direct plotting of graphics.

Core concepts behind **eom** are:

1. Provide useful built in astrodynamics models and analysis tools
2. Provide convenient access to SOFA functionality
3. Supply an interface to external legacy libraries that have fallen out
   of use simply because the GUI supporting them can no longer be
   maintained
4. Provide a tool that allows for easy comparison of results from
   different functions and/or libraries
5. Provide a framework in which short response analysis can be addressed
   or algorithms prototyped

C++ was chosen for development because:

1. It plays well with legacy C and FORTRAN libraries
2. A number of mature and fast C++ math libraries are available
3. Code designed to run on hardware can be tested

The Mozilla Public License 2.0 (MPL2) was chosen because it protects the
open source nature of this code while still allowing it to be used with
proprietary and closed source tools:
<https://www.mozilla.org/en-US/MPL/2.0/>

