             
	Title:   INTERPOLATING IERS EARTH ORIENTATION DATA
	Authors: Dennis McCarthy, Daniel Gambis, Ch. Bizouard, 
	
	Revised version of the IERS Gazette No 13, 29 January 1997


INTRODUCTION

     IERS Earth Orientation Parameters (EOP) are produced at regular intervals
(daily or longer) with the effects of semidiurnal and diurnal variations 
removed. Users who require high accuracy information may want to interpolate
the published data and include the semidiurnal/diurnal variations.
 
     This gazette provides a recommended procedure to follow in order to
determine the most accurate Earth orientation for a given instant.

     Currently all analysis centers contributing to the IERS employ a 
procedure using a priori values of the EOP which are then corrected through 
the analysis of observations.  These a priori values are best estimates based 
on the standard knowledge of the Earth's orientation, usually the interpolated 
tabular values plus a diurnal/semidiurnal model plus lower frequency tides 
(if appropriate). The corrections to the a priori estimates are determined 
from the analyses of data taken over some period of time, ranging from minutes 
to days. Thus, they represent the mean of the unmodelled variations over the 
length of the analyzed period of time. When the reported estimates (a priori 
+ mean estimated correction) are based on distinct time intervals, they may 
be referred to as normal points.
 
     The IERS provides to users polar motion and UT1 based on the combination 
of the analysis centers data, using either smoothed estimates in order to 
reduce observational errors (Bulletins A and B, EOP 97 C 04, see Explanatory 
Supplement to IERS Bulletins A and B), or normal points (EOP 97 C 01, 02, 03, 
see IERS 1995 Annual Report, part II.4). These results are provided without the 
effects of the diurnal/semidiurnal tides. The IERS Conventions (McCarthy 1996) 
for transformations between terrestrial and celestial frames, however, imply 
that IERS EOP data are exact estimates at the instant in time reported.

SEMIDIURNAL/DIURNAL VARIATIONS CAUSED BY OCEANIC TIDES

     The existence of diurnal and semidiurnal variations caused by ocean tides 
is well known (Eubanks 1993; Sovers et al. 1993; Herring & Dong 1994; Ray et 
al. 1994).  In recent years, through continuous, high-precision VLBI 
experiments (e.g. ERDE, EPOCH92, CONT93, etc.), analytical models have been 
derived for the prominent diurnal and semidiurnal tides (Herring and Dong 1991;
Brosche et al. 1991; Herring 1993; Herring and Dong 1994; Gipson 1996).  
This in turn has prompted a refinement in the theoretical models (Brosche et 
al. 1989; Seiler 1991; Dickman 1989, 1990, 1991, 1993; Gross 1993; Ray et al. 
1994; Ray 1995; Seiler and Wunsch 1995).

     As the observational data have improved in accuracy, the models (both 
empirical and theoretical) have quickly converged. Although there are 
differences between recent models (e.g. Ray et al. 1994; Ray 1995; Seiler 
and Wunsch 1995; Gipson 1996), these differences are less than 10%. The 
effects of these tides are on the order of 0.1 milliarcseconds (mas) in the 
polar motion coordinates x and y and 10**-6 s in UT1.  The inclusion of any 
of these models in high-precision, Earth orientation analysis software clearly 
improves the solution.

     If the data do not correctly model the diurnal/semidiurnal tide, the 
errors can reach up to +/-0.1 milliarcseconds (mas). Systematic errors of this 
magnitude may be unacceptable to high-accuracy users.  Possible ambiguities 
caused by unclear procedures concerning subdiurnal EOPs and the epoch of 
observation can be eliminated.  Although the problems are only on the fringes 
of detectability, they may cause significant systematic errors. It is expected 
that eventually subdiurnal observations of EOPs will be possible routinely from 
a reduction standpoint.

LUNISOLAR DIURNAL VARIATIONS CAUSED 

Another cause of diurnal variation in PM is the lunisolar effect associated with high degree of the geopotential. The IERS 2000 recommended table takes into account two models for non-rigid Earth, one by Mathews and Bretagnon (2002, Proc. JSR 2001) and the other one by Brzezinski (2001, Proc. JSR 2000), Brzezinski and Capitaine (2002, Proc. JSR 2001). The papers describing these theories are in preparation and should be soon submitted for publication in the reviewed international journals. After discussing several details and introducing all the necessary corrections, Mathews and Brzezinski could reach a consensus on the attached table. 

RECOMMENDATION

     The following software is recommended to interpolate the IERS polar motion 
and Universal Time products and account for the semidiurnal/diurnal variations 
in the Earths orientation. This procedure makes use of a Lagrangian 
interpolation scheme and applies the integral Ray model (71 tidal waves) and Brzezinski-Mathews-Bretagnon-Capitaine-Bizouard model (10 lunisolar waves) of the semidiurnal/diurnal variations in the Earth's orientation as recommended in the IERS 2000 Conventions
(McCarthy, 2002). Any equivalent interpolation scheme could, of course, be 
substituted.  

This software can be obtained in machine readable form on : 

http://hpiers.iers.obspm.fr/eop-pc/models/PM/interp.f 

REFERENCES

Brosche, P., Seiler, U., Sundermann, J., and Wunsch, J., 1989,
     "Periodic Changes in Earth's Rotation due to Oceanic Tides,"
     Astron. Astrophys., 220, pp. 318-320.

Brosche, P., Wunsch, J., Campbell, J., and Schuh, H., 1991,
     "Ocean Tide Effects in Universal Time detected by VLBI,"
     Astron. Astrophys., 245, pp. 676-682.

Dickman, S. R., 1989, "A Complete Spherical Harmonic Approach to
     Luni-Solar Tides," Geophys. J., 99, pp. 457-468.

Dickman, S. R., 1990, "Experiments in Tidal Mass Conservation,"
     (research note), Geophys. J., 102, pp. 257-262.

Dickman, S. R., 1991, "Ocean Tides for Satellite Geodesy," Mar.
     Geod., 14, pp. 21-56.

Dickman, S. R., 1993, "Dynamic ocean-tide effect on Earth's
     rotation," Geophys. J. Int., 112, pp. 448-470.

Eubanks, T. M., 1993, "Variations in the Orientation of the
     Earth," in Contributions of Space Geodesy to Geodynamics:
     Earth Dynamics Geodynamics, American Geophysical Union,
     Washington, DC, 24, pp. 1-54.

Gipson, J., 1996, "VLBI Determination of Neglected Tidal Terms in
     High-Frequency Earth Orientation Variation," submitted to J.
     Geophys. Res.

Gross, R. S., 1993, "The Effect of Ocean Tides on he Earth's
     Rotation as Predicted by the Results of an Ocean Tide
     Model," Geophys. Res. Lett., 20, pp. 293-296.

Herring, T. A., 1993, "Diurnal and semi-diurnal variations in
     Earth rotation," in The Orientation of the Planet Earth as
     Observed by Modern Space Techniques, M. Feissel (ed.),
     Pergamon Press, in press.

Herring, T. A. and Dong, D., 1991, "Current and future accuracy
     of Earth rotation measurements," in Proceedings of the
     Chapman conference on Geodetic VLBI: Monitoring Global
     Change, NOAA Technical Report NOS 137 NGS 49, pp. 306-324.

Herring, T. A. and Dong, D., 1994, "Measurement of Diurnal and
     Semidiurnal Rotational Variations and Tidal Parameters of
     Earth," J. Geophys. Res., 99, pp. 18051-18071.

IERS Annual Report, Available from the Central Bureau of the
     IERS, Paris Observatory, Paris.

McCarthy, D. D., 1996, IERS Conventions, IERS Technical Note, 21.

McCarthy, D. D., 2000, IERS Conventions 2000, in preparation, see http://maia.usno.navy.mil/conv2000.html

Ray, R., Steinberg, D. J., Chao, B. F., and Cartwright, D. E.,
     1994, " Diurnal and Semidiurnal Variations in the Earth's
     Rotation Rate Induced by Oceanic Tides," Science, 264, pp.
     830-832.

Ray, R., 1995, Personal Communication.

Seiler, U., 1991, "Periodic Changes of the Angular Momentum
     Budget due to the Tides of the World Ocean," J. Geophys.
     Res., 96, pp. 10287-10300.

Seiler, U. and Wunsch, J. 1995, "A refined model for the
     influence of ocean tides on UT1 and polar motion," Astron.
     Nachr., 316, pp. 419-423.

Sovers, O. J., Jacobs, C. S., and Gross, R. S., 1993, "Measuring
     rapid ocean tidal Earth orientation variations with VLBI, J.
     Geophys. Res., 98, 19959-19971.

