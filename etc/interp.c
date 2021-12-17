/* interp.f -- translated by f2c (version 20190311).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b4 = 1.296e6;

/* ---------------------------------------------------------------- */
/* Subroutine */ int interp_(doublereal *rjd, doublereal *x, doublereal *y, 
	doublereal *ut1, integer *n, doublereal *rjd_int__, doublereal *
	x_int__, doublereal *y_int__, doublereal *ut1_int__)
{
    extern /* Subroutine */ int pm_gravi__(doublereal *, doublereal *, 
	    doublereal *), pmut1_oceans__(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal cor_x__, cor_y__;
    extern /* Subroutine */ int lagint_(doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *);
    static doublereal cor_ut1__, cor_lod__;


/*     This subroutine takes a series of x, y, and UT1-UTC values */
/*     and interpolates them to an epoch of choice. This routine */
/*     assumes that the values of x and y are in seconds of */
/*     arc and that UT1-UTC is in seconds of time. At least */
/*     one point before and one point after the epoch of the */
/*     interpolation point are necessary in order for the */
/*     interpolation scheme to work. */

/*     parameters are : */
/*     RJD     - array of the epochs of data (given in mjd) */
/*     X       - array of x polar motion (arcsec) */
/*     Y       - array of y polar motion (arcsec) */
/*     UT1     - array of UT1-UTC (sec) */
/*     n       - number of points in arrays */
/*     rjd_int - epoch for the interpolated value */
/*     x_int   - interpolated value of x */
/*     y_int   - interpolated value of y */
/*     ut1_int - interpolated value of ut1-utc */

/*     CALLED SUBROUTINE : LAGINT (Lagrange interpolation) */
/*                         PMUT1_OCEANS (Diurnal and semidiurnal oceanic effects) */
/*                         PM_GRAVI (Diurnal and semidiurnal lunisolar effects) */

/*      coded by Ch. BIZOUARD (Observatoire de Paris) : November 2002 */
/*                                          Corrected : September 2007 */
    /* Parameter adjustments */
    --ut1;
    --y;
    --x;
    --rjd;

    /* Function Body */
    lagint_(&rjd[1], &x[1], n, rjd_int__, x_int__);
    lagint_(&rjd[1], &y[1], n, rjd_int__, y_int__);
    lagint_(&rjd[1], &ut1[1], n, rjd_int__, ut1_int__);
/* -------------- */
/* Oceanic effect */
/* -------------- */
    pmut1_oceans__(rjd_int__, &cor_x__, &cor_y__, &cor_ut1__, &cor_lod__);
    *x_int__ += cor_x__;
    *y_int__ += cor_y__;
    *ut1_int__ += cor_ut1__;
/* Lunisolar effect */
    pm_gravi__(rjd_int__, &cor_x__, &cor_y__);
    *x_int__ += cor_x__;
    *y_int__ += cor_y__;
    return 0;
} /* interp_ */


/* ---------------------------------------------------------------- */

/* Subroutine */ int lagint_(doublereal *x, doublereal *y, integer *n, 
	doublereal *xint, doublereal *yout)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal term;


/*     This subroutine performs lagrangian interpolation */
/*     within a set of (X,Y) pairs to give the y */
/*     value corresponding to xint. This program uses a */
/*     window of 4 data points to perform the interpolation. */
/*     if the window size needs to be changed, this can be */
/*     done by changing the indices in the do loops for */
/*     variables m and j. */

/*     PARAMETERS ARE : */
/*     X     - array of values of the independent variable */
/*     Y     - array of function values corresponding to x */
/*     n     - number of points */
/*     xint  - the x-value for which estimate of y is desired */
/*     yout  - the y value returned to caller */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    *yout = 0.;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*xint >= x[i__] && *xint < x[i__ + 1]) {
	    k = i__;
	}
    }
    if (k < 2) {
	k = 2;
    }
    if (k > *n - 2) {
	k = *n - 2;
    }
    i__1 = k + 2;
    for (m = k - 1; m <= i__1; ++m) {
	term = y[m];
	i__2 = k + 2;
	for (j = k - 1; j <= i__2; ++j) {
	    if (m != j) {
		term = term * (*xint - x[j]) / (x[m] - x[j]);
	    }
	}
	*yout += term;
    }
    return 0;
} /* lagint_ */

/* ---------------------------------------------------------------- */
/* Subroutine */ int pmut1_oceans__(doublereal *rjd, doublereal *cor_x__, 
	doublereal *cor_y__, doublereal *cor_ut1__, doublereal *cor_lod__)
{
    /* Initialized data */

    static integer narg[426]	/* was [71][6] */ = { 1,1,1,1,1,1,1,1,1,1,1,1,
	    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,
	    2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,-1,-2,-2,0,0,
	    -1,-1,1,1,0,0,0,-2,0,-1,1,1,-1,-1,1,0,0,0,0,0,0,0,0,0,0,0,-1,1,1,
	    0,2,0,0,0,1,1,-3,-1,-2,0,0,-1,-1,-1,-1,1,1,-2,0,0,0,0,-1,1,-1,-1,
	    0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,
	    0,1,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,1,0,1,0,
	    -1,0,0,1,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,-2,-2,-2,-2,-2,-2,-2,-2,-2,
	    -2,-2,-2,0,0,-2,-2,-2,0,0,0,-2,-2,-2,-2,0,0,0,0,0,0,2,0,0,0,0,0,2,
	    2,2,2,2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,0,0,
	    -2,-2,-2,0,0,0,0,0,0,2,-2,0,0,-2,-2,0,0,-2,-2,0,0,0,0,-2,2,0,0,0,
	    0,-2,2,2,2,2,0,0,0,0,0,0,-2,2,0,0,2,0,0,0,0,0,0,0,-2,0,-2,-2,0,0,
	    0,0,-2,-2,2,0,0,0,0,2,0,0,0,2,2,2,0,0,0,0,0,0,0,-2,-1,-2,-1,-2,-1,
	    -2,-1,-2,0,-1,-2,0,0,-2,-1,-2,0,-1,0,-2,-1,-2,-2,0,1,0,-1,-2,0,2,
	    0,0,-1,0,0,2,1,0,2,1,-2,-2,-2,-2,-2,-2,-1,-2,-2,-2,-2,-2,-2,-1,-2,
	    -2,-2,-2,0,-1,-2,-2,-2,1,0,-1,-2,0,-1,2 };
    static real xsin[71] = { -.05f,.06f,.3f,.08f,.46f,1.19f,6.24f,.24f,1.28f,
	    -.28f,9.22f,48.82f,-.32f,-.66f,-.42f,-.3f,-1.61f,-4.48f,-.9f,
	    -.86f,1.54f,-.29f,26.13f,-.22f,-.61f,1.54f,-77.48f,-10.52f,.23f,
	    -.61f,-1.09f,-.69f,-3.46f,-.69f,-.37f,-.17f,-1.1f,-.7f,-.15f,
	    -.03f,-.02f,-.49f,-1.33f,-6.08f,-7.59f,-.52f,.47f,2.12f,-56.87f,
	    -.54f,-11.01f,-.51f,.98f,1.13f,12.32f,-330.15f,-1.01f,2.47f,9.4f,
	    -2.35f,-1.04f,-8.51f,-144.13f,1.19f,.49f,-38.48f,-11.44f,-1.24f,
	    -1.77f,-.77f,-.33f };
    static real xcos[71] = { .94f,.64f,3.42f,.78f,4.15f,4.96f,26.31f,.94f,
	    4.99f,-.77f,25.06f,132.91f,-.86f,-1.72f,-.92f,-.64f,-3.46f,-9.61f,
	    -1.93f,-1.81f,3.03f,-.58f,51.25f,-.42f,-1.2f,3.f,-151.74f,-20.56f,
	    .44f,-1.19f,-2.11f,-1.43f,-7.28f,-1.44f,-1.06f,-.51f,-3.42f,
	    -2.19f,-.46f,-.59f,-.38f,-.04f,-.17f,-1.61f,-2.05f,-.14f,.11f,
	    .49f,-12.93f,-.12f,-2.4f,-.11f,.11f,.11f,1.f,-26.96f,-.07f,-.28f,
	    -1.44f,.37f,.17f,3.5f,63.56f,-.56f,-.25f,19.14f,5.75f,.63f,1.79f,
	    .78f,.62f };
    static real ysin[71] = { -.94f,-.64f,-3.42f,-.78f,-4.15f,-4.96f,-26.31f,
	    -.94f,-4.99f,.77f,-25.06f,-132.9f,.86f,1.72f,.92f,.64f,3.46f,
	    9.61f,1.93f,1.81f,-3.03f,.58f,-51.25f,.42f,1.2f,-3.f,151.74f,
	    20.56f,-.44f,1.19f,2.11f,1.43f,7.28f,1.44f,1.06f,.51f,3.42f,2.19f,
	    .46f,.59f,.38f,.63f,1.53f,3.13f,3.44f,.22f,-.1f,-.41f,11.15f,.1f,
	    1.89f,.08f,-.11f,-.13f,-1.41f,37.58f,.11f,-.44f,-1.88f,.47f,.21f,
	    3.29f,59.23f,-.52f,-.23f,17.72f,5.32f,.58f,1.71f,.75f,.65f };
    static real ycos[71] = { -.05f,.06f,.3f,.08f,.45f,1.19f,6.23f,.24f,1.28f,
	    -.28f,9.22f,48.82f,-.32f,-.66f,-.42f,-.3f,-1.61f,-4.48f,-.9f,
	    -.86f,1.54f,-.29f,26.13f,-.22f,-.61f,1.54f,-77.48f,-10.52f,.23f,
	    -.61f,-1.09f,-.69f,-3.46f,-.69f,-.37f,-.17f,-1.09f,-.7f,-.15f,
	    -.03f,-.02f,.24f,.68f,3.35f,4.23f,.29f,-.27f,-1.23f,32.88f,.31f,
	    6.41f,.3f,-.58f,-.67f,-7.31f,195.92f,.6f,-1.48f,-5.65f,1.41f,.62f,
	    5.11f,86.56f,-.72f,-.29f,23.11f,6.87f,.75f,1.04f,.45f,.19f };
    static real utsin[71] = { .396f,.195f,1.034f,.224f,1.187f,.966f,5.118f,
	    .172f,.911f,-.093f,3.025f,16.02f,-.103f,-.194f,-.083f,-.057f,
	    -.308f,-.856f,-.172f,-.161f,.315f,-.062f,5.512f,-.047f,-.134f,
	    .348f,-17.62f,-2.392f,.052f,-.144f,-.267f,-.288f,-1.61f,-.32f,
	    -.407f,-.213f,-1.436f,-.921f,-.193f,-.396f,-.253f,-.089f,-.224f,
	    -.637f,-.745f,-.049f,.033f,.141f,-3.795f,-.035f,-.698f,-.032f,
	    .05f,.056f,.605f,-16.195f,-.049f,.111f,.425f,-.106f,-.047f,-.437f,
	    -7.547f,.064f,.027f,-2.104f,-.627f,-.068f,-.146f,-.064f,-.049f };
    static real utcos[71] = { -.078f,-.059f,-.314f,-.073f,-.387f,-.474f,
	    -2.499f,-.09f,-.475f,.07f,-2.28f,-12.069f,.078f,.154f,.074f,.05f,
	    .271f,.751f,.151f,.137f,-.189f,.035f,-3.095f,.025f,.07f,-.171f,
	    8.548f,1.159f,-.025f,.065f,.111f,.043f,.187f,.037f,-.005f,-.005f,
	    -.037f,-.023f,-.005f,-.024f,-.015f,-.011f,-.032f,-.177f,-.222f,
	    -.015f,.013f,.058f,-1.556f,-.015f,-.298f,-.014f,.022f,.025f,.266f,
	    -7.14f,-.021f,.034f,.117f,-.029f,-.013f,-.019f,-.159f,0.f,-.001f,
	    .041f,.015f,.002f,.037f,.017f,.018f };

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double d_mod(doublereal *, doublereal *), cos(doublereal), sin(doublereal)
	    ;

    /* Local variables */
    static integer i__, j;
    static doublereal t, ag, dag, arg[6], darg[6], secrad, halfpi;


/*    This subroutine provides, in time domain, the diurnal/subdiurnal */
/*    tidal effets on polar motion ("), UT1 (s) and LOD (s). The tidal terms, */
/*    listed in the program above, have been extracted from the procedure */
/*    ortho_eop.f coed by Eanes in 1997. */

/*    N.B.:  The fundamental lunisolar arguments are those of Simon et al. */

/*    These corrections should be added to "average" */
/*    EOP values to get estimates of the instantaneous values. */

/*     PARAMETERS ARE : */
/*     rjd      - epoch of interest given in mjd */
/*     cor_x    - tidal correction in x (sec. of arc) */
/*     cor_y    - tidal correction in y (sec. of arc) */
/*     cor_ut1  - tidal correction in UT1-UTC (sec. of time) */
/*     cor_lod  - tidal correction in length of day (sec. of time) */

/*     coded by Ch. Bizouard (2002), initially coded by McCarthy and */
/*     D.Gambis(1997) for the 8 prominent tidal waves. */
/* Array of the tidal arguments */
/* Array of their time derivative */
    halfpi = 1.5707963267948966;
    secrad = halfpi * 2. / 6.48e5;
/*  Oceanic tidal terms present in x (microas),y(microas),ut1(microseconds) */
/*  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. */
    t = (*rjd - 51544.5) / 36525.;
/* Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega */
/* et leur derivee temporelle */
/* julian century */
/* Computing 2nd power */
    d__1 = t;
/* Computing 3rd power */
    d__2 = t;
    arg[0] = (t * 3164400184.8128662 + 67310.54841 + d__1 * d__1 * .093104 - 
	    d__2 * (d__2 * d__2) * 6.2e-6) * 15. + 6.48e5;
    arg[0] = d_mod(arg, &c_b4) * secrad;
/* Computing 2nd power */
    d__1 = t;
    darg[0] = (t * .18620800000000001 + 3164400184.8128662 - d__1 * d__1 * 
	    1.8599999999999998e-5) * 15.;
    darg[0] = darg[0] * secrad / 36525.;
/* rad/day */
/* Computing 4th power */
    d__1 = t, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = t;
/* Computing 2nd power */
    d__3 = t;
    arg[1] = d__1 * d__1 * -2.447e-4 + d__2 * (d__2 * d__2) * .051635 + d__3 *
	     d__3 * 31.8792 + t * 1717915923.2178 + 485868.249036;
    arg[1] = d_mod(&arg[1], &c_b4) * secrad;
/* Computing 3rd power */
    d__1 = t;
/* Computing 2nd power */
    d__2 = t;
    darg[1] = d__1 * (d__1 * d__1) * -9.7879999999999994e-4 + d__2 * d__2 * 
	    .15490500000000001 + t * 63.758400000000002 + 1717915923.2178;
    darg[1] = darg[1] * secrad / 36525.;
/* rad/day */
/* Computing 4th power */
    d__1 = t, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = t;
/* Computing 2nd power */
    d__3 = t;
    arg[2] = d__1 * d__1 * -1.149e-5 - d__2 * (d__2 * d__2) * 1.36e-4 - d__3 *
	     d__3 * .5532 + t * 129596581.0481 + 1287104.79305;
    arg[2] = d_mod(&arg[2], &c_b4) * secrad;
/* Computing 3rd power */
    d__1 = t;
/* Computing 2nd power */
    d__2 = t;
    darg[2] = d__1 * (d__1 * d__1) * -4.596e-5 - d__2 * d__2 * 4.08e-4 - t * 
	    1.1064000000000001 + 129596581.0481;
    darg[2] = darg[2] * secrad / 36525.;
/* rad/day */
/* Computing 4th power */
    d__1 = t, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = t;
/* Computing 2nd power */
    d__3 = t;
    arg[3] = d__1 * d__1 * 4.17e-6 - d__2 * (d__2 * d__2) * .001037 - d__3 * 
	    d__3 * 12.7512 + t * 1739527262.8478 + 335779.526232;
    arg[3] = d_mod(&arg[3], &c_b4) * secrad;
/* Computing 3rd power */
    d__1 = t;
/* Computing 2nd power */
    d__2 = t;
    darg[3] = d__1 * (d__1 * d__1) * 1.668e-5 - d__2 * d__2 * 
	    .0031110000000000001 - t * 25.502400000000002 + 1739527262.8478;
    darg[3] = darg[3] * secrad / 36525.;
/* rad/day */
/* Computing 4th power */
    d__1 = t, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = t;
/* Computing 2nd power */
    d__3 = t;
    arg[4] = d__1 * d__1 * -3.169e-5 + d__2 * (d__2 * d__2) * .006593 - d__3 *
	     d__3 * 6.3706 + t * 1602961601.209 + 1072260.70369;
    arg[4] = d_mod(&arg[4], &c_b4) * secrad;
/* Computing 3rd power */
    d__1 = t;
/* Computing 2nd power */
    d__2 = t;
    darg[4] = d__1 * (d__1 * d__1) * -1.2676000000000001e-4 + d__2 * d__2 * 
	    .019778999999999998 - t * 12.741199999999999 + 1602961601.209;
    darg[4] = darg[4] * secrad / 36525.;
/* rad/day */
/* Computing 4th power */
    d__1 = t, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = t;
/* Computing 2nd power */
    d__3 = t;
    arg[5] = d__1 * d__1 * -5.939e-5 + d__2 * (d__2 * d__2) * .007702 + d__3 *
	     d__3 * 7.4722 - t * 6962890.2665 + 450160.398036;
    arg[5] = d_mod(&arg[5], &c_b4) * secrad;
/* Computing 3rd power */
    d__1 = t;
/* Computing 2nd power */
    d__2 = t;
    darg[5] = d__1 * (d__1 * d__1) * -2.3756e-4 + d__2 * d__2 * 
	    .023105999999999998 + t * 14.9444 - 6962890.2665;
    darg[5] = darg[5] * secrad / 36525.;
/* CORRECTIONS */
/* rad/day */
    *cor_x__ = 0.;
    *cor_y__ = 0.;
    *cor_ut1__ = 0.;
    *cor_lod__ = 0.;
    for (j = 1; j <= 71; ++j) {
	ag = 0.;
	dag = 0.;
	for (i__ = 1; i__ <= 6; ++i__) {
	    ag += (doublereal) narg[j + i__ * 71 - 72] * arg[i__ - 1];
	    dag += (doublereal) narg[j + i__ * 71 - 72] * darg[i__ - 1];
	}
	d__1 = halfpi * 4.;
	ag = d_mod(&ag, &d__1);
	*cor_x__ = *cor_x__ + (doublereal) xcos[j - 1] * cos(ag) + (
		doublereal) xsin[j - 1] * sin(ag);
	*cor_y__ = *cor_y__ + (doublereal) ycos[j - 1] * cos(ag) + (
		doublereal) ysin[j - 1] * sin(ag);
	*cor_ut1__ = *cor_ut1__ + (doublereal) utcos[j - 1] * cos(ag) + (
		doublereal) utsin[j - 1] * sin(ag);
	*cor_lod__ -= (-((doublereal) utcos[j - 1]) * sin(ag) + (doublereal) 
		utsin[j - 1] * cos(ag)) * dag;
    }
    *cor_x__ *= 1e-6;
/* arcsecond (") */
    *cor_y__ *= 1e-6;
/* arcsecond (") */
    *cor_ut1__ *= 1e-6;
/* second (s) */
    *cor_lod__ *= 1e-6;
/* second (s) */
    return 0;
} /* pmut1_oceans__ */

/* ---------------------------------------------------------------- */
/* Subroutine */ int pm_gravi__(doublereal *rjd, doublereal *cor_x__, 
	doublereal *cor_y__)
{
    /* Initialized data */

    static integer narg[60]	/* was [10][6] */ = { 1,1,1,1,1,1,1,1,1,1,-1,
	    -1,1,0,0,-1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,-2,-2,-2,-2,-2,0,-2,0,0,0,
	    0,0,-2,0,0,0,2,0,0,0,-1,-2,-2,-1,-2,0,-2,0,-1,0 };
    static real xsin[10] = { -.44f,-2.31f,-.44f,-2.14f,-11.36f,.84f,-4.76f,
	    14.27f,1.93f,.76f };
    static real xcos[10] = { .25f,1.32f,.25f,1.23f,6.52f,-.48f,2.73f,-8.19f,
	    -1.11f,-.43f };
    static real ysin[10] = { -.25f,-1.32f,-.25f,-1.23f,-6.52f,.48f,-2.73f,
	    8.19f,1.11f,.43f };
    static real ycos[10] = { -.44f,-2.31f,-.44f,-2.14f,-11.36f,.84f,-4.76f,
	    14.27f,1.93f,.76f };

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double d_mod(doublereal *, doublereal *), cos(doublereal), sin(doublereal)
	    ;

    /* Local variables */
    static integer i__, j;
    static doublereal t, ag, arg[6], secrad, halfpi;


/*    This subroutine provides, in time domain, the diurnal */
/*    lunisolar effet on polar motion (") */

/*    N.B.:  The fundamental lunisolar arguments are those of Simon et al. */

/*    These corrections should be added to "average" */
/*    EOP values to get estimates of the instantaneous values. */

/*     PARAMETERS ARE : */
/*     rjd      - epoch of interest given in mjd */
/*     cor_x    - tidal correction in x (sec. of arc) */
/*     cor_y    - tidal correction in y (sec. of arc) */

/*     coded by Ch. Bizouard (2002) */
/* Array of the tidal arguments */
    halfpi = 1.5707963267948966;
    secrad = halfpi * 2. / 6.48e5;
/*  Diurnal lunisolar tidal terms present in x (microas),y(microas) */
/*  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. */
    t = (*rjd - 51544.5) / 36525.;
/* Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega */
/* et leur derivee temporelle */
/* julian century */
/* Computing 2nd power */
    d__1 = t;
/* Computing 3rd power */
    d__2 = t;
    arg[0] = (t * 3164400184.8128662 + 67310.54841 + d__1 * d__1 * .093104 - 
	    d__2 * (d__2 * d__2) * 6.2e-6) * 15. + 6.48e5;
    arg[0] = d_mod(arg, &c_b4) * secrad;
/* Computing 4th power */
    d__1 = t, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = t;
/* Computing 2nd power */
    d__3 = t;
    arg[1] = d__1 * d__1 * -2.447e-4 + d__2 * (d__2 * d__2) * .051635 + d__3 *
	     d__3 * 31.8792 + t * 1717915923.2178 + 485868.249036;
    arg[1] = d_mod(&arg[1], &c_b4) * secrad;
/* Computing 4th power */
    d__1 = t, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = t;
/* Computing 2nd power */
    d__3 = t;
    arg[2] = d__1 * d__1 * -1.149e-5 - d__2 * (d__2 * d__2) * 1.36e-4 - d__3 *
	     d__3 * .5532 + t * 129596581.0481 + 1287104.79305;
    arg[2] = d_mod(&arg[2], &c_b4) * secrad;
/* Computing 4th power */
    d__1 = t, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = t;
/* Computing 2nd power */
    d__3 = t;
    arg[3] = d__1 * d__1 * 4.17e-6 - d__2 * (d__2 * d__2) * .001037 - d__3 * 
	    d__3 * 12.7512 + t * 1739527262.8478 + 335779.526232;
    arg[3] = d_mod(&arg[3], &c_b4) * secrad;
/* Computing 4th power */
    d__1 = t, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = t;
/* Computing 2nd power */
    d__3 = t;
    arg[4] = d__1 * d__1 * -3.169e-5 + d__2 * (d__2 * d__2) * .006593 - d__3 *
	     d__3 * 6.3706 + t * 1602961601.209 + 1072260.70369;
    arg[4] = d_mod(&arg[4], &c_b4) * secrad;
/* Computing 4th power */
    d__1 = t, d__1 *= d__1;
/* Computing 3rd power */
    d__2 = t;
/* Computing 2nd power */
    d__3 = t;
    arg[5] = d__1 * d__1 * -5.939e-5 + d__2 * (d__2 * d__2) * .007702 + d__3 *
	     d__3 * 7.4722 - t * 6962890.2665 + 450160.398036;
    arg[5] = d_mod(&arg[5], &c_b4) * secrad;
/* CORRECTIONS */
    *cor_x__ = 0.;
    *cor_y__ = 0.;
    for (j = 1; j <= 10; ++j) {
	ag = 0.;
	for (i__ = 1; i__ <= 6; ++i__) {
	    ag += (doublereal) narg[j + i__ * 10 - 11] * arg[i__ - 1];
	}
	d__1 = halfpi * 4.;
	ag = d_mod(&ag, &d__1);
	*cor_x__ = *cor_x__ + (doublereal) xcos[j - 1] * cos(ag) + (
		doublereal) xsin[j - 1] * sin(ag);
	*cor_y__ = *cor_y__ + (doublereal) ycos[j - 1] * cos(ag) + (
		doublereal) ysin[j - 1] * sin(ag);
    }
    *cor_x__ *= 1e-6;
/* arcsecond (") */
    *cor_y__ *= 1e-6;
/* arcsecond (") */
    return 0;
} /* pm_gravi__ */

