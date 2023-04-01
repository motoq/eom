/******************************************************************************
 *  Copyright (C) 1996-97 Gim J. Der.  All rights reserved. 
 *  ------------------------------------------------------------------------  *
 *                                                                            *
 *  Vinti6()...Vinti Solution with Universal Variables.                       *
 *                                                                            *
 *  Originators:                Albert T. Monuki and Gim J. Der, 1996         *
 *  Fortran to C Conversion:    Herb Reynolds, August 5, 1997                 *
 *                                                                            *
 *  ------------------------------------------------------------------------  *
 *                                                                            *
 *  This function computes a Vinti analytic solution to the equations of      *
 *  motion for drag-free circular, elliptic, parabolic and hyperbolic         *
 *  satellite orbits and ballistic trajectories.                              *
 *                                                                            *
 *  ------------------------------------------------------------------------  *
 *                                                                            *
 *  The Keplerian solution is found by the routine Kepler1().                 *
 *                                                                            *
 *  The implementaion of Vinti theory is based on Getchell's method           *
 *  of mechanization (J. of Spacecraft, Vol. 7, N0. 4, April 1970).           *
 *  Improvements are made to give a fool-proof Vinti solution free of         *
 *  singularities by Der and Monuki.                                          *
 *                                                                            *
 *  The J2 and J3 vaules are those of the wgs84 Earth Gravity Model.          *
 *  Metric derived units (km, km/s) are used only for input / output,         *
 *  while Astronomical units (normalized units) are used internally           *
 *  to reduce numerical round-off errors.                                     *
 *                                                                            *
 *  Limitations: None, for trajectory at least 1,000,000 km from Earth.       *
 *               (mean radius of Moon's orbit is 384,000 km)                  *
 *               The limit is artificial since an unbounded trajectory        *
 *               passes the sphere of influence of the Earth.                 *
 *                                                                            *
 *  Input variables:                                                          *
 *    vt0   = initial time (seconds)                                          *
 *    x0[6] = state vector, position and velocity, (km, km/s)                 *
 *    vt1   = final time (seconds)                                            *
 *                                                                            *
 *  Output variables:                                                         *
 *    x1[6] = state vector, position and velocity, (km, km/s)                 *
 *    oe[6] = Vinti mean elements                                             *
 *                                                                            *
 *                                                                            *
 * Modifications history since release                                        *
 *   Kurt Motekew  December, 2021                                             *
 *     - Removed causes of compile time warnings (mostly unused variables)    *
 *     - Updated a few variables to be constants                              *
 *     - The goal was to make as few changes as possible - no algorithm or    *
 *       code structure changes were made, so sticking with the function      *
 *       names Vinti6 and Kepler1 with the hope these minor updates are not   *
 *       viewed as disrespectful towards the original codebase.                *
 *****************************************************************************/

#include <astro_vinti_prop.h>

#include <string>
#include <array>
#include <memory>

#include <Eigen/Dense>

#include <Vinti.h>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>

static void vinti_local(const double planet[4],
                        double t0, const double x0[6],
                        double tf, double x1[6], double oe[6]);

namespace eom {


VintiProp::VintiProp(const std::string& orbit_name,
                     const JulianDate& epoch,
                     const Eigen::Matrix<double, 6, 1>& xeci,
                     const std::shared_ptr<const EcfEciSys>& ecfeciSys)
{
  name = orbit_name;
  jd0 = epoch;
  ecfeci = ecfeciSys;
    // A true equator ECI frame is required for propagation
  Eigen::Matrix<double, 6, 1> xecf = ecfeci->eci2ecf(jd0,
                                                     xeci.block<3,1>(0,0),
                                                     xeci.block<3,1>(3,0));
  Eigen::Matrix<double, 6, 1> xteme = ecfeci->ecf2teme(jd0,
                                                       xecf.block<3,1>(0,0),
                                                       xecf.block<3,1>(3,0));
  x0[0] = xteme(0);
  x0[1] = xteme(1);
  x0[2] = xteme(2);
  x0[3] = xteme(3);
  x0[4] = xteme(4);
  x0[5] = xteme(5);
}


Eigen::Matrix<double, 6, 1> VintiProp::getStateVector(const JulianDate& jd,
                                                      EphemFrame frame) const
{
  std::array<double, 6> oe;
  std::array<double, 6> x1;
  double t1 {phy_const::tu_per_day*(jd - jd0)};
  vinti_local(planet.data(), 0.0, x0.data(), t1, x1.data(), oe.data());

  Eigen::Matrix<double, 6, 1> xteme;
  xteme(0) = x1[0];
  xteme(1) = x1[1];
  xteme(2) = x1[2];
  xteme(3) = x1[3];
  xteme(4) = x1[4];
  xteme(5) = x1[5];
  Eigen::Matrix<double, 6, 1> xecf = ecfeci->teme2ecf(jd,
                                                      xteme.block<3,1>(0,0),
                                                      xteme.block<3,1>(3,0));

  if (frame == EphemFrame::eci) {
    return ecfeci->ecf2eci(jd, xecf.block<3,1>(0,0), xecf.block<3,1>(3,0));
  }
  return xecf;
}


Eigen::Matrix<double, 3, 1> VintiProp::getPosition(const JulianDate& jd,
                                                   EphemFrame frame) const
{
  std::array<double, 6> oe;
  std::array<double, 6> x1;
  double t1 {phy_const::tu_per_day*(jd - jd0)};
  vinti_local(planet.data(), 0.0, x0.data(), t1, x1.data(), oe.data());

  Eigen::Matrix<double, 3, 1> xteme;
  xteme(0) = x1[0];
  xteme(1) = x1[1];
  xteme(2) = x1[2];
  Eigen::Matrix<double, 3, 1> xecf = ecfeci->teme2ecf(jd, xteme);

  if (frame == EphemFrame::eci) {
    return ecfeci->ecf2eci(jd, xecf);
  }
  return xecf;
}


}


static void vinti_local(const double planet[4],
                        double t0, const double x0[6],
                        double tf, double x1[6], double oe[6])
{
   int i, icf, icg, ick;   /* Loop indices */

   double pin[3], vin[3], pf[3], vf[3];
   double xhat0;
   double r1, r2;

   /* Steps 1 - 3 variables */
   double csq, d0, alph0 ;
   double r02, zpdelta, rhotemp, rho0, rho02, sigma0, rrd, csqsig0, rcs, v, v0;    
   double alph3, alph32, alph1, sqrf, sqrg, alph22, alph2;
   double gamma0, csgam0, p0, s0, pcsgam0, csqs0p0;
   double a1p, b1, a1;
   double gam1, pgam1, gamma, p, dela1;
   double e02;
   double smgam;
   double s1p, px, s1;
   double p0s1, q1, p1;
   double dels1;
   double h0[3], h02, r0mag, v0mag, alp0, a0, e0;

   /* Step 4 variables */
   double gams3, s, pxs, q;   
   double xinc;
   double q2, q4, betad, beta, betasq, g, a, b; 
   double g2, asq, bsq, ab, d5, xm, xk1, d1, ecc2, ecc, rho1;
   double a1sq, b1sq;

   double a2, a3, a4, a5, a6;
   double e3, e4, e5, e6;
   double psq;
   double p3, p4, p5, p6;

   double x21, x22;
   double x33, x32, x31;
   double x44, x43, x42, x41;
   double x55, x54, x53, x52, x51;
   double x66, x65, x64, x63, x62, x61;
   double x77, x76, x75, x74, x73, x72, x71;
   double gg1si, gg1psi;
   double cr11, cr12, cr13, cr14, cr15, cr16, cr17;
   double cr21, cr22, cr23, cr24, cr25, cr26, cr27;
   double cr31, cr32, cr33, cr34, cr35;
   double bmg, bpg, d1ma, d1pa, beta1, beta2, b12, b13, b22, b23, xmm1;
   double xmm2;
   double d3, d4;
   double d1md3, bmag;
   double d10, d20, dd2, dd3, dd4, dd5, dd6;
   double c15, c14, c13, c12, c11, c10, c25, c24, c23, c22, c21, c20;
   double b1q, b1q2, b1q4, b2q, b2q2, b2q4;
   double xk12, xk13;
   double sq, sq2, sq3, ucf1, ucf2, ucf3, denystt, denyst;
   double d1a2, cn11, cn12, cn13, cn14, cn15, cn16, cn17;
   double cn22, cn24, cn26;
   double d41ma, d41pa;
   double cn31, cn32, cn33, cn34, cn35, cn36;
   double u, d5sq, ut_n, ut_d;
   double csu, snu, snu2, snu4;
   double t1, t2, t3, t4, t5, t6;
   double xhat, shat, cacs, acs, chat, zz;  

   /* Step 5 variables */
   double eca;
   double sneca, cneca, s1mes, temp3, tra;
   double snw, cnw;
   double ecccnw;
   double dwdx, dw1dx, dw2dx, dw3dx, dw4dx;
   double v3, v4, v5, v6, v7; 
   double w1, w2, w3, w4, w5, w6;
   double uhat;
   double sstar, cstar, cb1qs, psi1;
   double cb2qs, psi2;
   double r3;
   double en3, en1, en2;
   double somega, capt, deltat, comega;

   /* Step 6 variables */
   double cn1r1, psixhat, dpsx1, dpsx, delx;
   double ystar;

   /* Step 7 variables */
   double qsnu, alphaf, rhof, sigmaf;
   double rhof2, sigf2;
   double rhocsq, df, snalp, csalp;
   double rhoqd, drhofp, drhof, temp12, dsigf;
   double temp20, temp21, dd, dalphf;

   if(fabs(tf - t0) < phy_const::epsdt)      /* 4/01/23 */  // kam
   {                                         /* 8/20/97 */
      for(i = 0; i < 6; i++) x1[i] = x0[i];  /* 8/20/97 */
                                             /* 8/20/97 */
      return;                                /* 8/20/97 */
   }                                         /* 8/20/97 */

   /*
    * Compute 
    *    a. The initial guess of the universal variable xhat at tf.
    *    b. The Kepler solution as a default for the focal circle case. 
    */
   Kepler1(planet, t0, x0, tf, x1, &xhat0);

   /*
    *  Check Vinti's forbidden zone 
    */
   h0[0] = x0[1]*x0[5] - x0[2]*x0[4];
   h0[1] = x0[2]*x0[3] - x0[0]*x0[5];
   h0[2] = x0[0]*x0[4] - x0[1]*x0[3];
   h02   = h0[0]*h0[0] + h0[1]*h0[1] + h0[2]*h0[2];

   r0mag = sqrt(x0[0]*x0[0] + x0[1]*x0[1] + x0[2]*x0[2]);
   v0mag = sqrt(x0[3]*x0[3] + x0[4]*x0[4] + x0[5]*x0[5]);

   alp0  = 2.0/r0mag - v0mag*v0mag/phy_const::gm;

   if (fabs(alp0) < 1.0e-15)
   {
     r1 = h02/(2.0*phy_const::gm);	  /* Parabolic */
   }
   else
   {
     a0 = 1.0/alp0;                   /* Others */
     e02 = 1 - alp0*h02/phy_const::gm;

     //if(e02 < 0) e0 = 0;   /* For most machines  this test is unnececessary */
     //else        e0 = sqrt(1.0 - alp0*h02/gm);
     
     // Above est replaced 9/5/97
     if(e02 <= 0) e0 = 0;    /* For most machines  this test is unnececessary */
     else         e0 = sqrt(e02);

     r1 = a0*(1.0 - e0);
   }

	 if (r1 < 210.0*phy_const::du_per_km)       // converted to du from km, kam
   {
      return;
   }

   for(i = 0; i < 3; i++)
   {
      pin[i] = x0[i];
      vin[i] = x0[i+3];
   }

   /*   
    * Step 1. Initial coordinate transformation
    */   
   csq   =  phy_const::j2;
   d0    =  sqrt( pin[0]*pin[0] + pin[1]*pin[1] ) ;
   alph0 =  atan2(pin[1], pin[0]);

   if(alph0 < 0)
   {
      alph0 = utl_const::tpi + alph0;
   }

   r02     =  d0*d0 + pin[2]*pin[2];
   zpdelta =  pin[2];
   rhotemp =  r02 - csq;
   rho0    =  sqrt(rhotemp + sqrt(rhotemp*rhotemp + 4*csq*(zpdelta*zpdelta)))/sqrt(2.0);
   rho02   =  rho0*rho0;
   sigma0  =  zpdelta/rho0;
   rrd     =  pin[0]*vin[0] + pin[1]*vin[1] + pin[2]*vin[2];
   csqsig0 =  csq*sigma0;
   rcs     =  rho02 + csqsig0*sigma0;
   v       = -rho0/rcs;
   v0      =  sqrt( vin[0]*vin[0] + vin[1]*vin[1] + vin[2]*vin[2] );

   /*
    * Step 2. The first three Jacobi constants, alph1, alph2, alph3
    */ 
   alph3   = pin[0]*vin[1] - pin[1]*vin[0];
   alph32  = alph3*alph3;
   alph1   = 0.5*v0*v0 + v;
   sqrf    = rho0*rrd + csqsig0*vin[2];
   sqrg    =-sigma0*rrd + rho0*vin[2];
   alph22  = 2*rho0 + 2*alph1*rho02 + (csq*alph32 - sqrf*sqrf)/(rho02 + csq);
   alph2   = sqrt(alph22);
   gamma0  = 2*alph1;
   csgam0  = csq*gamma0;
   p0      = alph22;
   s0      = 1 - alph32/alph22;
   pcsgam0 = p0 - csgam0;
   csqs0p0 = csq*s0*p0;

   /*
    *  Step 3. Factorizing the F and G quartics
    * 
    *  Vinti(1961) Equations (3.19) to (3.22)
    *  Getchell(1970) Equations (6) to (9)
    *  Factorization by numerical method is more accurate than
    *  analytic approximation.
    *
    *  alph1 < 0, gamma0 = 2*alph1 < 0  for ellipse and circle
    *  alph1 = 0, gamma0 = 0            for parabola
    *  alph1 > 0, gamma0 = 2*alph1 > 0  for hyperbola
    *  Note that even if the input state is a parabola,
    *  alph1 may not be zero.
    *
    *  Factorizing the F(rho) quartic
    */
   a1p = 0;
   b1 = csqs0p0/pcsgam0;
   a1 =(csq - b1)/pcsgam0;
   
   for(icf = 1; icf <= 5; icf++)
   {
      gam1  = 1 + gamma0*a1;
      pgam1 = pcsgam0 + b1*gamma0 - 4*a1*gam1;
      gamma = gamma0/gam1;
      p     = pgam1/gam1;
      b1    = csqs0p0/pgam1;
      a1    =(csq - gam1*b1)/pgam1;
      dela1 = a1 - a1p;
      
      if(fabs(dela1) < 1.0e-15)  break;
      a1p = a1;
   }

   oe[0] = p;     // Vinti mean element (conic parameter)

   /*   
    * Check on some critical parameters in SI derived units
    */   

   e02 = 1 + 2*alph1*alph22;

   /*	r1pr2_vinti = -2.d0/gamma	            // gamma=0 for parabola 
    *	r1xr2_vinti = -p/gamma
    *	rho1 = p/2.d0		 ! parabola
    *	ax = -ae/gamma
    * e0 = dsqrt(1.d0 + 2.d0*alph1*alph22)   // problem sometimes 
    */

   if (gamma < 0) smgam = sqrt(-gamma);
	else           smgam = sqrt(gamma);

   /*
    * Factorizing the G(sigma) quartic
    */
   s1p = 0;
   px = 0;
   s1 = 1;

   for(icg = 1; icg <= 5; icg++)
   {
      p0s1 = p0*s1;
      q1 = -csgam0/p0s1;
      p1 = 2.0*q1*px;
      px = s0*p1/(2.0*s1);
      s1 = (pcsgam0 - s0*p0*q1)/((1 - 2*px*p1)*p0);
      dels1 = s1 - s1p;
      
      if (fabs(dels1) < 1.0e-15)  break;
      s1p = s1;
   }

  /*
   * Step 4. Initialize the coefficients of the R and N integrals.
   * First, determine the constants after the factorization, inclination, 
   * eccentricity,rho1, etc.
   */
   gams3 = gamma*smgam;
   s = s0/s1;
	pxs = px*px + s;     oe[2] = pxs;   // Vinti mean element ( sin^2 (I) )

	if (pxs < 0) q = 0;
	else         q = sqrt(pxs);
 
   xinc = asin(q);

   if( alph3*cos(xinc) < 0 )  xinc = utl_const::pi - xinc;

   q2 = q*q;
   q4 = q2*q2;
   betad = 1 + p1*px - q1*px*px - q1*q2;
   beta  =(p1 - 2*q1*px)/betad;
   betasq = beta*beta;
   g = -beta/(1 + sqrt(1 - betasq*q2));
   a =  px + g*q2;
   b =  1 + g*px;

   g2  = g*g;
   asq = a*a;
   bsq = b*b;
   ab  = a*b;
   d5  = 1 + p1*a - q1*asq;
   xm  = (q1*bsq - p1*b*g - g2)/d5;
   xk1 = xm*q2;
   d1  = sqrt((1 - g2*q2)/(s1*d5));
   ecc2 = 1 + p*gamma;
   ecc  = sqrt(ecc2);      oe[1] = ecc;   // Vinti mean element (eccentricity)
   rho1 = p/(1 + ecc);

  /*
   * Coefficients for the R - integrals
   * A0 = 1, A1 from factorization, A2 to A6
   */
   a1sq = a1*a1;
   b1sq = b1*b1;

   a2 =(3*a1sq - b1)/2;
   a3 =(2.5*a1sq*a1 - 1.50*a1*b1);
   a4 = 0.375*(b1sq - 10*a1sq*b1);  
   a5 = 1.875*a1*b1sq;
   a6 = -0.3125*b1sq*b1;

  /*
   *  W1, W2 ... W6 for the R - integrals
   */
   e3 = ecc2*ecc;
   e4 = ecc2*ecc2;
   e5 = e3*ecc2;
   e6 = e3*e3;

   psq = p*p;

   p3 = psq*p;
   p4 = psq*psq;
   p5 = p3*psq;
   p6 = p4*psq;

   x21 = 1/p;
   x22 = ecc/p;

   x33 = 0.5*ecc2/psq;
   x32 = 2*ecc/psq;
   x31 = 1/psq + x33;

   x44 = e3/(3*p3);
   x43 = 1.5*ecc2/p3;
   x42 = 3*ecc/p3 + 2*x44;
   x41 = 1/p3 + x43;

   x55 = x33*x33;
   x54 = 4*x44*x21;
   x53 = 3*ecc2/p4 + 1.5*x55;
   x52 = 4*ecc/p4 + 2*x54;
   x51 = 1/p4 + x53;

   x66 = 0.2*e5/p5;
   x65 = 1.25*e4/p5;
   x64 = 10*e3/(3*p5) + 4*x66/3;
   x63 = 5*ecc2/p5 + 1.5*x65;
   x62 = 5*ecc/p5 + 2*x64;
   x61 = 1/p5 + x63;

   x77 = e6/(6*p6);
   x76 = 1.2*e5/p6;
   x75 = 3.75*e4/p6 + 1.25*x77;
   x74 = 20*e3/(3*p6) + x76/0.75;
   x73 = 7.5*ecc2/p6 + 1.5*x75;
   x72 = 6*ecc/p6 + 2*x74;
   x71 = 1/p6 + x73;	         
   /* Coefficients for the W's completed */


   gg1si = 1.0/sqrt(gam1);	   /* gam1 > 0 always */
   gg1psi = gg1si/sqrt(p);	   /* p > 0 always */

   /* R1 coefficients */
   cr11 = (rho1 + a1)*gg1si;
   cr12 = ecc*gg1si;
   cr13 = a2*gg1psi;
   cr14 = a3*gg1psi;
   cr15 = a4*gg1psi;
   cr16 = a5*gg1psi;
   cr17 = a6*gg1psi;

   /* R2 coefficients */
   cr21 = sqrt(p0/pgam1);
   cr22 = a1*cr21;
   cr23 = a2*cr21;
   cr24 = a3*cr21;
   cr25 = a4*cr21;
   cr26 = a5*cr21;
   cr27 = a6*cr21;

   /* R3 coefficients */
   cr31 = alph3*gg1psi;
   cr32 = a1*cr31;
   cr33 = (a2 - csq)*cr31;
   cr34 = (a3 - a1*csq)*cr31;
   cr35 = a4*cr31 - csq*cr33;
   /* R coefficients completed */


   /* Coefficients for the N - integrals */

   /*
    *  N3 coefficients
    */
   bmg   = b - g;
   bpg   = b + g;
   d1ma  = 1 - a;
   d1pa  = 1 + a;
   beta1 = bmg/d1ma;
   beta2 =-bpg/d1pa;
   b12   = beta1*beta1;
   b13   = b12*beta1;
   b22   = beta2*beta2;
   b23   = b22*beta2;
   xmm1  = sqrt(1 - b12*q2);
   
   if(xmm1*alph3 < 0) xmm1 = -xmm1;

   xmm2 = sqrt(1 - b22*q2);

   if(xmm2*alph3 < 0) xmm2 = -xmm2;

   d3 = q1;
   d4 = d1*alph3/(2.*alph2);

   d1md3 = 1 - d3;
   bmag  = b - a*g;

   d10 = bmag/bmg*sqrt(d1md3/(d5*(1.0 - xm/b12)));
   d20 = bmag/bpg*sqrt(d1md3/(d5*(1.0 - xm/b22)));

   dd2 = xm/2;
   dd3 = dd2*g;
   dd4 = 1.5*dd2*dd2;
   dd5 = dd4*g;
   dd6 = dd2*dd4/0.6;

   c15 = dd6/(b13*b13);
   c14 = c15 + dd5/(b12*b13);
   c13 = c14 + dd4/(b12*b12);
   c12 = c13 + dd3/b13;
   c11 = c12 + dd2/b12;
   c10 = c11 + g/beta1;

   c25 = dd6/(b23*b23);
   c24 = c25 + dd5/(b22*b23);
   c23 = c24 + dd4/(b22*b22);
   c22 = c23 + dd3/b23;
   c21 = c22 + dd2/b22;
   c20 = c21 + g/beta2;

   b1q  = beta1*q;
   b1q2 = b1q*b1q;
   b1q4 = b1q2*b1q2;
   b2q  = beta2*q;
   b2q2 = b2q*b2q;
   b2q4 = b2q2*b2q2;     /* N3 coefficients completed */

   /*
    *  N1 and N2 coefficients
    */
   xk12 = xk1*xk1;
   xk13 = xk12*xk1;

   /* Byrd and Friedman formula for elliptic integral */
   sq  = xk1/16 + xk12/32 + 21*xk13/1024;
   sq2 = sq*sq;
   sq3 = sq2*sq;

   ucf1 = 2*sq/(1 + sq2);
   ucf2 = sq2/(1 + sq2*sq2);
   ucf3 = 2.0*sq3/(3.0*(1 + sq3*sq3));

   denystt = 1 + sq + sq;
   denyst  = denystt*denystt*d1;

   /* N1 coefficients */
   d1a2 = d1/alph2;
   cn11 = d1a2*asq;
   cn12 = d1a2*2*ab*q;
   cn13 = d1a2*(bsq - 4*ab*g)*q2;
   cn14 = d1a2*(xm*ab - 2*bsq*g)*q2*q;
   cn15 = d1a2*(3*bsq*g2 + xm*bsq/2)*q4;
   cn16 =-d1a2*xm*bsq*g*q4*q;
   cn17 = 0.375*d1a2*xm*xm*bsq*q4*q2;

   /* N2 coefficients */
   cn22 = 0.5*xk1*d1;
   cn24 = 0.375*xk12*d1;
   cn26 = 0.3125*xk13*d1;

   /* N3 coefficients */
   d41ma = d4/d1ma;
   d41pa = d4/d1pa;

   cn31 = -d41ma*c10 - d41pa*c20;
   cn32 = -d41ma*c11*b1q - d41pa*c21*b2q;
   cn33 = -d41ma*c12*b1q2 - d41pa*c22*b2q2;
   cn34 = -d41ma*c13*b1q2*b1q - d41pa*c23*b2q2*b2q;
   cn35 = -d41ma*c14*b1q4 - d41pa*c24*b2q4;
   cn36 = -d41ma*c15*b1q4*b1q - d41pa*c25*b2q4*b2q;

   /*
    * This is to avoid the singularity at zero inclination
    * u must be determined exactly
    */
   if (q == 0.0)
   {
	   u = 0;
   }
	else
   {
      d5sq = sqrt( s1*p0*(1 + p1*sigma0 - q1*sigma0*sigma0) );
      ut_n = (sigma0 - a)*d5sq;
      ut_d = sqrg*sqrt(1 - g2*q2);

      u = atan2(ut_n, ut_d);
   }

   /*
    * Calculate the Tk = tk, k = 0, T0 = t0 = u.  Here k = 1, ... ,6
    */
   csu  = cos(u);
   snu  = sin(u);
   snu2 = snu*snu;
   snu4 = snu2*snu2;

   t1 = 1 - csu;
   t2 = (u - csu*snu)/2;
   t3 = (2*t1 - csu*snu2)/3;
   t4 = (3*t2 - csu*snu2*snu)/4;
   t5 = (4*t3 - csu*snu4)/5;
   t6 = (5*t4 - csu*snu4*snu)/6;

   /*
    * Compute xhat at the initial time which is used only once to get
    * the Jacobi's constants (beta1= - capt, beta2=somega, beta3=comega)
    * Note that xhat should not be computed from the initial osculating
    * position and velocity vectors which are Kepler's solutions.
    */
   if (gamma < -1.0e-14)      /* Ellipse */
   {
      //if (fabs(ecc) < 1.0e-10 || fabs(sqrf) < 1.0e-10)
      // Replaced 9/5/97
      if (fabs(ecc) < 1.0e-10 || fabs(sqrf) < 1.0e-10 || fabs(rho0 -rho1) < 1e-5 )
      {
         xhat = 0;
      }
      else
      {
         //shat = sign(1, sqrf)*sqrt(gamma*rho02 + 2*rho0 - p)/ecc;
         shat = ( fabs(sqrf)/sqrf )*sqrt(gamma*rho02 + 2*rho0 - p)/ecc;

         cacs =(rho0 - rho1)*gamma/ecc + 1;
         acs  = atan2(shat*smgam, cacs);
      
         if(acs <  0) acs = utl_const::tpi + acs;

         xhat = acs/smgam;
      }   
   }
	else if (gamma > 1.0e-14)  /* Hyperbola */
   {
      if (fabs(ecc) < 1.0e-10 || fabs(sqrf) < 1.0e-10)
      {   
          xhat = 0.0;
      }
      else
      {
        chat = (rho0 - rho1)/ecc;		   /* chat > 0, always */
	     zz = gamma*chat + 1;
        xhat = (  log( zz + sqrt(zz*zz - 1) )  )/sqrt(gamma);

        if (sqrf < 0)   xhat = -xhat;
      }
   }
	else								/* parabola */
   {
	    xhat = rrd;		      /* r*rdot */
   }

   /* .............. End initialization ...................*/


   /*
    * Step 5. The last three Jacobi constants, -capt, somega, comega
    *
    * Determine the true anomaly, tra, from xhat of ti
    */
   eca = smgam*xhat;

   if (gamma < -1.0e-14)      /* Ellipse */
   {
      sneca = sin(eca);
      cneca = cos(eca);
      s1mes = sqrt(1 - ecc2);
      temp3 = 2*atan( (ecc*sneca)/(1 + s1mes - ecc*cneca) );
      tra   = eca + temp3;
   }
   else if (gamma > 1.0e-14)  /* Hyperbola */
   {
      sneca = sinh(eca);
      cneca = cosh(eca);

      chat = (cneca - 1)/gamma;
      shat = sneca/smgam;
      tra = atan2( sqrt(p)*shat,( rho1 - chat) );
   }
   else								/* Parabola */
   {
      tra = 2*atan( xhat/sqrt(p) );
   }

   snw = sin(tra);
   cnw = cos(tra);
   ecccnw = ecc*cnw;

   /*
    *  Is there such a thing call parabolic trajectory in the OSC system? 
    */
   dwdx  = (1 + ecccnw)/sqrt(p); 
   dw1dx = (1 + ecccnw)/p;

   dw2dx = dw1dx*dw1dx;
   dw3dx = dw1dx*dw2dx;
   dw4dx = dw2dx*dw2dx;

   v3 = snw*cnw;
   v4 = v3*cnw;
   v5 = v4*cnw;
   v6 = v5*cnw;
   v7 = v6*cnw;

   w1 = x21*tra + x22*snw;
   w2 = x31*tra + x32*snw + x33*v3;
   w3 = x41*tra + x42*snw + x43*v3 + x44*v4;
   w4 = x51*tra + x52*snw + x53*v3 + x54*v4 + x55*v5;
   w5 = x61*tra + x62*snw + x63*v3 + x64*v4 + x65*v5 + x66*v6;
   w6 = x71*tra + x72*snw + x73*v3 + x74*v4 + x75*v5 + x76*v6 + x77*v7;

   if (fabs(gamma) >= 1.0e-14)	/* Circle, ellipse, hyperbloa */
   {  
      chat = (cneca - 1)/gamma;
      uhat = (sneca - eca)/gams3;
   }
	else								   /*  Parabola */
   {
      chat = xhat*xhat/2;
      uhat = xhat*xhat*xhat/6;
   }

   r1 = cr11*xhat + cr12*uhat + cr13*tra + cr14*w1 + cr15*w2 + cr16*w3 + cr17*w4;
   r2 = cr21*tra  + cr22*w1   + cr23*w2  + cr24*w3 + cr25*w4 + cr26*w5 + cr27*w6;

   sstar = sin(0.5*u);
   cstar = cos(0.5*u);
   cb1qs = cstar - b1q*sstar;
   psi1  = atan(xmm1*sstar/cb1qs);

   if(cb1qs*cos(psi1) < 0.0) psi1 = utl_const::pi + psi1;

   cb2qs = cstar - b2q*sstar;
   psi2  = atan(xmm2*sstar/cb2qs);

   if(cb2qs*cos(psi2) < 0.0) psi2 = utl_const::pi + psi2;

   r3 = cr31*w2 + cr32*w3 + cr33*w4 + cr34*w5 + cr35*w6;

   en3 = d10*psi1 + d20*psi2 + cn31*u + cn32*t1 + cn33*t2 + cn34*t3 + cn35*t4 + cn36*t5;
   en1 = cn11*u + cn12*t1 + cn13*t2 + cn14*t3 + cn15*t4 + cn16*t5 + cn17*t6;
   en2 = d1*u + cn22*t2 + cn24*t4 + cn26*t6;
   
   somega  = en2 - r2;		            /* beta2 */
   capt    = t0 - r1 - csq*en1;		   /* -beta1 */
   deltat  = tf - capt;			  
   comega  = alph0 + csq*r3 - en3;	   /* beta3 */

   oe[3] = -capt;                   /* Vinti mean element  "beta1" */
   oe[4] = somega;                  /* Vinti mean element  "beta2" */
   oe[5] = comega;                  /* Vinti mean element  "beta3" */

   /*
    * Step 6. At tf, solve the kinematical equations for the oblate
    * spheroidal coordinates (rho, sigma, phi)
    *
    * The generalized Kepler equation is solved by iteration
    * initial guess of xhat is that of the Kepler's solution at tf 
    */
	//xhat = xhat0;
	xhat += xhat0;  // replaced 9/5/97

   for(ick = 1; ick <= 10; ick++)
   {      
      eca = smgam*xhat;
      
      if (gamma < -1.0e-14)         /* Ellipse */
      {
         sneca = sin(eca);
         cneca = cos(eca);
         s1mes = sqrt(1 - ecc2);
         temp3 = 2*atan( (ecc*sneca)/(1 + s1mes - ecc*cneca) );
         tra   = eca + temp3;
      }
		else if (gamma > 1.0e-14)     /* Hyperbola */
      {
         sneca = sinh(eca);
         cneca = cosh(eca);
         chat  = (cneca - 1)/gamma;
         shat  = sneca/smgam;
         tra   = atan2( sqrt(p)*shat, (rho1 - chat) );
      }
		else									/* Parabola */
      {
         tra = 2*atan(xhat/sqrt(p));
      }
      
      snw = sin(tra);
      cnw = cos(tra);
      ecccnw = ecc*cnw;
      
      dwdx  = (1 + ecccnw)/sqrt(p);
      dw1dx = (1 + ecccnw)/p;
      dw2dx = dw1dx*dw1dx;
      dw3dx = dw1dx*dw2dx;
      dw4dx = dw2dx*dw2dx;
      
      v3 = snw*cnw;
      v4 = v3*cnw;
      v5 = v4*cnw;
      v6 = v5*cnw;
      v7 = v6*cnw;
      
      w1 = x21*tra + x22*snw;
      w2 = x31*tra + x32*snw + x33*v3;
      w3 = x41*tra + x42*snw + x43*v3 + x44*v4;
      w4 = x51*tra + x52*snw + x53*v3 + x54*v4 + x55*v5;
      w5 = x61*tra + x62*snw + x63*v3 + x64*v4 + x65*v5 + x66*v6;
      w6 = x71*tra + x72*snw + x73*v3 + x74*v4 + x75*v5 + x76*v6 + x77*v7;
      
      if(fabs(gamma) >= 1.0e-14)
      {
         
         chat =(cneca - 1)/gamma;
         uhat =(sneca - eca)/gams3;
      }
		else		/* Parabola */
      {
         chat = xhat*xhat/2;
         uhat = xhat*xhat*xhat/6;
      }
      
      r1 = cr11*xhat + cr12*uhat + cr13*tra + cr14*w1 + cr15*w2 + cr16*w3 + cr17*w4;
      r2 = cr21*tra  + cr22*w1   + cr23*w2  + cr24*w3 + cr25*w4 + cr26*w5 + cr27*w6;
      
      //if (iflag == 1) goto LABLE_60; 	/* solution converged */
      
      ystar = (r2 + somega)/denyst;
      u = ystar + ucf1*sin(2*ystar) + ucf2*sin(4*ystar) + ucf3*sin(6*ystar);
        
      csu  = cos(u);
      snu  = sin(u);
      snu2 = snu*snu;
      snu4 = snu2*snu2;
      
      t1 = 1 - csu;
      t2 = (u - csu*snu)/2;
      t3 = (2*t1 - csu*snu2)/3;
      t4 = (3*t2 - csu*snu2*snu)/4;
      t5 = (4*t3 - csu*snu4)/5;
      t6 = (5*t4 - csu*snu4*snu)/6;
      
      en1 = cn11*u + cn12*t1 + cn13*t2 + cn14*t3 + cn15*t4 + cn16*t5 + cn17*t6;
      
      cn1r1  = csq*en1 + r1;
      psixhat = cn1r1 - deltat;		   	    /* Function */
      
      dpsx1 = cr13 + cr14*dw1dx + cr15*dw2dx + cr16*dw3dx + cr17*dw4dx;
      dpsx = cr11 + cr12*chat + dpsx1*dwdx; 	 /* 1st derivative */
      delx = psixhat/dpsx;
      
      xhat = xhat - delx;
      
      if(fabs(delx) < 1.0e-14 )  break;
   }
      
   /*
    * Step 7. Final coordinate transformation, after solution converged
    * OSC state vector is (rhof, sigmaf, alphaf, drhof, dsigf, dalphf)
    * ECI state vector is (pf, vf)
    */
   sstar = sin(0.5*u);
   cstar = cos(0.5*u);
   cb1qs = cstar - b1q*sstar;
   psi1  = atan(xmm1*sstar/cb1qs);

   if(cb1qs*cos(psi1) < 0) psi1 = utl_const::pi + psi1;

   cb2qs = cstar - b2q*sstar;
   psi2  = atan(xmm2*sstar/cb2qs);

   if(cb2qs*cos(psi2) < 0) psi2 = utl_const::pi + psi2;

   r3 = cr31*w2 + cr32*w3 + cr33*w4 + cr34*w5 + cr35*w6;
   en3 = d10*psi1 + d20*psi2 + cn31*u + cn32*t1 + cn33*t2 + cn34*t3 + cn35*t4 + cn36*t5;
   qsnu = q*snu;
   alphaf = comega - csq*r3 + en3;
   rhof   = rho1 + ecc*chat;
   sigmaf = (a + b*qsnu)/(1 + g*qsnu);

   pf[2] = rhof*sigmaf;

   rhof2 = rhof*rhof;
   sigf2 = sigmaf*sigmaf;

   rhocsq = rhof2 + csq;
   df = sqrt(rhocsq*(1 - sigf2));
   snalp = sin(alphaf);
   csalp = cos(alphaf);

   pf[0] = df*csalp;
   pf[1] = df*snalp;

   shat = xhat + gamma*uhat;
   rcs  = rhof2 + csq*sigf2;

   rhoqd  = rhof2 - 2*a1*rhof + b1;
   drhofp = sqrt(gam1*rhoqd)/rcs;
   drhof  = ecc*shat*drhofp;
   temp12 =(1 + p1*sigmaf - q1*sigf2)*(1 - g2*q2)*s1;
   dsigf  = alph2*q*sqrt(temp12)/rcs*csu/(1 + g*qsnu);
   temp20 =(1 - sigf2)*(rhof*drhof);
   temp21 = rhocsq*sigmaf*dsigf;
   dd =(temp20 - temp21)/df;
   dalphf = alph3/(df*df);

   vf[0] = dd*csalp - pf[1]*dalphf;
   vf[1] = dd*snalp + pf[0]*dalphf;
   vf[2] = sigmaf*drhof + rhof*dsigf;

   /*
    * Change from astronomical units to SI (derived) units
    */
   for(i = 0; i < 3; i++)
   {
      x1[i] = pf[i];
      x1[i+3] = vf[i];
   }
}

