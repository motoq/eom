/*
 *  Significantly restructured version of Vinti6() from the Vinti6.c
 *  source.  Initialization has been localized in the constructor so it
 *  does not have to be repeated for each propagation call.  Support for
 *  non-canonical units has been dropped.  J3 effects have been removed
 *  from this implementation.
 *
 *  Original header file follows, although not completely applicable to
 *  this implementation.  See the included Vinti6.c source for original
 *  implementation.
 *
 *  C to C++ conversion:        Kurt A. Motekew, 2023
 */

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
 *****************************************************************************/

#include <astro_vinti_prop.h>

#include <string>
#include <array>
#include <memory>
#include <stdexcept>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>

namespace {
  constexpr double csq {phy_const::j2};
}

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
   for (int ii=0; ii<3; ++ii) {
     pin[ii] = xteme(ii);
     vin[ii] = xteme(ii+3);
   }

   m_kep = std::make_unique<KeplerProp>("", epoch, xeci, ecfeciSys);

   /*
    *  Check Vinti's forbidden zone 
    */
   std::array<double, 3> h0;
   h0[0] = pin[1]*vin[2] - pin[2]*vin[1];
   h0[1] = pin[2]*vin[0] - pin[0]*vin[2];
   h0[2] = pin[0]*vin[1] - pin[1]*vin[0];
   double h02 {h0[0]*h0[0] + h0[1]*h0[1] + h0[2]*h0[2]};
   double r0mag {sqrt(pin[0]*pin[0] + pin[1]*pin[1] + pin[2]*pin[2])};
   double v0mag {sqrt(vin[0]*vin[0] + vin[1]*vin[1] + vin[2]*vin[2])};
   double alp0 {2.0/r0mag - v0mag*v0mag/phy_const::gm};
   double rlim {};
   if (fabs(alp0) < 1.0e-15) {
     rlim = h02/(2.0*phy_const::gm);              // Parabolic
   } else {
     double a0 {1.0/alp0};                        // Others
     double e02 {1.0 - alp0*h02/phy_const::gm};
     double e0 {};
     if(e02 <= 0) {                               // For most machines
       e0 = 0;                                    // this test is unnececessary
     } else {
       e0 = sqrt(e02);
     }
     rlim = a0*(1.0 - e0);
   }
   if (rlim < 210.0*phy_const::du_per_km) {
     throw std::invalid_argument(
         "VintiProp::VintiProp(): Failure of Vinti's forbidden zone");
   }

   /*   
    * Step 1. Initial coordinate transformation
    */   
   double d0    =  sqrt( pin[0]*pin[0] + pin[1]*pin[1] ) ;
   double alph0 =  atan2(pin[1], pin[0]);

   if(alph0 < 0.0)
   {
      alph0 = utl_const::tpi + alph0;
   }

   double r02     =  d0*d0 + pin[2]*pin[2];
   double zpdelta =  pin[2];
   double rhotemp =  r02 - csq;
   double rho0    =  sqrt(rhotemp + sqrt(rhotemp*rhotemp +
                     4.0*csq*(zpdelta*zpdelta)))/sqrt(2.0);
   double rho02   =  rho0*rho0;
   double sigma0  =  zpdelta/rho0;
   double rrd     =  pin[0]*vin[0] + pin[1]*vin[1] + pin[2]*vin[2];
   double csqsig0 =  csq*sigma0;
   double rcs     =  rho02 + csqsig0*sigma0;
   double v       = -rho0/rcs;
   double v0      =  sqrt( vin[0]*vin[0] + vin[1]*vin[1] + vin[2]*vin[2] );

   /*
    * Step 2. The first three Jacobi constants, alph1, alph2, alph3
    */ 
   alph3   = pin[0]*vin[1] - pin[1]*vin[0];
   double alph32  = alph3*alph3;
   double alph1   = 0.5*v0*v0 + v;
   double sqrf    = rho0*rrd + csqsig0*vin[2];
   double sqrg    =-sigma0*rrd + rho0*vin[2];
   double alph22  = 2.0*rho0 + 2.0*alph1*rho02 + (csq*alph32 - sqrf*sqrf)/
                                                            (rho02 + csq);
   alph2   = sqrt(alph22);
   double gamma0  = 2.0*alph1;
   double csgam0  = csq*gamma0;
   double p0      = alph22;
   double s0      = 1.0 - alph32/alph22;
   double pcsgam0 = p0 - csgam0;
   double csqs0p0 = csq*s0*p0;

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
   double a1p = 0.0;
   b1 = csqs0p0/pcsgam0;
   a1 =(csq - b1)/pcsgam0;
   
   double pgam1 {};
   for(int icf = 1; icf <= 5; icf++)
   {
      gam1  = 1.0 + gamma0*a1;
      pgam1 = pcsgam0 + b1*gamma0 - 4.0*a1*gam1;
      gamma = gamma0/gam1;
      p     = pgam1/gam1;
      b1    = csqs0p0/pgam1;
      a1    =(csq - gam1*b1)/pgam1;
      double dela1 = a1 - a1p;
      
      if (fabs(dela1) < 1.0e-15)
      {
         break;
      }
      a1p = a1;
   }
   oe0[0] = p;                         // Vinti mean element (conic parameter)

   /*   
    * Check on some critical parameters
    */   

   if (gamma < 0.0)
   {
     smgam = sqrt(-gamma);
   }
   else
   {
     smgam = sqrt(gamma);
   }

   /*
    * Factorizing the G(sigma) quartic
    */
   double s1p = 0.0;
   double px = 0.0;
   s1 = 1.0;
   for(int icg = 1; icg <= 5; icg++)
   {
      double p0s1 = p0*s1;
      q1 = -csgam0/p0s1;
      p1 = 2.0*q1*px;
      px = s0*p1/(2.0*s1);
      s1 = (pcsgam0 - s0*p0*q1)/((1.0 - 2.0*px*p1)*p0);
      double dels1 = s1 - s1p;
      
      if (fabs(dels1) < 1.0e-15)
      {
         break;
      }
      s1p = s1;
   }

   /*
    * Step 4. Initialize the coefficients of the R and N integrals.
    * First, determine the constants after the factorization, inclination, 
    * eccentricity,rho1, etc.
    */

   gams3 = gamma*smgam;
   double s = s0/s1;
	 double pxs = px*px + s;
   oe0[2] = pxs;                       // Vinti mean element ( sin^2 (I) )

   if (pxs >= 0.0)
   {
      q = sqrt(pxs);
   }
 
   double xinc = asin(q);

   if (alph3*cos(xinc) < 0.0)
   {
      xinc = utl_const::pi - xinc;
   }

   q2 = q*q;
   double q4 = q2*q2;
   double betad = 1.0 + p1*px - q1*px*px - q1*q2;
   double beta  =(p1 - 2.0*q1*px)/betad;
   double betasq = beta*beta;
   g = -beta/(1.0 + sqrt(1.0 - betasq*q2));
   a =  px + g*q2;
   b =  1.0 + g*px;

   g2  = g*g;
   double asq = a*a;
   double bsq = b*b;
   double ab  = a*b;
   double d5  = 1.0 + p1*a - q1*asq;
   double xm  = (q1*bsq - p1*b*g - g2)/d5;
   double xk1 = xm*q2;
   double d1  = sqrt((1.0 - g2*q2)/(s1*d5));
   ecc2 = 1.0 + p*gamma;
   ecc  = sqrt(ecc2);
   oe0[1] = ecc;                       // Vinti mean element (eccentricity)
   rho1 = p/(1.0 + ecc);

  /*
   * Coefficients for the R - integrals
   * A0 = 1, A1 from factorization, A2 to A6
   */
   double a1sq = a1*a1;
   double b1sq = b1*b1;

   double a2 =(3.0*a1sq - b1)/2.0;
   double a3 =(2.5*a1sq*a1 - 1.50*a1*b1);
   double a4 = 0.375*(b1sq - 10.0*a1sq*b1);
   double a5 = 1.875*a1*b1sq;
   double a6 = -0.3125*b1sq*b1;

  /*
   *  W1, W2 ... W6 for the R - integrals
   */
   double e3 = ecc2*ecc;
   double e4 = ecc2*ecc2;
   double e5 = e3*ecc2;
   double e6 = e3*e3;

   double psq = p*p;

   double p3 = psq*p;
   double p4 = psq*psq;
   double p5 = p3*psq;
   double p6 = p4*psq;

   x21 = 1.0/p;
   x22 = ecc/p;

   x33 = 0.5*ecc2/psq;
   x32 = 2.0*ecc/psq;
   x31 = 1.0/psq + x33;

   x44 = e3/(3.0*p3);
   x43 = 1.5*ecc2/p3;
   x42 = 3.0*ecc/p3 + 2.0*x44;
   x41 = 1.0/p3 + x43;

   x55 = x33*x33;
   x54 = 4.0*x44*x21;
   x53 = 3.0*ecc2/p4 + 1.5*x55;
   x52 = 4.0*ecc/p4 + 2.0*x54;
   x51 = 1.0/p4 + x53;

   x66 = 0.2*e5/p5;
   x65 = 1.25*e4/p5;
   x64 = 10.0*e3/(3.0*p5) + 4.0*x66/3.0;
   x63 = 5.0*ecc2/p5 + 1.5*x65;
   x62 = 5.0*ecc/p5 + 2.0*x64;
   x61 = 1.0/p5 + x63;

   x77 = e6/(6.0*p6);
   x76 = 1.2*e5/p6;
   x75 = 3.75*e4/p6 + 1.25*x77;
   x74 = 20.0*e3/(3.0*p6) + x76/0.75;
   x73 = 7.5*ecc2/p6 + 1.5*x75;
   x72 = 6.0*ecc/p6 + 2.0*x74;
   x71 = 1.0/p6 + x73;
   /* Coefficients for the W's completed */

   double gg1si = 1.0/sqrt(gam1);      /* gam1 > 0 always */
   double gg1psi = gg1si/sqrt(p);      /* p > 0 always */

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
   double bmg   = b - g;
   double bpg   = b + g;
   double d1ma  = 1.0 - a;
   double d1pa  = 1.0 + a;
   double beta1 = bmg/d1ma;
   double beta2 =-bpg/d1pa;
   double b12   = beta1*beta1;
   double b13   = b12*beta1;
   double b22   = beta2*beta2;
   double b23   = b22*beta2;
   xmm1  = sqrt(1.0 - b12*q2);
   
   if (xmm1*alph3 < 0.0)
   {
      xmm1 = -xmm1;
   }

   xmm2 = sqrt(1.0 - b22*q2);

   if (xmm2*alph3 < 0.0)
   {
      xmm2 = -xmm2;
   }

   double d3 = q1;
   double d4 = d1*alph3/(2.*alph2);

   double d1md3 = 1 - d3;
   double bmag  = b - a*g;

   d10 = bmag/bmg*sqrt(d1md3/(d5*(1.0 - xm/b12)));
   d20 = bmag/bpg*sqrt(d1md3/(d5*(1.0 - xm/b22)));

   double dd2 = xm/2.0;
   double dd3 = dd2*g;
   double dd4 = 1.5*dd2*dd2;
   double dd5 = dd4*g;
   double dd6 = dd2*dd4/0.6;

   double c15 = dd6/(b13*b13);
   double c14 = c15 + dd5/(b12*b13);
   double c13 = c14 + dd4/(b12*b12);
   double c12 = c13 + dd3/b13;
   double c11 = c12 + dd2/b12;
   double c10 = c11 + g/beta1;

   double c25 = dd6/(b23*b23);
   double c24 = c25 + dd5/(b22*b23);
   double c23 = c24 + dd4/(b22*b22);
   double c22 = c23 + dd3/b23;
   double c21 = c22 + dd2/b22;
   double c20 = c21 + g/beta2;

   b1q  = beta1*q;
   double b1q2 = b1q*b1q;
   double b1q4 = b1q2*b1q2;
   b2q  = beta2*q;
   double b2q2 = b2q*b2q;
   double b2q4 = b2q2*b2q2;     /* N3 coefficients completed */

   /*
    *  N1 and N2 coefficients
    */
   double xk12 = xk1*xk1;
   double xk13 = xk12*xk1;

   /* Byrd and Friedman formula for elliptic integral */
   double sq  = xk1/16.0 + xk12/32.0 + 21.0*xk13/1024.0;
   double sq2 = sq*sq;
   double sq3 = sq2*sq;

   ucf1 = 2.0*sq/(1.0 + sq2);
   ucf2 = sq2/(1.0 + sq2*sq2);
   ucf3 = 2.0*sq3/(3.0*(1.0 + sq3*sq3));

   double denystt = 1.0 + sq + sq;
   denyst  = denystt*denystt*d1;

   /* N1 coefficients */
   double d1a2 = d1/alph2;
   cn11 = d1a2*asq;
   cn12 = d1a2*2.0*ab*q;
   cn13 = d1a2*(bsq - 4.0*ab*g)*q2;
   cn14 = d1a2*(xm*ab - 2.0*bsq*g)*q2*q;
   cn15 = d1a2*(3.0*bsq*g2 + xm*bsq/2.0)*q4;
   cn16 =-d1a2*xm*bsq*g*q4*q;
   cn17 = 0.375*d1a2*xm*xm*bsq*q4*q2;

   /* N2 coefficients */
   double cn22 = 0.5*xk1*d1;
   double cn24 = 0.375*xk12*d1;
   double cn26 = 0.3125*xk13*d1;

   /* N3 coefficients */
   double d41ma = d4/d1ma;
   double d41pa = d4/d1pa;

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
   double u, d5sq, ut_n, ut_d;
   if (q == 0.0)
   {
      u = 0.0;
   }
   else
   {
      d5sq = sqrt( s1*p0*(1.0 + p1*sigma0 - q1*sigma0*sigma0) );
      ut_n = (sigma0 - a)*d5sq;
      ut_d = sqrg*sqrt(1.0 - g2*q2);

      u = atan2(ut_n, ut_d);
   }

   /*
    * Calculate the Tk = tk, k = 0, T0 = t0 = u.  Here k = 1, ... ,6
    */
   double csu  = cos(u);
   double snu  = sin(u);
   double snu2 = snu*snu;
   double snu4 = snu2*snu2;

   double t1 = 1.0 - csu;
   double t2 = (u - csu*snu)/2.0;
   double t3 = (2.0*t1 - csu*snu2)/3.0;
   double t4 = (3.0*t2 - csu*snu2*snu)/4.0;
   double t5 = (4.0*t3 - csu*snu4)/5.0;
   double t6 = (5.0*t4 - csu*snu4*snu)/6.0;

   /*
    * Compute xhat0 at the initial time which is used only once to get
    * the Jacobi's constants (beta1= - capt, beta2=somega, beta3=comega)
    * Note that xhat0 should not be computed from the initial osculating
    * position and velocity vectors which are Kepler's solutions.
    */
   if (gamma < -1.0e-14)      /* Ellipse */
   {
      //if (fabs(ecc) < 1.0e-10 || fabs(sqrf) < 1.0e-10)
      // Replaced 9/5/97
      if (fabs(ecc) < 1.0e-10 ||
          fabs(sqrf) < 1.0e-10 || fabs(rho0 -rho1) < 1e-5 )
      {
         xhat0 = 0.0;
      }
      else
      {
         //shat = sign(1, sqrf)*sqrt(gamma*rho02 + 2*rho0 - p)/ecc;
         double shat = ( fabs(sqrf)/sqrf )*sqrt(gamma*rho02 + 2*rho0 - p)/ecc;

         double cacs =(rho0 - rho1)*gamma/ecc + 1.0;
         double acs  = atan2(shat*smgam, cacs);
      
         if (acs <  0.0)
         {
            acs = utl_const::tpi + acs;
         }

         xhat0 = acs/smgam;
      }
   }
   else if (gamma > 1.0e-14)                /* Hyperbola */
   {
      if (fabs(ecc) < 1.0e-10 || fabs(sqrf) < 1.0e-10)
      {   
          xhat0 = 0.0;
      }
      else
      {
         double chat = (rho0 - rho1)/ecc;   /* chat > 0, always */
         double zz = gamma*chat + 1.0;
         xhat0 = (log( zz + sqrt(zz*zz - 1.0)))/sqrt(gamma);
         if (sqrf < 0.0)
         {
            xhat0 = -xhat0;
         }
      }
   }
   else                                     /* parabola */
   {
      xhat0 = rrd;                          /* r*rdot */
   }

   /* .............. End initialization ...................*/

   /*
    * Step 5. The last three Jacobi constants, -capt, somega, comega
    *
    * Determine the true anomaly, tra, from xhat0 of ti
    */
   double sneca {};
   double cneca {};
   double tra {};
   double eca = smgam*xhat0;
   if (gamma < -1.0e-14)                    /* Ellipse */
   {
      sneca = sin(eca);
      cneca = cos(eca);
      double s1mes = sqrt(1.0 - ecc2);
      double temp3 = 2.0*atan( (ecc*sneca)/(1 + s1mes - ecc*cneca) );
      tra   = eca + temp3;
   }
   else if (gamma > 1.0e-14)                /* Hyperbola */
   {
      sneca = sinh(eca);
      cneca = cosh(eca);

      double chat = (cneca - 1.0)/gamma;
      double shat = sneca/smgam;
      tra = atan2( sqrt(p)*shat,( rho1 - chat) );
   }
   else                                     /* Parabola */
   {
      tra = 2.0*atan(xhat0/sqrt(p));
   }

   double snw = sin(tra);
   double cnw = cos(tra);
   double ecccnw = ecc*cnw;

   /*
    *  Is there such a thing call parabolic trajectory in the OSC system? 
    */
   //double dwdx  = (1.0 + ecccnw)/sqrt(p); 
   //double dw1dx = (1.0 + ecccnw)/p;

   //double dw2dx = dw1dx*dw1dx;
   //double dw3dx = dw1dx*dw2dx;                // kam
   //double dw4dx = dw2dx*dw2dx;

   double v3 = snw*cnw;
   double v4 = v3*cnw;
   double v5 = v4*cnw;
   double v6 = v5*cnw;
   double v7 = v6*cnw;

   double w1 = x21*tra + x22*snw;
   double w2 = x31*tra + x32*snw + x33*v3;
   double w3 = x41*tra + x42*snw + x43*v3 + x44*v4;
   double w4 = x51*tra + x52*snw + x53*v3 + x54*v4 + x55*v5;
   double w5 = x61*tra + x62*snw + x63*v3 + x64*v4 + x65*v5 + x66*v6;
   double w6 = x71*tra + x72*snw + x73*v3 + x74*v4 + x75*v5 + x76*v6 + x77*v7;

   double uhat {};
   double chat {};
   if (fabs(gamma) >= 1.0e-14)              /* Circle, ellipse, hyperbloa */
   {  
      chat = (cneca - 1.0)/gamma;
      uhat = (sneca - eca)/gams3;
   }
   else                                     /*  Parabola */
   {
      chat = xhat0*xhat0/2.0;
      uhat = xhat0*xhat0*xhat0/6.0;
   }

   double r1 = cr11*xhat0 + cr12*uhat + cr13*tra +
               cr14*w1 + cr15*w2 + cr16*w3 + cr17*w4;
   double r2 = cr21*tra +
               cr22*w1 + cr23*w2 + cr24*w3 + cr25*w4 + cr26*w5 + cr27*w6;

   double sstar = sin(0.5*u);
   double cstar = cos(0.5*u);
   double cb1qs = cstar - b1q*sstar;
   double psi1  = atan(xmm1*sstar/cb1qs);

   if (cb1qs*cos(psi1) < 0.0)
   {
      psi1 = utl_const::pi + psi1;
   }

   double cb2qs = cstar - b2q*sstar;
   double psi2  = atan(xmm2*sstar/cb2qs);

   if (cb2qs*cos(psi2) < 0.0)
   {
      psi2 = utl_const::pi + psi2;
   }

   double r3 = cr31*w2 + cr32*w3 + cr33*w4 + cr34*w5 + cr35*w6;

   double en3 = d10*psi1 + d20*psi2 + cn31*u +
                cn32*t1 + cn33*t2 + cn34*t3 + cn35*t4 + cn36*t5;
   double en1 = cn11*u +
                cn12*t1 + cn13*t2 + cn14*t3 + cn15*t4 + cn16*t5 + cn17*t6;
   double en2 = d1*u + cn22*t2 + cn24*t4 + cn26*t6;
   
   somega  = en2 - r2;                      /* beta2 */
   capt    = t0 - r1 - csq*en1;             /* -beta1 */
   comega  = alph0 + csq*r3 - en3;          /* beta3 */

   oe0[3] = -capt;                     // Vinti mean element  "beta1"
   oe0[4] = somega;                    // Vinti mean element  "beta2"
   oe0[5] = comega;                    // Vinti mean element  "beta3"
}


Eigen::Matrix<double, 6, 1> VintiProp::getStateVector(const JulianDate& jd,
                                                      EphemFrame frame) const
{
  std::array<double, 6> x1;
  vinti_local(jd, x1);

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
  Eigen::Matrix<double, 6, 1> pv = getStateVector(jd, frame);
  Eigen::Matrix<double, 3, 1> pos = pv.block<3,1>(0,0);

  return pos;
}


void VintiProp::vinti_local(const JulianDate& jd,
                            std::array<double, 6>& x1) const
{
   double tf {phy_const::tu_per_day*(jd - jd0)};


     // Return state at epoch for near zero propagation time
   if(fabs(tf - t0) < phy_const::epsdt)
   {
      for (int i = 0; i < 3; i++)
      {
         x1[i]   = pin[i];
         x1[i+3] = vin[i];
      }
      return;
   }

   /*
    * Step 6. At tf, solve the kinematical equations for the oblate
    * spheroidal coordinates (rho, sigma, phi)
    *
    * The generalized Kepler equation is solved by iteration
    * initial guess of xhat is that of the Kepler's solution at tf 
    */
     // universal variable xhat
   double xhat = xhat0 + m_kep->getX(jd);

   double u {};
   double chat {};
   double uhat {};
   double w2 {};
   double w3 {};
   double w4 {};
   double w5 {};
   double w6 {};
   double csu {};
   double snu {};
   double t1 {};
   double t2 {};
   double t3 {};
   double t4 {};
   double t5 {};
   double deltat  = tf - capt;
   for(int ick = 1; ick <= 10; ick++)
   {      
      double eca = smgam*xhat;
      
      double sneca {0.0};                   // kam - < vs. <= below
      double cneca {0.0};
      double tra {};
      if (gamma < -1.0e-14)                 /* Ellipse */
      {
         sneca = sin(eca);
         cneca = cos(eca);
         double s1mes = sqrt(1.0 - ecc2);
         double temp3 = 2.0*atan( (ecc*sneca)/(1 + s1mes - ecc*cneca) );
         tra   = eca + temp3;
      }
      else if (gamma > 1.0e-14)             /* Hyperbola */
      {
         sneca = sinh(eca);
         cneca = cosh(eca);
         chat  = (cneca - 1)/gamma;
         double shat  = sneca/smgam;
         tra   = atan2( sqrt(p)*shat, (rho1 - chat) );
      }
      else                                  /* Parabola */
      {
         tra = 2.0*atan(xhat/sqrt(p));
      }
      
      double snw = sin(tra);
      double cnw = cos(tra);
      double ecccnw = ecc*cnw;
      
      double dwdx  = (1.0 + ecccnw)/sqrt(p);
      double dw1dx = (1.0 + ecccnw)/p;
      double dw2dx = dw1dx*dw1dx;
      double dw3dx = dw1dx*dw2dx;
      double dw4dx = dw2dx*dw2dx;
      
      double v3 = snw*cnw;
      double v4 = v3*cnw;
      double v5 = v4*cnw;
      double v6 = v5*cnw;
      double v7 = v6*cnw;
      
      double w1 = x21*tra + x22*snw;
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
      else                                  /* Parabola */
      {
         chat = xhat*xhat/2.0;
         uhat = xhat*xhat*xhat/6.0;
      }
      
      double r1 = cr11*xhat + cr12*uhat + cr13*tra +
                  cr14*w1 + cr15*w2 + cr16*w3 + cr17*w4;
      double r2 = cr21*tra +
                  cr22*w1 + cr23*w2 + cr24*w3 + cr25*w4 + cr26*w5 + cr27*w6;
      
      //if (iflag == 1) goto LABLE_60; 	// solution converged
      
      double ystar = (r2 + somega)/denyst;
      u = ystar + ucf1*sin(2.0*ystar) +
                  ucf2*sin(4.0*ystar) + ucf3*sin(6.0*ystar);
        
      csu  = cos(u);
      snu  = sin(u);
      double snu2 = snu*snu;
      double snu4 = snu2*snu2;
      
      t1 = 1.0 - csu;
      t2 = (u - csu*snu)/2.0;
      t3 = (2.0*t1 - csu*snu2)/3.0;
      t4 = (3.0*t2 - csu*snu2*snu)/4.0;
      t5 = (4.0*t3 - csu*snu4)/5.0;
      double t6 = (5.0*t4 - csu*snu4*snu)/6.0;
      
      double en1 = cn11*u +
                   cn12*t1 + cn13*t2 + cn14*t3 + cn15*t4 + cn16*t5 + cn17*t6;
      
      double cn1r1  = csq*en1 + r1;
      double psixhat = cn1r1 - deltat;                /* Function */
      
      double dpsx1 = cr13 + cr14*dw1dx + cr15*dw2dx + cr16*dw3dx + cr17*dw4dx;
      double dpsx = cr11 + cr12*chat + dpsx1*dwdx;    /* 1st derivative */
      double delx = psixhat/dpsx;
      
      xhat = xhat - delx;
      
      if (fabs(delx) < 1.0e-14)
      {
        break;
      }
   }
      
   /*
    * Step 7. Final coordinate transformation, after solution converged
    * OSC state vector is (rhof, sigmaf, alphaf, drhof, dsigf, dalphf)
    * ECI state vector is (pf, vf)
    */

   double sstar = sin(0.5*u);
   double cstar = cos(0.5*u);
   double cb1qs = cstar - b1q*sstar;
   double psi1  = atan(xmm1*sstar/cb1qs);

   if (cb1qs*cos(psi1) < 0.0)
   {
      psi1 = utl_const::pi + psi1;
   }

   double cb2qs = cstar - b2q*sstar;
   double psi2  = atan(xmm2*sstar/cb2qs);

   if (cb2qs*cos(psi2) < 0.0)
   {
      psi2 = utl_const::pi + psi2;
   }

   double r3 = cr31*w2 + cr32*w3 + cr33*w4 + cr34*w5 + cr35*w6;
   double en3 = d10*psi1 + d20*psi2 +
                cn31*u + cn32*t1 + cn33*t2 + cn34*t3 + cn35*t4 + cn36*t5;
   double qsnu = q*snu;
   double alphaf = comega - csq*r3 + en3;
   double rhof   = rho1 + ecc*chat;
   double sigmaf = (a + b*qsnu)/(1 + g*qsnu);

     // Final position and velocity
   std::array<double, 3> pf;
   std::array<double, 3> vf;
   pf[2] = rhof*sigmaf;

   double rhof2 = rhof*rhof;
   double sigf2 = sigmaf*sigmaf;

   double rhocsq = rhof2 + csq;
   double df = sqrt(rhocsq*(1 - sigf2));
   double snalp = sin(alphaf);
   double csalp = cos(alphaf);

   pf[0] = df*csalp;
   pf[1] = df*snalp;

   double shat = xhat + gamma*uhat;
   double rcs  = rhof2 + csq*sigf2;

   double rhoqd  = rhof2 - 2.0*a1*rhof + b1;
   double drhofp = sqrt(gam1*rhoqd)/rcs;
   double drhof  = ecc*shat*drhofp;
   double temp12 =(1.0 + p1*sigmaf - q1*sigf2)*(1.0 - g2*q2)*s1;
   double dsigf  = alph2*q*sqrt(temp12)/rcs*csu/(1.0 + g*qsnu);
   double temp20 =(1.0 - sigf2)*(rhof*drhof);
   double temp21 = rhocsq*sigmaf*dsigf;
   double dd =(temp20 - temp21)/df;
   double dalphf = alph3/(df*df);

   vf[0] = dd*csalp - pf[1]*dalphf;
   vf[1] = dd*snalp + pf[0]*dalphf;
   vf[2] = sigmaf*drhof + rhof*dsigf;

     // Propagated state vector
   for(int i = 0; i < 3; i++)
   {
      x1[i] = pf[i];
      x1[i+3] = vf[i];
   }
}


}
