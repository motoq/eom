/******************************************************************************
 *  Orbital and Celestial Mechanics		                                       *
 *  ------------------------------------------------------------------------  *
 *                                                                            *
 *  Kepler1()...Kepler Solution with Universal Variables.                     *
 *                                                                            *
 *  Originator:                 Gim J. Der, 1995                              *
 *  Fortran to C Conversion:    Herbert B. Reynolds, 1997                     *
 *  ------------------------------------------------------------------------  *
 *                                                                            *
 *  This function computes the kepler solution (state vector, x1, at          *
 *  any time, t1) from a given initial state vector, x0, at an arbitrary      *
 *  time t0, using the universal variable formulation.                        *
 *                                                                            *
 *  It accommodates all orbits (circular, elliptic, parabolic,                *
 *  hyperbolic) with no discontinuity or singularity.  No problem of          *
 *  convergence was encountered, and it always converges in less than:        *
 *                                                                            *
 *      a. five iterations for circle and ellipse                             *
 *      b. seven iterations for parabola and hyperbola                        *
 *                                                                            *
 *  The "if(dx2 > 0)" loop is not for efficiency but to guarantee a good      *
 *  solution for parabolic and hyperbolic orbits up to 1,000,000 km from      *
 *  Earth.                                                                    *
 *                                                                            *
 *  Input variables:                                                          *
 *          x0[6] = initial state vector at t0, (ECI),     (km,km/s)          *
 *          t0    = time with respect to x0, initial time, (sec)              *
 *          t1    = time with respect to x1, final time,   (sec)              *
 *                                                                            *
 *          Note: unit (km, km/s) depends only on that of the data gmx, rex   *
 *                                                                            *
 *  Output variables:                                                         *
 *          x1[6] = state vector at t1:(t0+dt),   (km,km/s)                   *
 *          xxx   = universal variable, xhat in Getchell at t1  (Unitless)    *
 *                                                                            *
 *****************************************************************************/

#include <astro_kepler_prop.h>

#include <string>
#include <array>
#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <phy_const.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>

#include <Vinti.h>

namespace {
  constexpr int kn {10};
  constexpr double muqr {1.0};
}

namespace eom {


KeplerProp::KeplerProp(const std::string& orbit_name,
                       const JulianDate& epoch,
                       const Eigen::Matrix<double, 6, 1>& xeci,
                       const std::shared_ptr<const EcfEciSys>& ecfeciSys)
{
   name = orbit_name;
     // Initial state
   jd0 = epoch;
   ecfeci = ecfeciSys;
   for (int ii=0; ii<3; ++ii)
   {
     r0[ii] = xeci(ii);
     v0[ii] = xeci(ii+3);
   }
     // Fixed parameters based on initial state
   r0mag = sqrt( r0[0]*r0[0] + r0[1]*r0[1] + r0[2]*r0[2] );
   v0mag = sqrt( v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2] );
   d0 = r0[0]*v0[0] +r0[1]*v0[1] +r0[2]*v0[2];
}


  //kepler1_mod(0.0, x0.data(), t1, x1.data(), &x);
Eigen::Matrix<double, 6, 1> KeplerProp::getStateVector(const JulianDate& jd,
                                                       EphemFrame frame) const 
{
     // Time since initial state
   double dt {phy_const::tu_per_day*(jd - jd0)};
     // Propagated state vector
   std::array<double, 6> x1;
     // Universal variable for state transition
   double x {0.0};

   if(std::fabs(dt) < phy_const::epsdt)
   {
      for(int i = 0; i < 3; i++)
      {
          x1[i]   = r0[i];
          x1[i+3] = v0[i];
      }
   } else {
      double dfx, u1, u2, u3;
      x = f_and_g(dt, dfx, u1, u2, u3);
      /*
       *  Kepler solution  converged
       *  Coefficients and constants associated with the partials
       */
      double rmag {dfx};
      double f {1.0 - u2/r0mag};
      double g {dt - u3/muqr};

      /*
       *  Time derivatives associated with the partials
       */
      double df {-muqr*u1/(rmag*r0mag)};
      double dg {1.0 - u2/rmag};

      for(int i = 0; i < 3; i++)
      {
         x1[i] = ( f*r0[i] + g*v0[i] );
         x1[i+3] = ( df*r0[i] + dg*v0[i] );
      }
   }

   Eigen::Matrix<double, 6, 1> xeci;
   xeci(0) = x1[0];
   xeci(1) = x1[1];
   xeci(2) = x1[2];
   xeci(3) = x1[3];
   xeci(4) = x1[4];
   xeci(5) = x1[5];

   if (frame == EphemFrame::ecf)
   {
      return ecfeci->eci2ecf(jd, xeci.block<3,1>(0,0), xeci.block<3,1>(3,0));
   }
   return xeci;
}


Eigen::Matrix<double, 3, 1> KeplerProp::getPosition(const JulianDate& jd,
                                                    EphemFrame frame) const
{
   Eigen::Matrix<double, 6, 1> xeci = getStateVector(jd, frame);
   Eigen::Matrix<double, 3, 1> posi = xeci.block<3,1>(0,0);

   if (frame == EphemFrame::ecf)
   {
      return ecfeci->eci2ecf(jd, posi);
   }
   return posi;
}


double KeplerProp::getX(const JulianDate& jd) const
{
   double dfx {};
   double u1 {};
   double u2 {};
   double u3 {};
     // Time since initial state
   double dt {phy_const::tu_per_day*(jd - jd0)};

   return f_and_g(dt, dfx, u1, u2, u3);
}


double KeplerProp::f_and_g(double dt,
                           double& dfx,
                           double& u1, double& u2, double& u3) const
{
     // Nothing to do at epoch
   if(std::fabs(dt) < phy_const::epsdt)
   {
     dfx = 0.0;
     u1 = 0.0;
     u2 = 0.0;
     u3 = 0.0;
     return 0.0;
   }

   double sigma0 {d0 / muqr};
   double alp0 {2.0 / r0mag - v0mag*v0mag};

   /*
    * x = Flight angle (theta) computed by kepler's solution
    *     first guessed of:
    *
    *     x = alp0*muqr*(t - t0)          For circle and ellipse
    *     x = muqr*(t - t0)/(10*r0mag)    For parabola and hyperbola
    */
   
   double x {alp0*dt};

   if (alp0 <= 0)
   {
      //x = 0.1*dt/r0mag;
      x = 0.5*dt/r0mag; // replace 9/5/97
   }

   double y, yqr, cy, sy;
   double fx, dfx2, sdfx;
   double dx2;
   double dx;
   for(int k = 1; k <= kn; k++)
   {
      if(alp0 < 0)
      {
         y = alp0*x*x;
         yqr = sqrt(-y);
         cy = (1 - cosh(yqr))/y;
         sy = (sinh(yqr) - yqr)/(yqr*yqr*yqr);
      }
      else if(alp0 == 0)
      {
         y = 0;
         cy = 0.5;
         sy = 1.0/6.0;
      }
      else if(alp0 > 0)
      {
         y = alp0*x*x;
         yqr = sqrt(y);
         cy = ( 1 - cos(yqr) )/y;
         sy = (yqr - sin(yqr))/(yqr*yqr*yqr);
      }

      u1 = x*(1 - y*sy);
      u2 = x*x*cy;
      u3 = x*x*x*sy;

      fx = r0mag*u1 + sigma0*u2 + u3 - dt*muqr;
      dfx = sigma0*u1 + (1 - alp0*r0mag)*u2 + r0mag;
      dfx2 = sigma0*(1 - y*cy) + (1 - alp0*r0mag)*u1;
      sdfx = dfx /(fabs(dfx));
      
      dx2 = 16*dfx*dfx - 20*fx*dfx2;

      if(dx2 > 0)
      {
         dx = 5*fx/( dfx + sdfx*sqrt(dx2) );
      }
      else
      {
         dx = 0.5*x;
      }
      
      if(fabs(dx) < 1.0e-10) break;

      x = x - dx;
   }   

   return x;
}


}


