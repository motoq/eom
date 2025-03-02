#ifndef UTL_CPP_TMP_H
#define UTL_CPP_TMP_H

#include <limits>   

/*
 * Local functions that may not make sense outside of this namespace
 */
namespace {

/*
 * Applies a differential correction via Newton's method to determine
 * r given s where s = r*r.  Let f(r) = r*r - s = 0.  Then, f' = 2*r.
 * Therefore:
 *
 *   r_new = r - (r*r - s)/(2r)
 *         = 2r*r/(2r) - r*r/(2r) + s/(2*r)
 *         = r/2 + s/(2*r)
 *         = 0.5*(r + s/r)
 *
 * @param  x      Squared value (s from above) - fixed value across
 *                recursive calling sequence
 * @param  curr   Current estimate (r)
 * @param  previ  Previous estimate used to determine convergence
 *
 * @return  r_new
 *
 * From:  Alex Shtoff  2015/12/07
 *        https://stackoverflow.com/questions/8622256/
 *                in-c11-is-sqrt-defined-as-constexpr
 */
double constexpr nr_sqrt(double x, double curr, double prev)
{
    // Exit when there is no change in value (to numeric precision)
  return curr == prev ? curr
                      : nr_sqrt(x, 0.5*(curr + x/curr), curr);
}

}

namespace utl_cpp_tmp {

/**
 * Constexpr version of the square root function calling Newton's
 * method via recursion.  C++26 begins inclusion of a constexpr
 * sqrt() - until then, this implementation can be used.
 *
 * @param  x  Value for which to take the square root of
 *
 * @param  Square root of x if 0 <= x < inf, otherwise, NaN
 *
 * From:  Alex Shtoff  2015/12/07
 *        https://stackoverflow.com/questions/8622256/
 *                in-c11-is-sqrt-defined-as-constexpr
 */
double constexpr constexpr_sqrt(double x)
{
    // Perform recursion if x is positive and less than infinity
    // Use x as initial guess to sqrt(x)
  return x >= 0.0  &&  x < std::numeric_limits<double>::infinity()
      ? nr_sqrt(x, x, 0.0)
      : std::numeric_limits<double>::quiet_NaN();
}


}

#endif
