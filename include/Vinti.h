/*----------------------------------------------------------------------------
 * VINTI.H ... Header for Vinti and related functions
 *--------------------------------------------------------------------------*/

#if !defined(_VINTI_H_)
#define _VINTI_H_

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif 

double hmod360(double);
void VintToKep(double[4], double[6], double[6]);
void Kepler1(const double[4], double, const double[6],
                              double, double[6], double *);
void Vinti6 (const double[4], double, const double[6],
                              double, double[6], double[6]);

#ifdef __cplusplus
}
#endif

#endif
