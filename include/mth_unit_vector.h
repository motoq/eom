/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_UNIT_VECTOR_H
#define MTH_UNIT_VECTOR_H

#include <Eigen/Dense>

namespace eom {

/**
 * Class and functionality to generate unit vectors and their first two
 * derivatives.  Note, the derivatives of the unit vectors are
 * generated, not the trivial case of normalizing the derivatives of the
 * vectors.
 *
 * The vector (and its unit vector) are orthogonal to the first
 * derivative of the unit vector.
 *
 * @tparam  T  Vector component data type
 * @tparam  N  Dimension of input vectors
 *
 * @author  Kurt Motekew
 * @date    2024/12/16
 */
template<typename T, int N>
class UnitVector {
public:
  /**
   * Initialize position and time derivative of position
   *
   * @param  r      Position vector
   * @param  rdot   Rate of change (velocity) of position vector
   * @param  rddot  Second Rate of change (acceleration) of position vector
   *                Optional, defaults to zero
   */
  UnitVector(const Eigen::Matrix<T, N, 1>& r,
             const Eigen::Matrix<T, N, 1>& rdot,
             const Eigen::Matrix<T, N, 1>& rddot =
                 Eigen::Matrix<T, N, 1>::Zero());

  /**
   * @return  Unit vector
   */
  Eigen::Matrix<T, N, 1> getNormalized() const
  {
    return m_rhat;
  }

  /**
   * @return  First derivative of unit vector
   */
  Eigen::Matrix<T, N, 1> getNormalizedDot() const
  {
    return m_rhatdot;
  }

  /**
   * @return  Second derivative of unit vector
   */
  Eigen::Matrix<T, N, 1> getNormalizedDDot() const
  {
    return m_rhatddot;
  }

private:
  Eigen::Matrix<T, N, 1> m_rhat    {Eigen::Matrix<T, N, 1>::Zero()};
  Eigen::Matrix<T, N, 1> m_rhatdot {Eigen::Matrix<T, N, 1>::Zero()};
  Eigen::Matrix<T, N, 1> m_rhatddot {Eigen::Matrix<T, N, 1>::Zero()};
};


/**
 * Given a vector and its derivative, return the first derivative of the
 * associated unit vector unit vector.  This function is more convenient
 * and efficient when only the first derivative of the unit vector is
 * needed.
 *
 * @param  r      Position vector
 * @param  rdot   Rate of change (velocity) of position vector
 * 
 * @return  First derivative of unit vector
 */
template<typename T, int N>
Eigen::Matrix<T, N, 1> unit_vector_dot(const Eigen::Matrix<T, N, 1>& r,
                                       const Eigen::Matrix<T, N, 1>& rdot)
{
  T inv_rmag {1.0/r.norm()};
  Eigen::Matrix<T, N, 1> rhat = inv_rmag*r;
    // Projection onto normal of r
  Eigen::Matrix<T, N, N> proj = Eigen::Matrix<T, N, N>::Identity() - 
                                rhat*rhat.transpose();

  return inv_rmag*proj*rdot;
}

/*
 * Constructor generating rhat, rhat_dot, and rhat_ddot
 */
template<typename T, int N>
UnitVector<T,N>::UnitVector(const Eigen::Matrix<T, N, 1>& r,
                            const Eigen::Matrix<T, N, 1>& rdot,
                            const Eigen::Matrix<T, N, 1>& rddot)
{
  T inv_rmag {1.0/r.norm()};
  Eigen::Matrix<T, N, 1> rhat = inv_rmag*r;
    // Projection onto normal of r
  Eigen::Matrix<T, N, N> rhrht = rhat*rhat.transpose();
  Eigen::Matrix<T, N, N> proj = Eigen::Matrix<T, N, N>::Identity() - rhrht;

  m_rhat = rhat;
  m_rhatdot = inv_rmag*proj*rdot;
  m_rhatddot = inv_rmag*(proj*rddot + inv_rmag*(
                        rhat.dot(rdot)*(3.0*rhrht - 
                          Eigen::Matrix<T, N, N>::Identity())*rdot -
                          (rhat*rdot.transpose() + rdot*rhat.transpose())*rdot
                                               )
                        );
}


}

#endif
