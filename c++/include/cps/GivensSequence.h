#pragma once

#include <vector>

#include <cps/cps_api.h>
#include <cps/Givens.h>

namespace cps
{
  /** A sequence of Givens rotations G_0*G_1*...
    *
    * Note: it is ok to derive publicly here as we do not have any additional
    * data.
    */
  class CPS_DLLAPI GivensSequence final : public std::vector<Givens>
  {
  public:
    using std::vector<Givens>::vector;

    /** M = G_{n-1}^T G_{n-2}^T .... G_0^T M*/
    template<typename Derived>
    void applyTo(const Eigen::MatrixBase<Derived>& M) const;

    /** Computes M = M * G_0 G_1 ... G_{n-1}*/
    template<typename Derived>
    void applyOnTheRightTo(const Eigen::MatrixBase<Derived>& M) const;

    /** Perform G.extend(incr) on all Givens rotations in the sequence*/
    void extend(int incr);

    /** Return the corresponding nxn orthogonal matrix
      *
      * Use only for debugging purposes.
      */
    Eigen::MatrixXd matrix(int n);
  };


  template<typename Derived>
  inline void GivensSequence::applyTo(const Eigen::MatrixBase<Derived>& M) const
  {
    for (const auto& G : *this)
      G.applyTo(M);
  }

  template<typename Derived>
  inline void GivensSequence::applyOnTheRightTo(const Eigen::MatrixBase<Derived>& M) const
  {
    for (const auto& G : *this)
      G.applyOnTheRightTo(M);
  }
}
