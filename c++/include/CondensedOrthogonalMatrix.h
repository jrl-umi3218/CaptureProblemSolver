#pragma once

#include <vector>

#include <bms_api.h>
#include <GivensSequence.h>

namespace bms
{
  /** A class representing a product Q1 Q2 ... Qk P where the QI are orthogonal
    * matrices and P is a permutation matrix.
    * The matrices Qi are expressed as a product of Givens rotations G1 ... Gp.
    *
    * The goal of this class is twofold:
    * - providing an abstraction for easier manipulation
    * - preallocating memory to avoid memory allocation in critical code
    */
  class BMS_DLLAPI CondensedOrthogonalMatrix
  {
  public:
    /** Create an instance preallocating kmax Givens sequences with pmax Givens
      * rotation each. The represented matrix is nxn*/
    CondensedOrthogonalMatrix(int n, int kmax, int pmax);

    void reset();

    GivensSequence& Q(size_t i);
    Eigen::Transpositions<Eigen::Dynamic>& P();

    /** M = P^T Q_k^T Q_{k-1}^T .... Q_1^T M*/
    template<typename Derived>
    void applyTo(const Eigen::MatrixBase<Derived>& M) const;

    /** Computes M = M * Q_1 Q_2 ... G_n P */
    template<typename Derived>
    void applyOnTheRightTo(const Eigen::MatrixBase<Derived>& M) const;

    /** Return the corresponding nxn orthogonal matrix
    *
    * Use only for debugging purposes.
    */
    Eigen::MatrixXd matrix();

  private:
    Eigen::DenseIndex n_;
    std::vector<GivensSequence> sequences_;
    Eigen::Transpositions<Eigen::Dynamic> transpositions_;
  };

  template<typename Derived>
  inline void CondensedOrthogonalMatrix::applyTo(const Eigen::MatrixBase<Derived>& M) const
  {
    for (const auto& Q : sequences_)
      Q.applyTo(M);
    const_cast<MatrixBase<Derived>&>(M) = transpositions_.transpose()*M;
  }

  template<typename Derived>
  inline void CondensedOrthogonalMatrix::applyOnTheRightTo(const Eigen::MatrixBase<Derived>& M) const
  {
    for (const auto& Q : sequences_)
      Q.applyOnTheRightTo(M);
    const_cast<MatrixBase<Derived>&>(M) = M*transpositions_;
  }

  inline GivensSequence& CondensedOrthogonalMatrix::Q(size_t i)
  {
    return sequences_[i];
  }

  inline Eigen::Transpositions<Eigen::Dynamic>& CondensedOrthogonalMatrix::P()
  {
    return transpositions_;
  }
}