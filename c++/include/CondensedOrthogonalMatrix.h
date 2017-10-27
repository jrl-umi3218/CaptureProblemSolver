#pragma once

#include <vector>

#include <bms_api.h>
#include <GivensSequence.h>

#include <iostream>

namespace bms
{
  /** A class representing a product Q1 Q2 ... Qk P Qh where the Qi are 
    * orthogonal matrices and P is a permutation matrix.
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
      * rotation each. The represented matrix is nxn
      * If Ptranspose is true, P^T is stored instead of P.
      */
    CondensedOrthogonalMatrix(int n, int kmax, int pmax, bool Ptranspose = false);

    /**Empty version*/
    CondensedOrthogonalMatrix() {}

    /** This methods acts exactly as the operator=(), except that it supposes 
      * (and check by assert) that other will fit in the memory already allocated
      * to this object.
      */
    CondensedOrthogonalMatrix& copyNoAlloc(const CondensedOrthogonalMatrix& other);

    Eigen::DenseIndex size() const;
    void reset(bool Ptranspose = false);
    void resize(int n, int kmax, int pmax);

    std::vector<GivensSequence>& Q();
    GivensSequence& Q(size_t i);
    GivensSequence& Qh();
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
    bool ptranspose_;
    Eigen::DenseIndex n_;
    std::vector<GivensSequence> sequences_;
    GivensSequence Qh_;
    Eigen::Transpositions<Eigen::Dynamic> transpositions_;
    Eigen::VectorXi identityTransposition_;
  };

  inline CondensedOrthogonalMatrix& CondensedOrthogonalMatrix::copyNoAlloc(const CondensedOrthogonalMatrix& other)
  {
    assert(other.sequences_.size() <= sequences_.size());
    assert(other.n_ == n_);
    ptranspose_ = other.ptranspose_;
    Qh_ = other.Qh_;
    transpositions_ = other.transpositions_;
    //we want to keep the preallocated memory, so we just clear the GivensSequences we don't use
    size_t i;
    for (i = 0; i < other.sequences_.size(); ++i)
    {
      assert(other.sequences_[i].size() <= sequences_[i].capacity());
      sequences_[i] = other.sequences_[i];
    }
    size_t kmax = sequences_.size();
    for (; i < kmax; ++i)
      sequences_[i].clear();

    return *this;
  }

  inline Eigen::DenseIndex CondensedOrthogonalMatrix::size() const
  {
    return n_;
  }

  template<typename Derived>
  inline void CondensedOrthogonalMatrix::applyTo(const Eigen::MatrixBase<Derived>& M) const
  {
    for (const auto& Q : sequences_)
      Q.applyTo(M);
    if (ptranspose_)
      const_cast<MatrixBase<Derived>&>(M) = transpositions_*M;
    else
      const_cast<MatrixBase<Derived>&>(M) = transpositions_.transpose()*M;
    Qh_.applyTo(M);
  }

  template<typename Derived>
  inline void CondensedOrthogonalMatrix::applyOnTheRightTo(const Eigen::MatrixBase<Derived>& M) const
  {
    for (const auto& Q : sequences_)
      Q.applyOnTheRightTo(M);
    if (ptranspose_)
      const_cast<MatrixBase<Derived>&>(M) = M*transpositions_.transpose();
    else
      const_cast<MatrixBase<Derived>&>(M) = M*transpositions_;
    Qh_.applyOnTheRightTo(M);
  }

  inline std::vector<GivensSequence>& CondensedOrthogonalMatrix::Q()
  {
    return sequences_;
  }

  inline GivensSequence& CondensedOrthogonalMatrix::Q(size_t i)
  {
    return sequences_[i];
  }

  inline GivensSequence & CondensedOrthogonalMatrix::Qh()
  {
    return Qh_;
  }

  inline Eigen::Transpositions<Eigen::Dynamic>& CondensedOrthogonalMatrix::P()
  {
    return transpositions_;
  }
}