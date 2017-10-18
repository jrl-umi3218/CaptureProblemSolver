#pragma once

#include <vector>

#include <Eigen/Core>

#include <bms_api.h>
#include <defs.h>

namespace bms
{
  enum class Activation
  {
    None,
    Lower,
    Upper,
    Equal
  };


  /** A class to manipulate the linear constraints of the problem.
    * l_i <= x_i-x_{i-1} <= u_i for i=0..n-1  (x_{-1} = 0)
    * xln <= x_{n-1} <= xun
    * Constraints are considered in this order
    */
  class BMS_DLLAPI LinearConstraints
  {
  private:
    /** Intermediate temporary structure for shifting the constraints. See shift*/
    struct Shift
    {
      LinearConstraints& c;
      const Eigen::VectorXd& x;
    };

  public:
    LinearConstraints(const Eigen::VectorXd& l, const Eigen::VectorXd& u, double xln, double xun);

    /** Build the constraints on p such that on s.x+p verifies the original constraints s.c*/
    LinearConstraints& operator&= (const Shift& s);

    Shift shift(const Eigen::VectorXd& x) const;

    /** Activate the i-th constraint
      *
      * If a = Activation::None, this effectively deactivate the i-th constraint
      */
    void activate(size_t i, Activation a);
    /** Deactivate the i-th constraint*/
    void deactivate(size_t i);

    Activation activationStatus(size_t i) const;

    const std::vector<bool>& activeSet() const;
    Eigen::DenseIndex numberOfActiveConstraints() const;

    /** Y = X*N_A */
    void applyNullSpaceOnTheRight(MatrixRef Y, const MatrixConstRef& X) const;
    /** Y = N_A X */
    void applyNullSpaceOnTheLeft(MatrixRef Y, const MatrixConstRef& X) const;

    /** Y = C_A * X */
    void mult(MatrixRef Y, const MatrixConstRef& X) const;
    /** Y = C_A^T * X */
    void transposeMult(MatrixRef Y, const MatrixConstRef& X) const;
    /** Y = pinv(C_A^T)*X */
    void pinvTransposeMult(MatrixRef Y, const MatrixConstRef& X);

    /** Move the furthest possible on the segment x+ap for 0<=a<=1. Stop at the
      * first constraint encountered and activate it. Returns x+ap.
      */
    Eigen::VectorXd performQPstep(const Eigen::VectorXd& x, const Eigen::VectorXd& p);

    /** The matrix of active constraints. For debugging purposes*/
    Eigen::MatrixXd matrix() const;

  private:
    void computeIdx() const;

    //data
    Eigen::DenseIndex n_;
    Eigen::VectorXd l_;
    Eigen::VectorXd u_;
    double xln_;
    double xun_;

    //active-set management
    Eigen::DenseIndex na_; //number of activated constraint
    std::vector<Activation> activationStatus_; //type of activation of the constraints
    std::vector<bool> activeSet_; //Invariant: activeSet_[i] == (activationStatus_[i] != Activation::None) for all i

    //internal intermediate data
    mutable bool validIdx_;                       // indicate if idx_ is valid
    mutable std::vector<Eigen::DenseIndex> idx_;  // indices used for the premultiplication by N_A
  };
}