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

  /** An information returned by the search of a feasible point
    *  - Found: the point was found
    *  - Infeasible: the constraints are not compatible
    *  - NumericalWarningXXX: the lower bound on x_n, its upper bound or both
    *    are numerically redundant with the other constraints. Please consider
    *    moving them.
    *  - TooSmallZonotope: the n first constraints define a zonotope too small.
    *  - TooSmallFeasibilityZone: the feasible zone is too small.
    */
  enum class FeasiblePointInfo
  {
    Found,
    Infeasible,
    NumericalWarningLower,
    NumericalWarningUpper,
    NumericalWarningBoth,
    TooSmallZonotope,
    TooSmallFeasibilityZone
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
      const LinearConstraints& c;
      const Eigen::VectorXd& x;
      bool feasible;
    };

  public:
    LinearConstraints(int n);
    LinearConstraints(const Eigen::VectorXd& l, const Eigen::VectorXd& u, double xln, double xun);

    /** See shift*/
    LinearConstraints& operator= (const Shift& s);

    /** Build the constraints on p such that on s.x+p verifies the original constraints s.c.
      * If feasible is true, try to correct numerical errors by making:
      * - l<=0
      * - u>=0
      * - l[i] = u[i] = 0 for equality constraints
      * It is the user responsibility to ensure that x is indeed feasible in this case.
      */
    Shift shift(const Eigen::VectorXd& x, bool feasible=false) const;

    const Eigen::VectorXd& l() const;
    const Eigen::VectorXd& u() const;

    /** Return a feasible point for the constraints*/
    std::pair<FeasiblePointInfo, Eigen::VectorXd> initialPoint() const;

    Eigen::DenseIndex size() const;
    Eigen::DenseIndex nullSpaceSize() const;

    /** Activate the i-th constraint
      *
      * If a = Activation::None, this effectively deactivate the i-th constraint
      */
    void activate(size_t i, Activation a);
    /** Deactivate the i-th constraint*/
    void deactivate(size_t i);
    /** Deactivate all non-equality constraints*/
    void resetActivation();

    Activation activationStatus(size_t i) const;

    const std::vector<bool>& activeSet() const;
    const std::vector<Activation>& activationStatus() const;
    void setActivationStatus(const std::vector<Activation>& act);
    /** Index of active constraints*/
    const std::vector<Eigen::DenseIndex>& activeSetIdx() const;
    Eigen::DenseIndex numberOfActiveConstraints() const;

    void changeBounds(Eigen::DenseIndex i, double l, double u);

    /** Y = X*N_A */
    void applyNullSpaceOnTheRight(MatrixRef Y, const MatrixConstRef& X) const;
    /** Y = N_A X */
    void applyNullSpaceOnTheLeft(MatrixRef Y, const MatrixConstRef& X) const;

    /** Y = C * X */
    void mult(MatrixRef Y, const MatrixConstRef& X) const;
    /** Y = C^T * X */
    void transposeMult(MatrixRef Y, const MatrixConstRef& X) const;
    /** Y = pinv(C_A^T)*X */
    void pinvTransposeMultAct(MatrixRef Y, const MatrixConstRef& X) const;

    /** fill y such that y(!active) = 0 and y(active) = x*/
    void expandActive(VectorRef y, const VectorConstRef& x) const;

    /** Move the furthest possible on the segment x+ap for 0<=a<=1. Stop at the
      * first constraint encountered and activate it. Returns x+ap.
      *
      * If fullstep is passed by the user, set it to true when a full step is
      * performed, false otherwise
      */
    Eigen::VectorXd performQPstep(const Eigen::VectorXd& x, const Eigen::VectorXd& p, bool* fullstep = nullptr);

    /** Deactivate the constraint with largest lambda violation. */
    void deactivateMaxLambda(const VectorConstRef& lambda);

    /** The matrix of all constraints. For debugging purposes.*/
    Eigen::MatrixXd matrix() const;
    /** The matrix C_A of active constraints. For debugging purposes.*/
    Eigen::MatrixXd matrixAct() const;

    /** Check if x verifies the constraints*/
    bool checkPrimal(const VectorConstRef& x, double eps = 1e-15) const;
    /** Check the signs of lambda*/
    bool checkDual(const VectorConstRef& lambda, double eps = 1e-15) const;

  private:
    void computeIdx() const;
    void computeActIdx() const;

    //data
    Eigen::DenseIndex n_;
    Eigen::VectorXd l_;   //[l;xln]
    Eigen::VectorXd u_;   //[u;xun]

    //active-set management
    Eigen::DenseIndex na_; //number of activated constraint
    std::vector<Activation> activationStatus_; //type of activation of the constraints
    std::vector<bool> activeSet_; //Invariant: activeSet_[i] == (activationStatus_[i] != Activation::None) for all i

    //internal intermediate data
    mutable bool validIdx_;                         // indicate if idx_ is valid
    mutable bool validActIdx_;                      // indicate if idx_ is valid
    mutable std::vector<Eigen::DenseIndex> idx_;    // indices used for the premultiplication by N_A
    mutable std::vector<Eigen::DenseIndex> actIdx_; // indices of active constraints, ordered.
    mutable Eigen::VectorXd Cx_;                    // to store the results of C*x
    mutable Eigen::VectorXd Cp_;                    // to store the results of C*p
  };
}