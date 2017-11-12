#pragma once
#include <Eigen/Core>

#include <bms_api.h>
#include <defs.h>
#include <LinearConstraints.h>
#include <QuadraticObjective.h>

namespace bms
{
  enum class SolverStatus
  {
    Converge,                       /** A solution was found*/
    MaxIteration,                   /** The maximum number of iteration was reached*/
    LineSearchFailed,               /** Step in the line search became too small*/
    NoLinearlyFeasiblePoint,        /** Linear inequality constraints are incompatible*/
    NumericallyEquivalentIterates,  /** Previous and new iterates are equal*/
    Fail
  };

  class BMS_DLLAPI LeastSquare
  {
  public:
    LeastSquare(int n);

    SolverStatus solve(const LeastSquareObjective& obj, const VectorConstRef& Jx0, const VectorConstRef& j, double c, LinearConstraints& lc);
    SolverStatus solveFeasibility(const VectorConstRef& j, double c, LinearConstraints& lc);

    /** Retrieve the solution (after call to solve() or solveFeasibility())*/
    const Eigen::VectorXd& x() const;
    /** Retrieve the Lagrange multipliers (after call to solve() or solveFeasibility())*/
    const Eigen::VectorXd& lambda() const;

  private:
    Eigen::DenseIndex n_;

    //solver parameters;
    int maxIter_;

    //computation data
    Eigen::VectorXd x_;
    Eigen::VectorXd p_;
    Eigen::VectorXd z_;
    Eigen::VectorXd lambda_;
    Eigen::VectorXd lambdaAct_;
    Eigen::VectorXd jN_;
    Eigen::MatrixXd A_;
    Eigen::VectorXd b_;
    Eigen::VectorXd Jx_;
    Eigen::VectorXd JtJx_;
    Eigen::VectorXd tmp_;
    GivensSequence Qg_;
    CondensedOrthogonalMatrix Q_;
  };
}