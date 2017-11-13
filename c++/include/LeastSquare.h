#pragma once
#include <Eigen/Core>

#include <bms_api.h>
#include <defs.h>
#include <LinearConstraints.h>
#include <QuadraticObjective.h>
#include <Statistics.h>

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
    class Parameters
    {
    public:
      int    maxIter()       const { return maxIter_; }
      double minNorm_p()     const { return minNorm_p_; }
      double dualEps()       const { return dualEps_; }
      double rankThreshold() const { return rankThreshold_; }

      Parameters& maxIter      (int i)    { maxIter_ = i;       return *this; }
      Parameters& minNorm_p    (double d) { minNorm_p_ = d;     return *this; }
      Parameters& dualEps      (double d) { dualEps_ = d;       return *this; }
      Parameters& rankThreshold(double d) { rankThreshold_ = d; return *this; }

    private:
      int maxIter_ = -1;              //maximum number of active set iterations. If <= 0, defaulted to 10*n
      double minNorm_p_ = 1e-10;      //bound on the L-inf norm of p under which p is treated as 0
      double dualEps_ = 1e-15;        //margin on the 0 for the Lagrange multipliers
      double rankThreshold_ = 1e-12;  //threshold on the last diagonal value in the QR of A, to detect rank loss
    };


    LeastSquare(int n);

    void parameters(const Parameters& param);
    const Parameters& parameters() const;

    SolverStatus solve(const LeastSquareObjective& obj, const VectorConstRef& Jx0, const VectorConstRef& j, double c, LinearConstraints& lc);
    SolverStatus solveFeasibility(const VectorConstRef& j, double c, LinearConstraints& lc);

    /** Retrieve the solution (after call to solve() or solveFeasibility())*/
    const Eigen::VectorXd& x() const;
    /** Retrieve the Lagrange multipliers (after call to solve() or solveFeasibility())*/
    const Eigen::VectorXd& lambda() const;

    /** Get stats on the last run. Only meaningful if USE_STATS is defined*/
    const stats::LSStats& statistics() const;

  private:
    Eigen::DenseIndex n_;

    //solver parameters;
    Parameters params_;

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

    //stats
    stats::LSStats stats_;
  };
}