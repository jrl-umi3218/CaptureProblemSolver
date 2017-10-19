#pragma once
#include <Eigen/Core>

#include <bms_api.h>
#include <defs.h>
#include <LinearConstraints.h>

namespace bms
{
  enum class LSStatus
  {
    Converge,
    MaxIteration,
    Fail
  };

  class BMS_DLLAPI LeastSquare
  {
  public:
    LeastSquare(int n);

    LSStatus solveFeasibility(const VectorConstRef& j, double c, LinearConstraints& lc);

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
  };
}