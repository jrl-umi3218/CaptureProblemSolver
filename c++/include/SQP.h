#pragma once

#include <Eigen/Core>

#include <bms_api.h>
#include <LeastSquare.h>
#include <Problem.h>

namespace bms
{
  class BMS_DLLAPI SQP
  {
  public:
    SQP(int n);

    SolverStatus solveFeasibility(const Problem& pb);

    const Eigen::VectorXd& x() const;
    const Eigen::VectorXd& lambda() const;

  private:
    bool checkKKTFeasibility(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, double f, const Eigen::VectorXd& g, const LinearConstraints& lc) const;

    Eigen::DenseIndex n_;
    LeastSquare ls_;

    //optimization parameters
    int maxIter_ = 100;
    double beta = 0.9;            //backtracking multiplier in line search
    double c1 = 0.01;             //gradient coefficient in line search
    double smallestLSStep = 1e-8; 
    double tau_p = 1e-6;          //precision parameter on primal condition
    double tau_d = 1e-6;          //precision parameter on dual condition

    //optimization data
    LinearConstraints shiftedLC_;
    Eigen::VectorXd x_;
    Eigen::VectorXd xa_;
    Eigen::VectorXd j_;
    Eigen::VectorXd lambda_;
    mutable Eigen::VectorXd Cx_;
    mutable Eigen::VectorXd Cl_;
    mutable Eigen::VectorXd gradL_;
  };
}