#pragma once

#include <Eigen/Core>

#include <bms_api.h>
#include <LeastSquare.h>
#include <Problem.h>
#include <QuadraticObjective.h>

namespace bms
{
  class BMS_DLLAPI SQP
  {
  public:
    SQP(int n);

    SolverStatus solve(const Problem& pb);
    SolverStatus solveFeasibility(const Problem& pb);

    const Eigen::VectorXd& x() const;
    const Eigen::VectorXd& lambda() const;
    const std::vector<Activation>& activeSet() const;
    int numberOfIterations() const;

  private:
    bool checkKKT(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, double f, const Eigen::VectorXd& g, 
                  const LinearConstraints& lc, const LeastSquareObjective* const obj = nullptr) const;

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
    LinearConstraints shiftedLC_;                 //constraints passed to the LS
    Eigen::VectorXd x_;                           //value of the iterate
    Eigen::VectorXd xa_;                          //value of the current point in the line search
    Eigen::VectorXd j_;                           //gradient of the non-linear function
    Eigen::VectorXd lambda_;                      //Lagrange multipliers
    mutable Eigen::VectorXd Cx_;                  //value of C*x
    mutable Eigen::VectorXd Cl_;                  //value of C^T*lambda
    mutable Eigen::VectorXd Jx_;                  //value of J*x
    mutable Eigen::VectorXd JtJx_;                //value of J^T * J *x
    mutable Eigen::VectorXd Jp_;                  //value of J*p
    mutable Eigen::VectorXd gradL_;               //gradient of the Lagrangian
    std::vector<Activation> currentActiveSet_;    //current active set
    std::vector<Activation> previousActiveSet_;   //active set backup
    int k_;                                       //iteration number
  };
}
