#pragma once

#include <Eigen/Core>

#include <bms_api.h>
#include <LeastSquare.h>
#include <Problem.h>
#include <QuadraticObjective.h>
#include <Statistics.h>

namespace bms
{
  class BMS_DLLAPI SQP
  {
  public:
    class Parameters
    {
    public:
      int maxIter()           const { return maxIter_; }
      double mu()             const { return mu_; }
      double beta()           const { return beta_; }
      double c1()             const { return c1_; }
      double smallestLSStep() const { return smallestLSStep_; }
      double tau_p()          const { return tau_p_; }
      double tau_d()          const { return tau_d_; }
      double feasibilityEps() const { return feasibilityEps_; }

      Parameters& maxIter       (int i)    { maxIter_ = i;        return *this; }
      Parameters& mu            (double d) { mu_ = d;             return *this; }
      Parameters& beta          (double d) { beta_ = d;           return *this; }
      Parameters& c1            (double d) { c1_ = d;             return *this; }
      Parameters& smallestLSStep(double d) { smallestLSStep_ = d; return *this; }
      Parameters& tau_p         (double d) { tau_p_ = d;          return *this; }
      Parameters& tau_d         (double d) { tau_d_ = d;          return *this; }
      Parameters& feasibilityEps(double d) { feasibilityEps_ = d; return *this; }

    private:
      int maxIter_ = 100;
      double mu_ = 100000;           //penalty parameter
      double beta_ = 0.9;            //backtracking multiplier in line search
      double c1_ = 0.01;             //gradient coefficient in line search
      double smallestLSStep_ = 1e-8;
      double tau_p_ = 1e-6;          //precision parameter on primal condition
      double tau_d_ = 1e-6;          //precision parameter on dual condition
      double feasibilityEps_ = 1e-8; //margin of feasibility for which a potential numerical issue is detected
    };

    SQP(int n);
    void SQPParameters(const Parameters& p);
    void LSParameters(const LeastSquare::Parameters& p);
    const Parameters& SQPParameters() const;
    const LeastSquare::Parameters& LSParameters() const;

    SolverStatus solve(const Problem& pb);
    SolverStatus solveFeasibility(const Problem& pb);

    const Eigen::VectorXd& x() const;
    const Eigen::VectorXd& lambda() const;
    const std::vector<Activation>& activeSet() const;
    int numberOfIterations() const;

    const stats::SQPStats& statistics() const;

  private:
    bool checkKKT(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, double f, const Eigen::VectorXd& g, 
                  const LinearConstraints& lc, const LeastSquareObjective* const obj = nullptr) const;

    Eigen::DenseIndex n_;
    LeastSquare ls_;

    //optimization parameters
    Parameters params_;

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

    //stats
    stats::SQPStats stats_;
  };
}
