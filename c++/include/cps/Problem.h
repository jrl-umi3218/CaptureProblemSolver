#pragma once

#include <string>

#include <Eigen/Core>

#include <cps/cps_api.h>
#include <cps/BoundenessConstraint.h>
#include <cps/LinearConstraints.h>
#include <cps/QuadraticObjective.h>


namespace cps
{
  /** Description of a capture problem parameters.
    * Fields name correspond to the notations in the paper.
    */
  struct CPS_DLLAPI RawProblem
  {
    /** Read from file (see examples in test/data/.*/
    void read(const std::string& filepath);

    double g = 9.80665; // ISO 80000-3
    double lambda_min;
    double lambda_max;
    Eigen::VectorXd delta;
    double init_omega_min;
    double init_omega_max;
    double init_zbar;         //z_i
    double init_zbar_deriv;   //dz_i/dt
    double target_height;     //z_f

    //solution (optional)
    Eigen::VectorXd Phi_;
  };


  /** A class representing a capture problem to be passed to the SQP solver.
    */
  class CPS_DLLAPI Problem
  {
  public:
    /** Build from a RawProblem.*/
    Problem(const RawProblem& pb);

    /** Get the underlying LeastSquareObjective instance.*/
    const LeastSquareObjective& objective() const;
    /** Get the underlying LeastSquareObjective instance.*/
    LeastSquareObjective& objective();
    /** Get the underlying BoundenessConstraint instance.*/
    const BoundenessConstraint& nonLinearConstraint() const;
    /** Get the underlying BoundenessConstraint instance.*/
    BoundenessConstraint& nonLinearConstraint();
    /** Get the underlying LinearConstraints instance.*/
    const LinearConstraints& linearConstraints() const;
    /** Get the underlying LinearConstraints instance.*/
    LinearConstraints& linearConstraints();
    /** Get the size n of the problem.*/
    Eigen::VectorXd::Index size() const;

    /** Precompute QR decompositions of the projected-nullspace jacobians of
      * the least-square objective for all possible active set value.
      *
      * Note that these pre-computations depend only on the `delta` vector of
      * the internal RawProblem. All other parameters (`init_zbar`,
      * `target_height`, etc.) can be changed with no need to call this function
      * again.
      */
    void precompute();

    /** Change omega_i_min and omega_i_max.*/
    void set_init_omega(double init_omega_min, double init_omega_max);
    /** Change omega_i_max.*/
    void set_init_omega_max(double init_omega_max);
    /** Change omega_i_min.*/
    void set_init_omega_min(double init_omega_min);
    /** Change z_i.*/
    void set_init_zbar(double init_zbar);
    /** Change dz_i/dt.*/
    void set_init_zbar_deriv(double init_zbar_deriv);
    /** Change lambda_max.*/
    void set_lambda_max(double lambda_max);
    /** Change lambda_min.*/
    void set_lambda_min(double lambda_min);
    /** Change lambda_min and lambda_max.*/
    void set_lambdas(double lambda_min, double lambda_max);
    /** Change z_f*/
    void set_target_height(double target_height);

    /** Get delta.*/
    const Eigen::VectorXd& delta() const { return raw_.delta; }
    /** Get omega_i_max.*/
    double init_omega_max() const { return raw_.init_omega_max; }
    /** Get omega_i_min.*/
    double init_omega_min() const { return raw_.init_omega_min; }
    /** Get z_i.*/
    double init_zbar() const { return raw_.init_zbar; }
    /** Get dz_i/dt.*/
    double init_zbar_deriv() const { return raw_.init_zbar_deriv; }
    /** Get lambda_max.*/
    double lambda_max() const { return raw_.lambda_max; }
    /** Get lambda_min.*/
    double lambda_min() const { return raw_.lambda_min; }
    /** Get z_f.*/
    double target_height() const { return raw_.target_height; }

  private:
    void computeAndSetBounds0();
    void computeAndSetZonotopeBounds();
    void computeAndSetBoundsN();
    void computeAndSetAlpha();
    void computeAndSetb();

    LeastSquareObjective lso_;
    LinearConstraints lc_;
    BoundenessConstraint bc_;

    RawProblem raw_;
  };
}
