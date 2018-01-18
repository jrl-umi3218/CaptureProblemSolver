/* Copyright 2018 CNRS-AIST JRL, CNRS-UM LIRMM
 *
 * This file is part of CPS.
 *
 * CPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with CPS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cps/ProblemMatrices.h>

using namespace Eigen;

namespace cps
{
  Eigen::MatrixXd buildJ(const Eigen::VectorXd& d)
  {
    auto n = d.size();
    MatrixXd J(n - 1, n);
    J.setZero();
    J.diagonal() = -d.head(n - 1) - d.tail(n - 1);
    J.diagonal<-1>() = d.segment(1, n - 2);
    J.diagonal<+1>() = d.segment(1, n - 1);

    return J;
  }

  Eigen::MatrixXd buildJ0(const Eigen::VectorXd& e)
  {
    auto n = e.size();
    MatrixXd J0(n+1,n);
    J0.setZero();
    J0.diagonal() = e;
    J0.diagonal<-1>().head(n - 1) = -e.head(n - 1) - e.tail(n - 1); 
    J0(n, n - 1) = -e(n - 1);
    J0.diagonal<-2>() = e.tail(n - 1);

    return J0;
  }

  Eigen::MatrixXd buildJj(const Eigen::VectorXd& e)
  {
    auto n = e.size();
    MatrixXd Jj(n + 1, n + 1);
    Jj.setZero();
    Jj.diagonal().head(n) = -e;
    Jj.diagonal().tail(n) -= e;
    Jj.diagonal<-1>() = e;
    Jj.diagonal<+1>() = e;

    return Jj;
  }

  Eigen::MatrixXd buildJpm1(const Eigen::VectorXd& e)
  {
    auto n = e.size();
    MatrixXd Jp(n, n);
    Jp.setZero();
    Jp.diagonal() = -e;
    Jp.diagonal().tail(n - 1) -= e.head(n - 1);
    Jp.diagonal<-1>() = e.head(n - 1);
    Jp.diagonal<+1>() = e.head(n - 1);
    
    return Jp;
  }

  Eigen::MatrixXd buildCZ(int n)
  {
    MatrixXd Cz(n, n);
    Cz.setZero();
    Cz.diagonal().setOnes();
    Cz.diagonal<-1>().setConstant(-1);

    return Cz;
  }

  Eigen::MatrixXd buildC(int n)
  {
    MatrixXd C(n+1, n);
    C.setZero();
    C.diagonal().setOnes();
    C.diagonal<-1>().setConstant(-1);
    C(n, n - 1) = 1;

    return C;
  }
}
