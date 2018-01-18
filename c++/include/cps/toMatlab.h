#pragma once
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

#include <Eigen/Dense>
#include <iosfwd>


/** A small utility class to write Eigen matrices in a stream with a matlab-readable format.
  *
  * Example of use:
  * Eigen::MatrixXd M = Eigen::MatrixXd::Random(5,6);
  * std::cout << (toMatlab)M << std::endl;
  *
  * Inspired from a code from Nicolas Mansard.
  */
class toMatlab
{
public:
  template<typename Derived>
  toMatlab(const Eigen::DenseBase<Derived>& M)
    : mat(M)
  {}

private:
  Eigen::MatrixXd mat;

  friend std::ostream& operator<< (std::ostream&, const toMatlab&);
};

inline std::ostream& operator<< (std::ostream& o, const toMatlab& tom)
{
  if (tom.mat.cols() == 1)
  {
    Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";", "", "", "[", "]");
    o << tom.mat.format(fmt);
  }
  else
  {
    Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    o << tom.mat.format(fmt);
  }
  return o;
}