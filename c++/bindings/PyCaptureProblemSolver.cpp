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

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include <cps/LeastSquare.h>
#include <cps/Problem.h>
#include <cps/QuadraticObjective.h>
#include <cps/SQP.h>

#include "converters.h"

namespace py = boost::python;
namespace np = boost::python::numpy;

BOOST_PYTHON_MODULE(PyCaptureProblemSolver)
{
  Py_Initialize();
  np::initialize();
  pygen::convert<double>(pygen::Converters::Vector);

  using namespace cps;

  py::class_<LeastSquare>("LeastSquare", py::init<int>());

  py::class_<LeastSquareObjective>("LeastSquareObjective", py::init<const Eigen::VectorXd&>());

  LeastSquareObjective& (Problem::*pyObjective)() = &Problem::objective;
  BoundenessConstraint& (Problem::*pyNonLinearConstraint)() = &Problem::nonLinearConstraint;
  LinearConstraints& (Problem::*pyLinearConstraints)() = &Problem::linearConstraints;
  py::class_<Problem>("Problem", py::init<RawProblem>())
    .add_property("size", &Problem::size)
    .def("linear_constraints", pyLinearConstraints, py::return_internal_reference<>())
    .def("nonlinear_constraint", pyNonLinearConstraint, py::return_internal_reference<>())
    .def("objective", pyObjective, py::return_internal_reference<>())
    .def("precompute", &Problem::precompute)
    .def("set_init_omega", &Problem::set_init_omega)
    .def("set_init_omega_max", &Problem::set_init_omega_max)
    .def("set_init_omega_min", &Problem::set_init_omega_min)
    .def("set_init_zbar", &Problem::set_init_zbar)
    .def("set_init_zbar_deriv", &Problem::set_init_zbar_deriv)
    .def("set_lambda_max", &Problem::set_lambda_max)
    .def("set_lambda_min", &Problem::set_lambda_min)
    .def("set_lambdas", &Problem::set_lambdas)
    .def("set_target_height", &Problem::set_target_height);

  py::class_<RawProblem>("RawProblem")
    .def("read", &RawProblem::read)
    .add_property("delta", 
        py::make_getter(&RawProblem::delta, 
          py::return_value_policy<py::copy_non_const_reference>()),
        py::make_setter(&RawProblem::delta))
    .add_property("Phi_", 
        py::make_getter(&RawProblem::Phi_, 
          py::return_value_policy<py::copy_non_const_reference>()),
        py::make_setter(&RawProblem::Phi_))
    .def_readwrite("g", &RawProblem::g)
    .def_readwrite("lambda_min", &RawProblem::lambda_min)
    .def_readwrite("lambda_max", &RawProblem::lambda_max)
    .def_readwrite("init_omega_min", &RawProblem::init_omega_min)
    .def_readwrite("init_omega_max", &RawProblem::init_omega_max)
    .def_readwrite("init_zbar", &RawProblem::init_zbar)
    .def_readwrite("init_zbar_deriv", &RawProblem::init_zbar_deriv)
    .def_readwrite("target_height", &RawProblem::target_height);

  py::enum_<SolverStatus>("SolverStatus")
    .value("Converge", SolverStatus::Converge)
    .value("MaxIteration", SolverStatus::MaxIteration)
    .value("LineSearchFailed", SolverStatus::LineSearchFailed)
    .value("NoLinearlyFeasiblePoint", SolverStatus::NoLinearlyFeasiblePoint)
    .value("NumericallyEquivalentIterates", SolverStatus::NumericallyEquivalentIterates)
    .value("Fail", SolverStatus::Fail);

  py::class_<SQP>("SQP", py::init<int>())
    .add_property("nb_iter", &SQP::numberOfIterations)
    .def("lambda_", &SQP::lambda, py::return_value_policy<py::copy_const_reference>())
    .def("numberOfIterations", &SQP::numberOfIterations)
    .def("solve", &SQP::solve)
    .def("x", &SQP::x, py::return_value_policy<py::copy_const_reference>());
}
