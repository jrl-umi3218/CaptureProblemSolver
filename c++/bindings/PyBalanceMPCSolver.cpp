#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "converters.h"

#include "LeastSquare.h"
#include "Problem.h"
#include "QuadraticObjective.h"
#include "SQP.h"

namespace py = boost::python;
namespace np = boost::python::numpy;


BOOST_PYTHON_MODULE(PyBalanceMPCSolver)
{
  Py_Initialize();
  np::initialize();
  pygen::convert<double>(pygen::Converters::Vector);

  using namespace bms;

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
    .def("set_dzi", &Problem::set_dzi)
    .def("set_lambda_max", &Problem::set_lambda_max)
    .def("set_lambda_min", &Problem::set_lambda_min)
    .def("set_lambdas", &Problem::set_lambdas)
    .def("set_wi", &Problem::set_wi)
    .def("set_wi_max", &Problem::set_wi_max)
    .def("set_wi_min", &Problem::set_wi_min)
    .def("set_zf", &Problem::set_zf)
    .def("set_zi", &Problem::set_zi);

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
    .def_readwrite("lambda_min", &RawProblem::lmin)
    .def_readwrite("lambda_max", &RawProblem::lmax)
    .def_readwrite("omega_i_min", &RawProblem::wi_min)
    .def_readwrite("omega_i_max", &RawProblem::wi_max)
    .def_readwrite("z_bar", &RawProblem::zi)
    .def_readwrite("zd_bar", &RawProblem::dzi)
    .def_readwrite("z_f", &RawProblem::zf);

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
    .def("solve", &SQP::solve)
    .def("x", &SQP::x, py::return_value_policy<py::copy_const_reference>());
}
