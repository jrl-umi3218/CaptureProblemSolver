#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "converters.h"

#include "Problem.h"

namespace py = boost::python;
namespace np = boost::python::numpy;


BOOST_PYTHON_MODULE(PyBalanceMPCSolver)
{
  Py_Initialize();
  np::initialize();
  pygen::convert<double>(pygen::Converters::Vector);

  using namespace bms;

  py::class_<RawProblem>("RawProblem")
    .def("read", &RawProblem::read)
    .add_property("delta", 
        py::make_getter(&RawProblem::delta, 
          py::return_value_policy<py::copy_non_const_reference>()),
        py::make_setter(&RawProblem::delta))
    .def_readwrite("g", &RawProblem::g)
    .def_readwrite("g", &RawProblem::g)
    .def_readwrite("lmin", &RawProblem::lmin)
    .def_readwrite("lmax", &RawProblem::lmax)
    .def_readwrite("wi_min", &RawProblem::wi_min)
    .def_readwrite("wi_max", &RawProblem::wi_max)
    .def_readwrite("zi", &RawProblem::zi)
    .def_readwrite("dzi", &RawProblem::dzi)
    .def_readwrite("zf", &RawProblem::zf);
}
