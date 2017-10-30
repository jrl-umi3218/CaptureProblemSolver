#include <boost/python.hpp>
// #include <boost/python/numpy.hpp>

#include "Problem.h"

namespace py = boost::python;


BOOST_PYTHON_MODULE(PyBalanceMPCSolver)
{
  // Py_Initialize();
  using namespace bms;

  py::class_<RawProblem>("RawProblem")
    .def("read", &RawProblem::read);
}
