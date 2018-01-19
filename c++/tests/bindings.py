# Copyright 2018 CNRS-AIST JRL, CNRS-UM LIRMM
#
# This file is part of CPS.
#
# CPS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CPS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with CPS.  If not, see <http://www.gnu.org/licenses/>.
#

import os

from CaptureProblemSolver import Problem, RawProblem, SQP


def load_problem(number):
    dirname = os.path.dirname(os.path.realpath(__file__)) + "/data/"
    basename = "Problem%02d.txt" % number
    raw_pb = RawProblem()
    raw_pb.read(dirname + basename)
    return Problem(raw_pb)


if __name__ == "__main__":
    pb = load_problem(1)
    sqp = SQP(pb.size)
    sqp.solve(pb)
    print "SQP solution:", sqp.x()
