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
