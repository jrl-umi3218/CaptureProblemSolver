from BalanceMPCSolver import Problem, RawProblem, SQP


def load_problem(number):
    dirname = "/home/tastalian/Projects/BalanceMPCSolver/c++/tests/data/"
    basename = "Problem%02d.txt" % number
    raw_pb = RawProblem()
    raw_pb.read(dirname + basename)
    return Problem(raw_pb)


if __name__ == "__main__":
    pb = load_problem(1)
    sqp = SQP(pb.size)
    sqp.solve(pb)
