from block.iterative import LGMRES


class Solver:
    def __init__(self, tolerance: float, maxiter: int, debug: int):
        self.tol = tolerance
        self.maxiter = maxiter
        self.debug = debug

    def solve(self, lhs, rhs, prec):
        raise NotImplementedError


class LGMRESSolver(Solver):
    def solve(self, lhs, rhs, prec):
        lhs_inv = LGMRES(lhs, precond=prec,
                         tolerance=self.tol,
                         maxiter=self.maxiter,
                         show=self.debug)

        velocity, pressure = lhs_inv * rhs

        return velocity, pressure

