import dolfin.cpp as cpp
from dolfin import DOLFIN_EPS
from dolfin import*

BOUNDARY_MARGIN = 1.e-3 + DOLFIN_EPS


class Boundary(cpp.mesh.SubDomain):
    def inside(self, x, on_boundary):
        raise NotImplemented

class BoundaryConditions:
    def __init__(self, mesh: Mesh, velocity_space, pressure_space):
        self.mesh = mesh
        self.velocity_space = velocity_space
        self.pressure_space = pressure_space

    def velocity_boundary_conditions(self):
        raise NotImplementedError

    def pressure_boundary_conditions(self):
        raise NotImplementedError

    @property
    def get_boundary_conditions(self):
        return self.velocity_boundary_conditions(), self.pressure_boundary_conditions()


class BackwardFacingStepBoundaryConditions(BoundaryConditions)

    def create_inflow(self):
        class InFlowBC(Boundary):
            def inside(self, x, on_boundary):
                return on_boundary and x[0] < 0 + bmarg
        # Create inflow boundary condition for velocity with parabolic inflow
        g0 = Expression(("-amp*(1.0 - x[1])*(2.0 - x[1])", "0.0"), amp=0.25)
        inflowboundary = OutflowBoundary()
        inflow_u = DirichletBC(V, g0, inflowboundary)

    def velocity_boundary_conditions(self):



def backward_facing_step_boundary_conditions(mesh: Mesh,
                                             velocity_space: VectorFunctionSpace,
                                             pressure_space: FunctionSpace):


    class InflowBoundary(Boundary):
            def inside(self, x, on_boundary):
                return on_boundary and x[0] > 8 - BOUNDARY_MARGIN

    # No-slip boundary
    class NoslipBoundary(Boundary):
            def inside(self, x, on_boundary):
                return on_boundary \
                and (x[1] < 0 + BOUNDARY_MARGIN \
                or x[1] > 2 - BOUNDARY_MARGIN \
                or (x[0] < 1. + BOUNDARY_MARGIN and x[1] < 1. + BOUNDARY_MARGIN) \
                or (x[1] < 1. + BOUNDARY_MARGIN and x[0] < 1. + BOUNDARY_MARGIN))

    # Outflow boundary
    class OutflowBoundary(Boundary):
            def inside(self, x, on_boundary):
                return on_boundary and x[0] < 0 + bmarg


def boundary_condition_wrapper(mesh):