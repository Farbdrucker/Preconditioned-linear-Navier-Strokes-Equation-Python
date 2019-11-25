from util import mesh_loader
from boundaries import
class Simulation:
    def __init__(self):
        self.mesh = None
        self.boundary_conditions = None
        self.space = None


class SimulationBuilder:
    def __init__(self):
        self.simulation = Simulation()

    def load_mesh(self, path: str):
        self.simulation.mesh = mesh_loader(path)

    def load_spaces(self):
        self.simulation.space = function_space_wrapper(self.simulation.mesh)

    def initiate_boundary_conditions(self):
        self.simulation.boundary_conditions = boundary_conditions

    @property
    def get_trial_functions(self):
        velocity = TrialFunction(self.simulation.space.vector)
        pressure = TrialFunction(self.simulation.space.function)
        return velocity, pressure

    @property
    def get_test_functiosn(self):
        velocity = TestFunction(self.simulation.space.vector)
        pressure = TestFunction(self.simulation.space.function)
        return velocity, pressure

    @property
    def get_velocity_function(self):
        return Function(self.simulation.space.vector)

    @property
    def get_pressure_function(self):
        return Function(self.simulation.space.function)

    def create_environment(self, mesh_dir):
        self.load_mesh(mesh_dir)
        self.load_spaces()
        self.initiate_boundary_conditions()
        return self.simulation
