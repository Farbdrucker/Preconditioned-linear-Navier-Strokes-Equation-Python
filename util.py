from  dolfin import VectorFunctionSpace, FunctionSpace
from dolfin.cpp.mesh import Mesh


def mesh_loader(path: str)->Mesh:
    if path.endswith('.xml'):
        return Mesh(path)
    elif path.endswith('.xml.gz'):
        return Mesh(path)
    else:
        raise NotImplementedError


def function_space_wrapper(mesh: Mesh):
    """

    Args:
        mesh:

    Returns:

    """
    class Spaces:
        def __init__(self,
                     space_mesh,
                     space_optimizer,
                     vector_space_degree,
                     function_space_degree):

            self.vector = VectorFunctionSpace(space_mesh, space_optimizer, vector_space_degree)
            self.function = FunctionSpace(space_mesh, space_optimizer, function_space_degree)

    return Spaces(mesh, 'GC', 2, 1)




