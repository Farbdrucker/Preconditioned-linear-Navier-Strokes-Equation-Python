from block import *

from block.algebraic.petsc import *
from block.dolfin_util import *

from instant import*
from flufl import*

def least_square_commutator(lhs, velocity_trial_function, velocity_test_function):
    [[At, B],
     [C, _]] = lhs
    # LSC - Least Square Commutator
    # build Mf
    Mf = ILU(At)

    # build Ms, schur complement
    Mv = assemble(inner(velocity_trial_function, velocity_test_function) * dx)
    MvInvDiag = InvDiag(Mv)

    BTB = ML(collapse(C * MvInvDiag * B))
    BFB = collapse(C * MvInvDiag * At * MvInvDiag * B)
    Ms = BTB * BFB * BTB

    preconditioner = block_mat([[Mf, 0],
                      [C, -Ms]]).scheme('gs')

    return preconditioner


def simple_preconditioner(lhs):
    [[At, B],
     [C, _]] = lhs
    Mf = ILU(At)
    Ms = ML(collapse(C * InvDiag(At) * B))

    preconditioner = block_mat([[Mf, 0],
                      [C, -Ms]]).scheme('sgs')
    return preconditioner


def preconditioner(scheme: str):
    # Least Square Commutator
    if scheme == 'lsc':
        return least_square_commutator
    if scheme == 'simple'
        return simple_preconditioner