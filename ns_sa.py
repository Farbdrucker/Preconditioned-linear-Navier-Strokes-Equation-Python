from typing import List
from dolfin import *
from block import *
from block.iterative import *
from block.algebraic.petsc import *
from block.dolfin_util import *
from block.iterative import *
import numpy as np

import os
from instant import*
from flufl import*
import scipy.io as sio

import time
import ns_make as m











class FlowSimulation:
    def __init__(self):
        builder = SimulationBuilder()
        self.environment = builder.create_environment()

        self.velocity

    def step_scheme(self):
        """
        implement discretization scheme to solve linearized Navier-Stokes equation
        return corresponding blocks A of mass matrix
            [[A B],
            [B.T 0]]
        Returns:

        """
        raise NotImplementedError


    def build_block(self):
        a11 = self.step_scheme()
        lhs = block_assemble([[at, self.a12],
                                [self.a21, 0]],
                                  bcs=self.environment.boundary_conditions)

        bb2_e = block_assemble([L2, 0],
                               bcs=self.environment.boundary_conditions)



    def solve(self):
        raise NotImplementedError





# specify simulation time and step size

# Inital Step

# Build Block
    # Simo Armero Scheme

    # Backward Euler Scheme

    # Time Stepping Scheme

# Build Preconditioner
    # Simple
    # LSC

# Solve


def main(refinement_level,nu,prec_choice,bplot,T):
	# Print log messages only from the root process in parallel
	parameters["std_out_all_processes"] = True; 	# False


	# load refined mesh
	if refinement_level > 8:
		raise RuntimeError, "No mesh available for refinement level %d" % refinement_level

	mesh = Mesh("data/cylinder_%d.xml.gz" % refinement_level)

	if bplot == 1:
		plot(mesh, title="Finite element mesh")


	# Define function spaces (P2-P1)
	V = VectorFunctionSpace(mesh, "CG", 2)
	Q = FunctionSpace(mesh, "CG", 1)

	# create boundary conditions
	bcs = m.cy_bc(mesh,V,Q)

	###############################
	# discretization in time scheme
	method = 'SA'	# Simo Armero

	# choice of preconditioning, 
	# choose lsc for constant good prec over different nu and rfnmlvl
	#prec_choice = 'lsc'

	# set navier stokes parameters
	f = Constant((0, 0))
	#T  = 5.    		# total simulation time
	dt = 0.0125 		# time stepsize - should be smaller then a cell diameter

	#nu = 0.001  		# viscosity 
				# viscosity  nu = 0.001 (p.429 Fenics boook) or Turek (1996).

	ttictoc = []
	# create a path - /scratch
	path = m.saveall(method,prec_choice,refinement_level,nu,T,dt)

	# initial array for residuals
	res=[]

	# add dir for res
	if not os.path.exists(path+'/residuals'):
        	os.makedirs(path+'/residuals')


	print "\n\nTimestep dt = %e\n" %dt
	print "Viscosity nu = %e\n\n" %nu


	# count steps
	step = 0
	t = dt	


	########################################
	# Define trial and test functions
	u = TrialFunction(V)
	p = TrialFunction(Q)
	v = TestFunction(V)
	q = TestFunction(Q)
	
	a12 = 0.5*dt*inner(div(v),p)*dx
	a21 = 0.5*dt*inner(div(u),q)*dx

	###################################
	# mass matrices for preconditioning

	Mv = assemble(inner(u,v)*dx)		
	MvInvDiag = InvDiag(Mv)
	####################################	
	ufile=File(path + "/velocity.pvd")

	u_1 = Function(V)
	u_2 = Function(V)

	p1 = Function(Q)

	# only one image
	uplot = Function(V)
	pplot = Function(Q)

	# Solution function
	uu = Function(V)
	pp = Function(Q)

	# save velocity vectors
	vis_f = File(path + "/solution_navier.pvd")

	##############################################
	# NOT SELF STARTING --> EULER STEP (u_1 = 0 .> Stokes problem)
	dt_e = 0.1*dt

	at = inner(u,v)*dx + dt_e*nu*inner( grad(u), grad(v) )*dx + dt_e*inner( grad(u)*Function(V,u_1), v )*dx 
	L2 = inner(Function(V,u_1),v)*dx + dt_e*inner(f,v)*dx

	AA_e = block_assemble([[at, a12],
			[a21,  0 ]], bcs=bcs)

	bb2_e  = block_assemble([L2, 0], bcs=bcs)

	[[At, B],
	[C, _]] = AA_e

	#################	
	# build precond
	# build Mf
	Mf  = ML(At)

	# build Ms
	MvInvDiag = InvDiag(Mv)
	BTB = ML(collapse(C * MvInvDiag * B))
	BFB = C * collapse(MvInvDiag * At * MvInvDiag) * B
	Ms = (BTB*BFB*BTB)


	prec = block_mat([[Mf, 0],
                 [C, -Ms]]).scheme('gs') 
		
	################

	AAinv_e = LGMRES(AA_e,precond=prec, tolerance=1e-8, maxiter=50, show=2)

	uu,pp = AAinv_e * bb2_e

	p1.assign(Function(Q,pp))
	u_1.assign(Function(V,u_2))
	u_2.assign(Function(V,uu))
	u_2.rename("u","u")

	##################################################


	while t <= T + DOLFIN_EPS:   

		step +=1

		u_st  = 1.5*u_2 - 0.5*u_1
		u_bar = 0.5*(u + u_2)

		# lhs / rhs splitting method
		F = inner(u-u_2,v)*dx + dt*nu*inner(grad(u_bar),grad(v))*dx + dt*inner(grad(u_bar)*u_st,v)*dx - 0.5*dt*inner(div(v),p1)*dx
		at = lhs(F); L2 = rhs(F)


		# build block matrices	
		AA = block_assemble([[at, a12],
			     	     [a21,  0 ]], bcs=bcs)

		bb2  = block_assemble([L2, 0], bcs=bcs)
	
		# assign blocks
		[[At, B],
	 	[C, _]] = AA
	

	
		##################################	
		# build precond
		if prec_choice == 'LSC':
			# LSC - Least Square Commutator
			# build Mf
			Mf  = ILU(At)

		# build Ms, schur complement
	
			BTB = ML(collapse(C * MvInvDiag * B))
			BFB =  collapse(C *MvInvDiag * At * MvInvDiag* B) 
			Ms = BTB*BFB*BTB


			prec = block_mat([[Mf, 0],
                  		[C, -Ms]]).scheme('gs')

		if prec_choice == 'SIMPLE':
			# SIMPLE
			Mf = ILU(At)
			Ms = ML(collapse(C*InvDiag(At)*B))

		
			prec = block_mat([[Mf, 0],
                  		[C, -Ms]]).scheme('sgs') 



		############################


		# Create the block inverse
		Ainv = BiCGStab(AA, precond=prec, tolerance=1e-8, maxiter=100, show=2)
				#show = 3 -> shows the residual in iteration-> very slow
		
		tic = time.clock()

		uu,pp = Ainv * bb2


		toc = time.clock()
		ttictoc.append(toc - tic)
		tt = np.array(ttictoc)
		sio.savemat(path + '/tictoc.mat', {'tt':tt})
	
		# save matrices - only one time - just to take a look at sparse
		# -> matlab mmread.m
		if step == 3:

			mfile = File (path + "/matrices/At.mtx")
			mfile << At
			mfile = File(path + "/matrices/B.mtx")
			mfile << B
			mfile = File (path + "/matrices/C.mtx")
			mfile << C
		

		# save residuals in an array and save it as .mat for matlab
		res = Ainv.residuals
		sio.savemat(path + '/residuals/residuals_%d.mat'%step, {'res':res})

		# assign velocity
		u_1.assign(Function(V,u_2))   
		u_2.assign(Function(V,uu))	
		p1.assign(Function(Q,pp))
    
		print 'time =', t
		t += dt	   


		##############
		# naming problem for visit
		u_2.rename("u","u")
    		Ftmp = Function(V,u_2)
    		Ftmp.rename("Visit Velocity", "Visit Velocity")
    		vis_f << Ftmp


		ela_time = t/T * 100
		print "\n \nelapsed time: %.2f in percent with %i steps\n \n" % (ela_time, step)
	
		if bplot== 1:
			# plot velocity
			uplot.assign(Function(V,uu))
			plot(uplot,title = "Velocity (progress %.2f %%)" %ela_time)
	
			# plot pressure
			pplot.assign(Function(Q,p1))
			plot(pplot, title = "Pressure")

	



