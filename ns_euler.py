from dolfin import *
from block import *
from block.iterative import *
from block.algebraic.petsc import *
from block.dolfin_util import *
from block.iterative import *
import numpy as np
import subprocess
import os
from instant import*
from flufl import*
import numpy as np
import scipy.io as sio
import scipy
import matplotlib.pyplot as p
import sys
import time
import ns_make as m

def main(refinement_level,nu,prec_choice,bplot,T):


	# change instant cache dir - dolfin matrices will be saved here
	os.environ["INSTANT_CACHE_DIR"] = "/scratch/instant1/.instant"
	os.environ["INSTANT_ERROR_DIR"] = "/scratch/instant1/.instant"
	# Print log messages only from the root process in parallel
	parameters["std_out_all_processes"] = False;



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

	##########################################################
	method = 'Euler'
	# set navier stokes parameters
	f = Constant((0, 0))
	#T  = 5.    	# total simulation time
	dt = 0.0125   		# time step
	#nu = 0.001  		# viscosity 
		
		# steady state with nu = 1/150 
		# nu = 1/300 is below a critival values -> not stable

	step = 0
	t = dt	
	ttictoc = []
	
	# define trial and test functions
	u = TrialFunction(V)
	p = TrialFunction(Q)
	v = TestFunction(V)
	q = TestFunction(Q)


	a12 = dt*div(v)*p*dx
	a21 = dt*div(u)*q*dx


	path = m.saveall(method,prec_choice,refinement_level,nu,T,dt)
	res=[]
	if not os.path.exists(path+'/residuals'):
        	os.makedirs(path+'/residuals')


	###############################################
	# velocity mass matrix
	Mv = assemble(inner(u,v)*dx)		
	MvInvDiag = InvDiag(Mv)
	###############################################	
	ufile=File(path + "/velocity.pvd")
	vis_f = File(path + "/solution_navier.pvd")
	uold = Function(V)


	# result
	ueuler = Function(V)
	peuler = Function(Q)


	# only one image
	uplot = Function(V)
	pplot = Function(Q)



	while t <= T + DOLFIN_EPS:   
		step +=1
	
		at = inner(u,v)*dx + dt*nu*inner( grad(u), grad(v) )*dx + dt*inner( grad(u)*Function(V,uold), v )*dx 
		L2 = inner(Function(V,uold),v)*dx + dt*inner(f,v)*dx
	
		
		# build block matrices	
		AA = block_assemble([[at, a12],
				[a21,  0 ]], bcs=bcs)
		bb2  = block_assemble([L2, 0], bcs=bcs)
		[[At, B],
		[C, _]] = AA
		

	
		####################################	
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


		
		#################################

		
		# Create the block inverse, using the BiCGStab method 
		AAinv = BiCGStab(AA, precond=prec, tolerance=1e-8, maxiter=100, show=2, callback = True)
		
		tic = time.clock()
		ueuler,peuler = AAinv * bb2
	
		toc = time.clock()
		ttictoc.append(toc - tic)
		tt = np.array(ttictoc)
		sio.savemat(path + '/tictoc.mat', {'tt':tt})		

		# save residuals
		res = AAinv.residuals	
		sio.savemat(path + '/residuals/residuals_%d.mat'%step, {'res':res})

	    
		uold = ueuler	    
		print 'time =', t
		t += dt	   
		
		##############
		# naming problem for visit
		ueuler.rename("u","u")
    		Ftmp = Function(V,ueuler)
    		Ftmp.rename("Visit Velocity", "Visit Velocity")
    		vis_f << Ftmp

		ela_time = t/T * 100
		print "\n \nelapsed time: %.2f in percent with %i steps\n \n" % (ela_time, step)
		
		if bplot == 1:
			# plot velocity
			uplot.assign(Function(V,uold))
			plot(uplot,title = "Velocity", interactive = False)
	
			# plot pressure
			pplot.assign(Function(Q,peuler))
			plot(pplot, title = "Pressure")
	



