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
import matplotlib.pyplot as p
import sys
import time
import ns_make as m
# compare [KGGS10] and ifiss3.4 timestepping stabtrNS.m


def main(refinement_level,nu,prec_choice,bplot,tol,T):
	
	# Print log messages only from the root process in parallel
	parameters["std_out_all_processes"] = False

	# deeper recursive function possibility
	sys.getrecursionlimit()
	sys.setrecursionlimit(1000000)
	
	print("Instant dir:", get_instant_dir())
	print("Default cache dir:", get_default_cache_dir())


	os.environ["INSTANT_CACHE_DIR"] = "/scratch/instant1" #/.instant
	os.environ["INSTANT_ERROR_DIR"] = "/scratch/instant1"	#/.instant
	#export INSTANT_CACHE_DIR=/scratch/instant1/.instant

	print("Instant dir:", get_instant_dir())
	print("Default cache dir:", get_default_cache_dir())


	flag = 0
	#load mesh
	#bf8 - backwardfacing step with length 8 for a natural outflow
	mesh = Mesh("data/bf8f%d.xml" %refinement_level)
	
	if bplot ==1:
		plot(mesh, title="Finite element mesh")

	# Define function spaces (P2-P1)
	V = VectorFunctionSpace(mesh, "CG", 2)
	Q = FunctionSpace(mesh, "CG", 1)


	# create boundary conditions
	bcs = m.bf_bc(mesh,V,Q)

	##########################################################
	# Set data for the problem, create a folder dict and save data in a file
	method = "TS bf"	# adaptive time stepping
	reason =''
	# set navier stokes parameters
	f = Constant((0, 0))	# inner force, set zero
	#T  = 150     	# total simulation time
	#nu = 1./50	# viscosity 
		
		
	# count steps
	step = 0

	#initial time stepping
	dt0 = 1.e-8		# SIAM ADAPTIVE TIMESTEPPING [KGGS10] dt0 = 1.e-8
	dt  = dt0
	t   = 0.

	# time tolerance
	#tol = 1.e-2		# SIAM ADAPTIVE TIMESTEPPING [KGGS10]tol = 1.e-4	
				# too small time steps for this setting-> tol = 1.e-2

	# periodic averaging
	nstar = 5000		# SIAM ADAPTIVE TIMESTEPPING [KGGS10] n = 10 periodic averaging

	# array for elapsed time in an BiCGStab
	ttictoc = []
	# array for time steps -> save for visualization
	time_steps = []
	failedstep = 0


	mode="standard"

	# create a path
	path = m.saveall(method,prec_choice,refinement_level,nu,T,dt0,tol)

	# initialize array for the residuals
	res=[]

	# a dir to save all res
	if not os.path.exists(path+'/residuals'):
        	os.makedirs(path+'/residuals')


##############################################

	# Define trial and test functions
	u = TrialFunction(V)
	p = TrialFunction(Q)
	v = TestFunction(V)
	q = TestFunction(Q)

	#for solve function
	d = TrialFunction(V)	#udot

###############################################
	# velocity mass matrices for preconditioning
	Mv = assemble(inner(u,v)*dx)		
	MvInvDiag = InvDiag(Mv)
###############################################	
	ufile = File(path + "/velocity.pvd")

	u1 = Function(V)
	u2 = Function(V)

# velocity estimation AB2 vs TR2
	e_h = Function(V)

	ud	= Function(V)
	udotb 	= Function(V)
	udiff	= Function(V)
	Mv_u	= Function(V)

	p1  	= Function(Q)


# only one image
	uplot 	= Function(V)
	pplot 	= Function(Q)

# Solution function
	dd 	= Function(V)
	uu 	= Function(V)
	pp 	= Function(Q)


	vis_f = File(path + "/solution_navier.pvd")

#########################################
# not self starting							
# solve (2.19) & (2.20)							
#	(2.9)  & (2.10)	from SIAM Adaptive time stepping [KGGS10]					
############################################

# NOT SELF STARTING --> EULER STEP (u_1 = 0 -> Stokes problem)
	dt_e = 0.1*dt

	at = inner(u,v)*dx + dt_e*nu*inner( grad(u), grad(v) )*dx + dt_e*inner( grad(u)*Function(V,u1), v )*dx 
	a12 = 0.5*dt*inner(div(v),p)*dx
	a21 = 0.5*dt*inner(div(u),q)*dx

	L2 = inner(Function(V,u1),v)*dx + dt_e*inner(f,v)*dx

	AA_e = block_assemble([[at, a12],
			[a21,  0 ]], bcs=bcs)

	bb2_e  = block_assemble([L2, 0], bcs=bcs)

	[[At, B],
	[C, _]] = AA_e

#################	
# build precond

# build Mf
	Mf  = ILU(At)

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


	u1.assign(Function(V,uu))


########################################
#(2.19).discrete (potential flow) problem
# solve
# (d,v) - (p0, div(v)) 	= -nu*(grad(u0),grad(v) - (u0*grad(u0), v)
#		(grad(d),q)  = 0.

	a11 = inner(d,v)*dx
	a12 = -div(v)*p*dx	 
	a21 = q*div(d)*dx

	F   = -nu*inner(grad(u1),grad(v))*dx - inner(grad(u1)*u1,v)*dx



##########
#block form
	AA1 = block_assemble([[a11, a12],
		      [a21,  0 ]], bcs = bcs)



	b1 = block_assemble([F, 0], bcs = bcs)


# Extract the individual submatrices
	[[A, B],
 	[C, _]] = AA1

# Create preconditioners: An ILU preconditioner for A, and an ML inverse of the
# Schur complement approximation for the (2,2) block.
	Ap = ILU(A)
	Dp = ML(collapse(C*InvDiag(A)*B))

	prec = block_mat([[Ap, B],
                  [C, -Dp]]).scheme('sgs')

	AA1inv = LGMRES(AA1,precond = prec, tolerance=1e-8, maxiter= 50, show=2)

	dd,pp = AA1inv*b1
	print "discrete potential flow problem solved..\n"
##############

	udotb.assign(Function(V,dd))

 

	ww = u1 + dt * udotb

#####################################
#(2.9)..								
# solve Oseen problem							
####################################
	a11 = 2.*inner(u,v)*dx + dt*nu*inner(grad(u),grad(v))*dx + dt*inner(grad(u)*ww,v)*dx 
	a12 = - dt*inner(p,div(v))*dx
	a21 = dt*inner(div(u),q)*dx


	F2 = 2.*inner(u1,v)*dx + dt*inner(udotb,v)*dx



################################
#block form
	AA2 = block_assemble([[a11 , a12],
		      [a21, 0]], bcs = bcs)

	b2 = block_assemble([F2, 0], bcs = bcs)

# Extract the individual submatrices
	[[A, B],
 	[C, _]] = AA2

# Create preconditioners: An ILU preconditioner for A, and an ML inverse of the
# Schur complement approximation for the (2,2) block.
	Ap = ILU(A)
	Dp = ML(collapse(C*InvDiag(A)*B))

	prec = block_mat([[Ap, B],
                  [C, -Dp]]).scheme('sgs')


	AA2inv = BiCGStab(AA2,precond = prec, tolerance=1e-8, maxiter=100, show=2)

	dd,pp = AA2inv*b2

	print "Oseen problem solved...\n"
##############################


	u2.assign(Function(V,dd))
#######################


	t = dt
	udot  = 2./dt*(u2-u1) - udotb 

	udd = 1./dt0*(udot - udotb)




	a12 = -inner(div(v),p)*dx
	a21 = dt*inner(div(d),q)*dx


	while (t <= T + DOLFIN_EPS):   
		if step > 5000:
			reason = 'too many steps'
			break

		
	
	##########################	
	# wind vector
		ww  = (1.0+ dt/dt0)*u2 - dt/dt0*u1
 
	#2.11 discrete Oseen problem

		at = 2.0*inner(d,v)*dx + nu*dt*inner(grad(d),grad(v))*dx + dt*inner(grad(d)*ww,v)*dx			

		F1 = (inner(udotb,v) - nu*inner(grad(u2),grad(v)) - inner(grad(u2)*ww,v))*dx
		F2 = 0.0



		# build block matrices	
		AA = block_assemble([[at, a12],
			     [a21,  0 ]], bcs=bcs)

		bb2  = block_assemble([F1, F2], bcs=bcs)
	
		# assign blocks
		[[At, B],
	 	[C, _]] = AA
		

	
	#################################	
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

	########################
	

	# Create the block inverse
		
		Ainv = BiCGStab(AA, precond=prec, tolerance=1e-8, maxiter=100, show=2)
				#show = 3 -> shows the residual in iteration-> very slow
		tic = time.clock()
		# solve LGS
		dd,pp = Ainv * bb2

	
		toc = time.clock()
		ttictoc.append(toc - tic)
		tt = np.array(ttictoc)
		sio.savemat(path + '/tictoc.mat', {'tt':tt})

	# save residuals in an array and save it as .mat for 
		res = Ainv.residuals
		sio.savemat(path + '/residuals/residuals_%d.mat'%step, {'res':res})

		ud.assign(Function(V,dd))
		p1.assign(Function(Q,pp))

	######################################################
	#ww = (1+dt/dt0)*u2 - (dt/dt0)*u2
	

		w = udot + .5*dt*udd
		udiff.assign( ud - w)

		upTR = u2 + dt*ud


	# local truncation error estimate	
	# u_st is the AB2 velocity solution
	# u2
		u_st 	= u2 + dt/2*((2+dt/dt0)*udot - (dt/dt0)*udotb)
	# usind the standard estimate	
		e_h1  	= 1./(3*(1+dt0/dt))*u2
		e_h2	= 1./(3*(1+dt0/dt))*u_st
		e_h.assign(e_h1 -e_h2)
	
	# local truncation error is estimated by comparing the TR sol and AB2 sol
		eps_h = norm(e_h)
		tol_h = ((1/0.7)**3)*tol
	
	
		if eps_h < tol_h:
		# accepted step
			print "accepted time step\n"
		# nstar
			if (step % nstar == 0):				.
		# smoothing every n* step 
		# stabilized TR - AB2 with periodic averaging
			# smoothing by averaging:
				u1.assign(5*(u2 + u1))		# smooth by averaging
				udotb 	= .5*(udot + udotb)
				dt0 	= .5*(dt + dt0)		# correct the old timestep
				u2.assign(u2 + 0.5*dt*ud)
				udot 	= ud
				t 	= t + .5*dt		# leave dt unchanged

			
			
			else:	
			# regular step
				dt0 	= dt
				t 	= t + dt0
				u1.assign(u2)
				u2.assign(u2 + dt*ud)
				udotb 	= udot		
				udot 	= 2*ud - udot
			
			#save dt in a file as .mat for matlab
				vec = np.array(time_steps)
				sio.savemat(path + '/time_step.mat', {'vec':vec})
	
		
			udd = (udot - udotb)/dt0
	
	# time step has been rejected
		else:
			failedstep +=1
	
		print "failed steps = ", failedstep

	# compute next timestep 
		dt = dt *(tol/eps_h)**(1./3.)
		step +=1

	# add current time step to array
		time_steps.append(dt)
	
	
		print "\ndt = " ,dt


	###############################
    
	   	
	##############
	# naming problem
		u2.rename("u","u")
    		Ftmp = Function(V,u2)
    		Ftmp.rename("Visit Velocity", "Visit Velocity")
    		vis_f << Ftmp
	############

	
		ela_time = t/T * 100
		print "\n \nelapsed time: %.2f in percent with %i steps\n \n" % (ela_time, step)
		if bplot == 1:
			# plot velocity
			uplot.assign(Function(V,u2))
			plot(uplot,title = "Velocity (step %d )" %step, interactive = False)
	
			# plot pressure
			pplot.assign(Function(Q,p1))
			plot(pplot, title = "Pressure")
	
	return reason
