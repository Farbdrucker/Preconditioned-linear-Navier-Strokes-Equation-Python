from dolfin import*
import numpy as np
import os
from time import*

def bf_bc(mesh,V,Q):
	bmarg = 1.e-3 + DOLFIN_EPS
	# define boundarys
	# Note in- and out flow boundary definition are switched!!
	# Inflow boundary
	class InflowBoundary(SubDomain):
    		def inside(self, x, on_boundary):
        		return on_boundary and x[0] > 8 - bmarg

	# No-slip boundary
	class NoslipBoundary(SubDomain):
    		def inside(self, x, on_boundary):
        		return on_boundary \
				and (x[1] < 0 + bmarg \
				or x[1] > 2 - bmarg \
				or (x[0] < 1. + bmarg and x[1] < 1. + bmarg) \
				or (x[1] < 1. + bmarg and x[0] < 1. + bmarg))

	# Outflow boundary
	class OutflowBoundary(SubDomain):
    		def inside(self, x, on_boundary):
        		return on_boundary and x[0] < 0 + bmarg


	# Create no-slip boundary condition
	noslipboundary  = NoslipBoundary()
	g1  	= Constant((0, 0))
	noslip 	= DirichletBC(V, g1, noslipboundary)	#no slip

	# Create inflow boundary condition for velocity with parabolic inflow
	g0 = Expression(("-amp*(1.0 - x[1])*(2.0 - x[1])", "0.0"),amp = 0.25)
	inflowboundary = OutflowBoundary()
	inflow_u = DirichletBC(V, g0, inflowboundary)


	# Create outflow boundary condition for pressure
	outflowboundary = InflowBoundary()
	g2  		= Constant(0.0)
	outflow_p 	= DirichletBC(Q, g2, outflowboundary)	#outflow

	# collect bcs
	bcu = [inflow_u,noslip]
	bcp = [outflow_p]

	bcs = [bcu, bcp]

	return bcs


def cy_bc(mesh,V,Q):
	# Constants related to the geometry of the mesh
	bmarg   = 1.e-3 + DOLFIN_EPS
	xmin    = 0.0
	xmax    = 2.2
	ymin    = 0.0
	ymax    = 0.41
	xcenter = 0.2
	ycenter = 0.2
	radius  = 0.05

	
	# Inflow boundary
	class InflowBoundary(SubDomain):
    		def inside(self, x, on_boundary):
        		return on_boundary and x[0] < xmin + bmarg

	# No-slip boundary
	class NoslipBoundary(SubDomain):
    		def inside(self, x, on_boundary):
        		dx = x[0] - xcenter
        		dy = x[1] - ycenter
        		r = sqrt(dx*dx + dy*dy)
        		return on_boundary and \
               			(x[1] < ymin + bmarg or x[1] > ymax - bmarg or \
                		r < radius + bmarg)

	# Outflow boundary
	class OutflowBoundary(SubDomain):
    		def inside(self, x, on_boundary):
        		return on_boundary and x[0] > xmax - bmarg

	# Create no-slip boundary condition
	noslipboundary  = NoslipBoundary()
	g1  	= Constant((0, 0))
	noslip 	= DirichletBC(V, g1, noslipboundary)	#no slip

	# Create inflow velocity
	g0 = Expression(("4.*Um*(x[1]*(ymax-x[1]))/(ymax*ymax)", "0.0"),Um = 1.5, ymax=ymax, t=0.0)
	inflowboundary = InflowBoundary()
	inflow_u = DirichletBC(V, g0, inflowboundary)

	#create outflow velocity (same as inlet)
	outflowboundary = OutflowBoundary()
	outflow_u 	= DirichletBC(V, g0, outflowboundary)


	# Create outflow boundary condition for pressure
	outflowboundary = OutflowBoundary()
	g2  		= Constant(0.0)
	outflow_p 	= DirichletBC(Q, g2, outflowboundary)	#outflow
 
	bcu = [inflow_u, noslip]
	bcp = [outflow_p]

	bcs = [bcu, bcp]

	return bcs





def saveall(method,prec,rfnlvl,nu,T,dt,tol = 'default'):
    	#create path
	lt = localtime()
	path =  "/scratch/results/ml_vs_ilu/"+method+"/"+prec+"/rfnlvl_%d/nu_%.3f/T_%d/dt_%.8f"%(rfnlvl,nu,T,dt)

	if method == 'TS' or method == 'TS bf':
		path += "/tol_%.4f"%tol



    	# if folder structure does not exist yet make direction
	if not os.path.exists(path):
        	os.makedirs(path)
    
    	# create a file with all data in in
	f = open(path + "/data.txt","w")
	f.write("method " + method + "\nrfnlvl %d\nnu %f\ndt %f\n" %(rfnlvl,nu,dt))
	if method == 'TS' or method == 'TS bf':
		f.write("tol %f" %tol)
	f.close() 
	
	return path
