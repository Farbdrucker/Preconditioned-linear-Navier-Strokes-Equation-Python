# import python script
import ns_sa
import ns_euler
import ns_bf_ts

dt = 1
reason = ''
#########################
count = 0	
# boolean if you want to show the plot of velocity, pressure or mesh
bplot = 0

py = ['ns_sa.py','ns_euler.py','ns_bf_sa2.py','ns_bf_ts.py']
	
# main programm to execute ns_***.py severel times with different parameters

# meshes - Cylinder mesh with refinement level from 0 to 8
cy_meshs = [0,1,2,3,4,5,6,7, 8]

# viscosity - Cylinder viscosity from 1:/200 to 1./2000
#cy_nu = [0.001]
cy_nu = [0.005000, 0.003333, 0.002500, 0.001667, 0.001250, 0.001000]


# preconditioning strategies
# LSC and simple available
cy_prec = ['LSC']


# simulation time T
cy_T = 5

#############################
bf_meshs = [8]

# viscosity
bf_nu = [0.02, 0.01, 0.005]

# preconditioning strategies
bf_prec = ['LSC']

# time stepping tol
bf_tol = [1.e-2]

# simulation time T
bf_T = 50

#############################
# create a file with all data in in
lt = localtime()
mastertxt = "log_%02d_%02d_(%02d-%02d).txt" %(lt[2], lt[1], lt[3], lt[4])
f = open(mastertxt,"a")
f.write("Auto genereted log file")
f.write("\n\nNavier-Stokes simulation of flow around an obstacle")
f.write("\nParameters:")
f.write("\n\tMeshs: ")
for m in cy_meshs:
	f.write("%d, "%m)

f.write("\n\tViscosity:")
for nu in cy_nu:
	 f.write("%f, "%nu)

f.write("\n\tPreconditioning: ")
for p in cy_prec:
	f.write(p+"," )

f.write("\n\tTime: %d" %cy_T)

f.write("\n\tShow plots: %d" %bplot)

f.write("\n\n---------------------------------------------")
f.close()



for m in cy_meshs:
	for nu in cy_nu:
		for p in cy_prec:
			# Simo-Armero method	
			writetxt(count,mastertxt, py[0] ,m,nu,p,cy_T)
			ns_sa.main(m,nu,p,bplot,cy_T)			
			count = count+1
			writeend(mastertxt)

		

			# Lagged Euler method
			writetxt(count,mastertxt, py[1] ,m,nu,p,cy_T)
			ns_euler.main(m,nu,p,bplot,cy_T)
			count = count+1
			writeend(mastertxt)
			


########################
# backwardfacing step


f = open(mastertxt,"a")
f.write("\n\n\n\n#############################################")
f.write("\n#############################################")
f.write("\n#############################################")
f.write("\n#############################################")

f.write("\n\nNavier-Stokes simulation of backward facing step")
f.write("\nParameters:")
f.write("\n\tMeshs: ")
for m in bf_meshs:
	f.write("%d, "%m)

f.write("\n\tViscosity:")
for nu in bf_nu:
	 f.write("%f, "%nu)

f.write("\n\tPreconditioning: ")
for p in bf_prec:
	f.write(p+"," )

f.write("\n\ttol:")
for t in bf_tol:
	 f.write("%f, "%t)

f.write("\n\tTime: %d" %bf_T)

f.write("\n\tShow plots: %d" %bplot)

f.write("\n\n---------------------------------------------")
f.close()



for m in bf_meshs:
	for nu in bf_nu:
		for p in bf_prec:
			# bf with Simo-Armero	
			#writetxt(count,mastertxt, py[2] ,m,nu,p,bf_T)
			#ns_sa_bf3.main(m,nu,p,bplot,bf_T)
			#count = count +1
			#f = open(mastertxt,'a')
			#f.write("\ndt =\t\t\t%f"%dt)
			#f.close()
			#writeend(mastertxt)
			

			# bf with adaptive timestepping
			for t in bf_tol:
				
				writetxt(count,mastertxt, py[3] ,m,nu,p,bf_T,t)
				ns_bf_ts.main(m,nu,p,bplot,t,bf_T)
				count = count+1
				writeend(mastertxt)

f = open(mastertxt,'a')

f.write("\n\nENDE\n\n")
f.close()
