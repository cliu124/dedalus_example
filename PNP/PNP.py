import numpy as np
from scipy.ndimage import gaussian_filter
from mpi4py import MPI
import time

from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.core import operators

import logging
logger = logging.getLogger(__name__)
"""
This file contains all the parameters for Nernst Planck Poisson Equation

And here are couple resources from the internet:

https://www.comsol.com/blogs/how-to-model-ion-exchange-membranes-and-donnan-potentials/
https://en.wikipedia.org/wiki/Nernst–Planck_equation

"""

class flag:
    pass
# Parameters
flag=flag()

L_total=25 # unit of mm
L_membrane= 0.1 # unit of mm

flag.Lx = (L_total-L_membrane)/2 #domain size of one side
flag.scaling= flag.Lx/L_membrane #scaling factor for dx
flag.Nx=128
flag.stop_sim_time=100
flag.initial_dt=0.01
flag.post_store_dt=0.1

# Anode Potential: V
flag.E_an = 2.1

# Cathode Potential: V
flag.E_ca = -2.1

# kB (Boltzmann Constant): J/K
flag.kB = 1.3806452e-23

# elemental charge: C
flag.e0 = 1.602176634e-19

# charge of K+ and H+
flag.e_K = 1
flag.e_H = 1

# Diffusion coefficient of K+ in water: m2/s
flag.D_K_water = 1.96e-9
# Diffusion coefficient of K+ in the membrane: m2/s
flag.D_K_membrane = 3.72e-11

# Diffusion coefficient of H+ in water: m2/s
flag.D_H_water = 9.31e-9
# Diffusion coefficient of H+ in the membrane: m2/s
flag.D_H_membrane = 2.86e-10

# water permitivity: F/m or (s4*A2)/(kg*m3)
flag.eps_water = 80.2 * 8.8541878128e-12

# Space charge density:
flag.rho = 3.95

# Faradaic constant: C/mol
flag.F = 9.64853321233100184e4

# Idea gas constant: J/(mol*K)
flag.R = 8.314

# initial K+ concentration: mol
flag.c0_K = 1

# initial H+ concentration: mol
flag.c0_H = 3e-1

# initial room temperature, in unit of K?
flag.T=300

# Create bases and domain
x_basis = de.Chebyshev('x', flag.Nx, interval=(0,flag.Lx), dealias=1)
domain = de.Domain([x_basis], grid_dtype=np.float64)
x = domain.grid(0)

#define variables
problem = de.IVP(domain, variables=['c_K_left','c_H_left','phi_left',\
                                    'd_c_K_left','d_c_H_left','d_phi_left',\
                                    'c_K_mid','c_H_mid','phi_mid',\
                                    'd_c_K_mid','d_c_H_mid','d_phi_mid',\
                                    'c_K_right','c_H_right','phi_right',\
                                    'd_c_K_right','d_c_H_right','d_phi_right'])

#define parameters
problem.parameters['E_an']=flag.E_an
problem.parameters['E_ca']=flag.E_ca
problem.parameters['e_kBT']=flag.e0/flag.kB/flag.T
problem.parameters['e_K']=flag.e_K
problem.parameters['e_H']=flag.e_H
problem.parameters['D_K_water']=flag.D_K_water
problem.parameters['D_K_membrane']=flag.D_K_membrane
problem.parameters['D_H_water']=flag.D_H_water
problem.parameters['D_H_membrane']=flag.D_H_membrane
problem.parameters['eps_water']=flag.eps_water
problem.parameters['rho']=flag.rho
problem.parameters['F']=flag.F
problem.parameters['R']=flag.R
problem.parameters['c0_K']=flag.c0_K
problem.parameters['c0_H']=flag.c0_H
problem.parameters['scaling']=flag.scaling

#left in water
problem.add_equation('dt(c_K_left)-D_K_water*dx(d_c_K_left)=D_K_water*e_K*e_kBT*dx(c_K_left*d_phi_left)')
problem.add_equation('dt(c_H_left)-D_H_water*dx(d_c_H_left)=D_H_water*e_H*e_kBT*dx(c_H_left*d_phi_left)')
problem.add_equation('-eps_water*dx(d_phi_left)=F*(e_K*c_K_left+e_H*c_H_left)')
problem.add_equation('dx(c_K_left)-d_c_K_left=0')
problem.add_equation('dx(c_H_left)-d_c_H_left=0')
problem.add_equation('dx(phi_left)-d_phi_left=0')

#middle in membrane, #need scaling factor for dx
problem.add_equation('dt(c_K_mid)-D_K_membrane*scaling*dx(d_c_K_mid)=D_K_membrane*e_K*e_kBT*scaling*dx(c_K_mid*d_phi_mid)')
problem.add_equation('dt(c_H_mid)-D_H_membrane*scaling*dx(d_c_H_mid)=D_H_membrane*e_H*e_kBT*scaling*dx(c_H_mid*d_phi_mid)')
problem.add_equation('-eps_water*scaling*dx(d_phi_mid)=F*(e_K*c_K_mid+e_H*c_H_mid)')
problem.add_equation('scaling*dx(c_K_mid)-d_c_K_mid=0')
problem.add_equation('scaling*dx(c_H_mid)-d_c_H_mid=0')
problem.add_equation('scaling*dx(phi_mid)-d_phi_mid=0')

#right in water
problem.add_equation('dt(c_K_right)-D_K_water*dx(d_c_K_right)=D_K_water*e_K*e_kBT*dx(c_K_right*d_phi_right)')
problem.add_equation('dt(c_H_right)-D_H_water*dx(d_c_H_right)=D_H_water*e_H*e_kBT*dx(c_H_right*d_phi_right)')
problem.add_equation('-eps_water*dx(d_phi_right)=F*(e_K*c_K_right+e_H*c_H_right)')
problem.add_equation('dx(c_K_right)-d_c_K_right=0')
problem.add_equation('dx(c_H_right)-d_c_H_right=0')
problem.add_equation('dx(phi_right)-d_phi_right=0')

#BC

#left side
problem.add_bc('left(d_c_K_left)=0')
problem.add_bc('left(d_c_H_left)=0')
problem.add_bc('left(phi_left)=E_an')

#interface between left water and membrane, assume all variable and their x derivative are continuous
problem.add_bc('right(c_K_left)-left(c_K_mid)=0')
problem.add_bc('right(c_H_left)-left(c_H_mid)=0')
problem.add_bc('right(phi_left)-left(phi_mid)=0')
problem.add_bc('right(d_c_K_left)-left(d_c_K_mid)=0')
problem.add_bc('right(d_c_H_left)-left(d_c_H_mid)=0')
problem.add_bc('right(d_phi_left)-left(d_phi_mid)=0')

#interface between membrane and right water, assume all variable and their x derivative are continuous
problem.add_bc('right(c_K_mid)-left(c_K_right)=0')
problem.add_bc('right(c_H_mid)-left(c_H_right)=0')
problem.add_bc('right(phi_mid)-left(phi_right)=0')
problem.add_bc('right(d_c_K_mid)-left(d_c_K_right)=0')
problem.add_bc('right(d_c_H_mid)-left(d_c_H_right)=0')
problem.add_bc('right(d_phi_mid)-left(d_phi_right)=0')

#right size
problem.add_bc('right(d_c_K_right)=0')
problem.add_bc('right(d_c_H_right)=0')
problem.add_bc('right(phi_right)=E_ca')

# Build solver
solver = problem.build_solver(de.timesteppers.RK111)
logger.info('Solver built')

# Integration parameters
solver.stop_sim_time = flag.stop_sim_time
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf

# Analysis
analysis = solver.evaluator.add_file_handler('analysis', sim_dt=flag.post_store_dt)
analysis.add_system(solver.state)

def print_screen(flag,logger):
    #print the flag onto the screen
    flag_attrs=vars(flag)
    #print(', '.join("%s: %s, \n" % item for item in flag_attrs.items()))
    logger.info(', Attributes: Value,\n,')
    logger.info(', '.join("%s: %s, \n" % item for item in flag_attrs.items()))

def print_file(flag):
    #print the flag onto file
    flag_text=open('./analysis'+'/flag.txt','w+')
    flag_attrs=vars(flag)
    print(', Attributes: 123,\n ------\n-------\n------',file=flag_text)
    print(', test: 123,',file=flag_text)
    print(', '+', '.join("%s: %s, \n" % item for item in flag_attrs.items()),file=flag_text)
    flag_text.close()
# Main loop

try:
    logger.info('Starting loop')
    start_time = time.time()
    print_screen(flag,logger)
    print_file(flag)
    dt=flag.initial_dt
    while solver.ok:
        solver.step(dt)
        if (solver.iteration-1) % 100 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
