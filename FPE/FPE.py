"""
Dedalus script for Fokker-Planck equations.

This script uses a Hermite basis in the x direction with vanishing boundary
conditions.  

This script can be ran serially only because this is a 1D problem. and uses the built-in analysis
framework
"""

import numpy as np
from scipy.ndimage import gaussian_filter
from mpi4py import MPI
import time

from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.core import operators

import logging
logger = logging.getLogger(__name__)

class flag:
    pass
# Parameters
flag=flag()
flag.Lx = 5
flag.Nx=128
flag.stop_sim_time=100
flag.initial_dt=0.01
flag.post_store_dt=0.1

# Create bases and domain
x_basis = de.Hermite('x', flag.Nx, interval=(-flag.Lx, flag.Lx), dealias=3/2)
domain = de.Domain([x_basis], grid_dtype=np.float64)
x = domain.grid(0)

#governing equations
problem = de.IVP(domain, variables=['q'])

problem.add_equation("dt(q) = 1/(np.floor(t)+21)*dx((2*x+0.3*x**2+0.04*x**3)*q)+2/(np.floor(t)+21)**2*dx(qx)")
problem.add_equation("qx-dx(q)=0")

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



#gshape = domain.dist.grid_layout.global_shape(scales=3/2)
#slices = domain.dist.grid_layout.slices(scales=3/2)
#tmp    = rand.standard_normal(gshape)[slices]
#mask   = np.ones_like(tmp['c'])
# mask = np.ones([round(Nx/2),round(Ny-1)],dtype=bool)
# for i in range(len(mask[:,0])):
#     kx=i*2*np.pi/Lx
#     for j in range(len(mask[0,:])):
#         ky=(j-np.floor((Ny-1)/2))*2*np.pi/Ly
#         if ((k1 <= np.sqrt(kx**2 + ky**2) ) and (np.sqrt(kx**2 + ky**2)<=k2)):
#             mask[i,j] = False
            