"""
Dedalus script for 2D double diffusive convection with settling velocity.

This script uses a Fourier basis in the x direction with periodic boundary
conditions.  The equations are scaled in units of the buoyancy time (Fr = 1).

This script can be ran serially or in parallel, and uses the built-in analysis
framework to save data snapshots in HDF5 files.  The `merge_procs` command can
be used to merge distributed analysis sets from parallel runs, and the
`plot_slices.py` script can be used to plot the snapshots.

To run, merge, and plot using 4 processes, for instance, you could use:
    $ mpiexec -n 4 python3 double_diffusive_settling.py
    $ mpiexec -n 4 python3 -m dedalus merge_procs snapshots
    $ mpiexec -n 4 python3 plot_slices.py snapshots/*.h5

This script can restart the simulation from the last save of the original
output to extend the integration.  This requires that the output files from
the original simulation are merged, and the last is symlinked or copied to
`restart.h5`.

To run the original example and the restart, you could use:
    $ mpiexec -n 4 python3 double_diffusive_settling.py
    $ mpiexec -n 4 python3 -m dedalus merge_procs snapshots
    $ ln -s snapshots/snapshots_s2.h5 restart.h5
    $ mpiexec -n 4 python3 double_diffusive_settling.py

The simulations should take a few process-minutes to run.

"""

import numpy as np
from mpi4py import MPI
import time
import pathlib

from dedalus import public as de
from dedalus.extras import flow_tools

import logging
logger = logging.getLogger(__name__)

class flag:
    pass
# Parameters
flag=flag()

# Parameters
flag.Lx, flag.Lz = (4., 4.)
flag.Nx, flag.Nz = (256,256)
flag.Prandtl = 1.
#Rayleigh = 1e6
flag.tau=0.01 #diffusivity ratio
flag.R_rho=2 #density ratio
flag.d_T_z=0 #background temperature gradient
flag.d_C_z=1 #background salinity gradient
flag.W_st=0 #non-dimensional settling velocity
flag.initial_dt = 0.125
flag.stop_sim_time = 300
flag.post_store_dt = 0.25
flag.A_noise=1e-3
flag.phi=1
flag.phi_C=0

# Create bases and domain
x_basis = de.Fourier('x', flag.Nx, interval=(0, flag.Lx), dealias=3/2)
z_basis = de.Fourier('z', flag.Nz, interval=(0, flag.Lz), dealias=3/2)
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)

# 2D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['p','T','C','u','w'])
problem.parameters['Pr'] = flag.Prandtl
problem.parameters['tau'] = flag.tau
problem.parameters['R_rho'] = flag.R_rho
problem.parameters['d_T_z'] = flag.d_T_z
problem.parameters['d_C_z'] = flag.d_C_z
problem.parameters['W_st'] = flag.W_st
problem.parameters['phi'] = flag.phi
problem.parameters['phi_C'] = flag.phi_C

problem.add_equation("dx(u) + dz(w) = 0",condition="(nx!=0) or (nz!=0)") #\nabla \cdot \boldsymbol{u}=0
problem.add_equation("p=0", condition="(nx==0) and (nz==0)") #pressure gauge condition
problem.add_equation("dt(T) - (dx(dx(T)) + dz(dz(T)))   +phi*u+d_T_z*w          = -(u*dx(T) + w*dz(T))") #temperature equation,
problem.add_equation("dt(C) - tau*(dx(dx(C)) + dz(dz(C)))-W_st*dz(C) +d_C_z/R_rho*w = -(u*dx(C) + w*dz(C))") #concentration equation
problem.add_equation("dt(u) - Pr*(dx(dx(u)) + dz(dz(u))) + Pr*dx(p)     = -(u*dx(u) + w*dz(u))") #horizontal velocity equation
problem.add_equation("dt(w) - Pr*(dx(dx(w)) + dz(dz(w))) + Pr*dz(p)-(T-C) = -(u*dx(w) + w*dz(w))") #vertical velocity equation

# Build solver
solver = problem.build_solver(de.timesteppers.RK222)
logger.info('Solver built')

# Initial conditions or restart
if not pathlib.Path('restart.h5').exists():

    # Initial conditions
    x, z = domain.all_grids()
    T = solver.state['T']
    C = solver.state['C']
    # Random perturbations, initialized globally for same results in parallel
    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand = np.random.RandomState(seed=42)
    noise = rand.standard_normal(gshape)[slices]

    # Linear background + perturbations damped at walls
    pert =  flag.A_noise * noise 
    T['g'] = pert
    C['g'] = pert

    # Timestepping and output
    dt = flag.initial_dt
    stop_sim_time = flag.stop_sim_time
    fh_mode = 'overwrite'

else:
    # Restart
    write, last_dt = solver.load_state('restart.h5', -1)

    # Timestepping and output
    dt = last_dt
    stop_sim_time = flag.stop_sim_time
    fh_mode = 'append'

# Integration parameters
solver.stop_sim_time = stop_sim_time

# Analysis
analysis = solver.evaluator.add_file_handler('analysis', sim_dt=flag.post_store_dt)
analysis.add_system(solver.state)
analysis.add_task("C",layout='c',name='C_coeff')
analysis.add_task("T",layout='c',name='T_coeff')
analysis.add_task("u",layout='c',name='u_coeff')
analysis.add_task("w",layout='c',name='w_coeff')

# CFL
CFL = flow_tools.CFL(solver, initial_dt=flag.initial_dt, cadence=10, safety=0.5,
                     max_change=1.5, min_change=0.5, max_dt=0.125, threshold=0.05)
CFL.add_velocities(('u', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("sqrt(u*u + w*w)/2", name='KE')
#flow_out = flow_tools.GlobalFlowProperty(solver, cadence=1)
#flow_out.add_property('w*b',name='wb')
               

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
    print_screen(flag,logger)
    print_file(flag)
    while solver.proceed:
        dt = CFL.compute_dt()
        dt = solver.step(dt)
        if (solver.iteration-1) % 10 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('kinetic energy = %f' %flow.max('KE'))
            #dy_T_mean_q=flow_out.volume_average('wb')-1                
            #logger.info('dy_T_mean_q: {}'.format(dy_T_mean_q))
            #logger.info('Nu: {}'.format(-1/dy_T_mean_q))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()