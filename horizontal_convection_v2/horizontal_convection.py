"""
Dedalus script for 2D Rayleigh-Benard convection.

This script uses a Fourier basis in the x direction with periodic boundary
conditions.  The equations are scaled in units of the buoyancy time (Fr = 1).

This script can be ran serially or in parallel, and uses the built-in analysis
framework to save data snapshots in HDF5 files.  The `merge_procs` command can
be used to merge distributed analysis sets from parallel runs, and the
`plot_slices.py` script can be used to plot the snapshots.

To run, merge, and plot using 4 processes, for instance, you could use:
    $ mpiexec -n 4 python3 inclined_porous_convection.py
    $ mpiexec -n 4 python3 -m dedalus merge_procs analysis --cleanup
    $ mpiexec -n 4 python3 plot_slices.py analysis/*.h5

This script can restart the simulation from the last save of the original
output to extend the integration.  This requires that the output files from
the original simulation are merged, and the last is symlinked or copied to
`restart.h5`.

To run the original example and the restart, you could use:
    $ mpiexec -n 4 python3 inclined_porous_convection.py
    $ mpiexec -n 4 python3 -m dedalus merge_procs analysis
    $ ln -s analysis/analysis_s1.h5 restart.h5
    $ mpiexec -n 4 python3 inclined_porous_convection.py

The simulations should take a few process-minutes to run.

"""

import numpy as np
from mpi4py import MPI
import time
import pathlib
from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.tools import post

import logging
logger = logging.getLogger(__name__)

class flag:
    pass
# Parameters
flag=flag()

# Parameters
flag.Lx, flag.Lz = (10., 1.) #domain size
flag.phi = 35/180*np.pi #inclination angle
flag.Rayleigh = 1e2 #Rayleigh number
flag.Nx=256 #grid point number in x
flag.Nz=64 #grid point number in z

#a parameter determine the boundary condition, kappa=0 is Dirichlet, and kappa=1 for Neumann
#The top layer boundary condition reads as (1-kappa)*T(z=1)+kappa dT/dz(z=1)=0
flag.kappa=0.05

#parameter to control simulation and storage time
flag.initial_dt=0.001 #the initial time step
flag.stop_sim_time=100 #The simulation time to stop
flag.post_store_dt=0.4 #The time step to store the data

#paramter for the initial guess
flag.A_noise=1e-3 #random noise magnitude in the initial condition
flag.A_LS=1 #The magnitude of initial localized structure guess
flag.modulation='gaussian'# The modulation function shape, either 'sin' or 'gaussian'
flag.gaussian_sigma=1 #The sigma parameter in the Gaussian modulation
flag.restart_t0=1 #if 1, the simulation time will start from zero. Otherwise, will continue the previous one 


# Create bases and domain
x_basis = de.Fourier('x', flag.Nx, interval=(0, flag.Lx), dealias=3/2)
z_basis = de.Chebyshev('z', flag.Nz, interval=(0, flag.Lz), dealias=3/2)
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)

# 2D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['p','T','u','w','Tz','wz'])
problem.parameters['Ra'] = flag.Rayleigh
problem.parameters['sin_phi'] = np.sin(flag.phi)
problem.parameters['cos_phi'] = np.cos(flag.phi)
problem.parameters['kappa'] = flag.kappa

#The base state is U=Ra sin(phi)*[(kappa-1)*z-(kappa-1)/2] 
#T=(kappa-1)*z+1. This will satisfies the boundary conditions. 
#All variables here are perturbations around this base state. 
problem.add_equation("dx(u) + wz = 0")
problem.add_equation("dt(T) - (dx(dx(T)) + dz(Tz))+w*(kappa-1)+Ra*sin_phi*((kappa-1)*z-(kappa-1)/2)*dx(T) = -((u)*dx(T) + w*Tz)")
problem.add_equation(" u + dx(p) - Ra*sin_phi*T = 0")
problem.add_equation(" w + dz(p) - Ra*cos_phi*T = 0")
problem.add_equation("Tz - dz(T) = 0")
problem.add_equation("wz - dz(w) = 0")
problem.add_bc("T(z='left') = 0")
problem.add_bc("(1-kappa)*T(z='right')+kappa*Tz(z='right') = 0")
problem.add_bc("w(z='left') = 0")
problem.add_bc("w(z='right') = 0", condition="(nx != 0)")
problem.add_bc("integ(p) = 0", condition="(nx == 0)")

# Build solver
solver = problem.build_solver(de.timesteppers.RK222)
logger.info('Solver built')

# Initial conditions or restart
if not pathlib.Path('restart.h5').exists():
    print('Set up initial condition!')
    # Initial conditions
    x, z = domain.all_grids()
    
    T = solver.state['T']
    Tz = solver.state['Tz']
    w = solver.state['w']
    wz = solver.state['wz']
    u = solver.state['u']
    # Random perturbations, initialized globally for same results in parallel
    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand = np.random.RandomState(seed=42)
    noise = rand.standard_normal(gshape)[slices]

    # Linear background + perturbations damped at walls
    #zb, zt = z_basis.interval
    pert = flag.A_noise * noise * z * (1 - z)
    
    #flag.A_LS=1 gives 3 wave LS...
    #flag.A_LS=5 gives 2 wave LS... 
    #flag.A_LS=3 gives steady convection roll, not LS..
    
    if flag.modulation=='sin':
        T['g'] = flag.A_LS*np.sin(2*np.pi/(2*flag.Lx)*x)*np.sin(2*np.pi/2*x)*(1-z)*z +pert
        #pert
        T.differentiate('z', out=Tz)
        w['g'] = flag.A_LS*np.sin(2*np.pi/(2*flag.Lx)*x)*np.sin(2*np.pi/2*x)*(1-z)*z +pert
        u['g'] = flag.A_LS*np.sin(2*np.pi/(2*flag.Lx)*x)*np.sin(2*np.pi/2*x)*(1-z)*z +pert
    elif flag.modulation == 'gaussian':
        T['g'] = flag.A_LS*np.exp(-(x-flag.Lx/2)**2/2/flag.gaussian_sigma**2)*np.sin(2*np.pi/2*x)*(1-z)*z +pert
        #pert
        T.differentiate('z', out=Tz)
        w['g'] = flag.A_LS*np.exp(-(x-flag.Lx/2)**2/2/flag.gaussian_sigma**2)*np.sin(2*np.pi/2*x)*(1-z)*z +pert
        u['g'] = flag.A_LS*np.exp(-(x-flag.Lx/2)**2/2/flag.gaussian_sigma**2)*np.sin(2*np.pi/2*x)*(1-z)*z +pert
    
    # Timestepping and output
    dt = flag.initial_dt
    stop_sim_time = flag.stop_sim_time
    fh_mode = 'overwrite'

else:
    # Restart
    print('Restart')
    write, last_dt = solver.load_state('restart.h5', -1)

    # Timestepping and output
    dt = last_dt
    stop_sim_time = flag.stop_sim_time
    fh_mode = 'append'
    if flag.restart_t0:
        solver.sim_time=0
        fh_mode='overwrite'

# Integration parameters
solver.stop_sim_time = stop_sim_time

# Analysis
analysis = solver.evaluator.add_file_handler('analysis', sim_dt=flag.post_store_dt)
analysis.add_system(solver.state)

# CFL
CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=0.5,
                     max_change=1.5, min_change=0.5, max_dt=0.125, threshold=0.05)
CFL.add_velocities(('u', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("integ(sqrt(u*u + w*w))/2", name='TKE')
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
            logger.info('TKE = %f' %flow.max('TKE'))

    #add check point, only store one snapshot
    checkpoint=solver.evaluator.add_file_handler('checkpoint')
    checkpoint.add_system(solver.state)
    end_world_time = solver.get_world_time()
    end_wall_time = end_world_time - solver.start_time
    solver.evaluator.evaluate_handlers([checkpoint], timestep = flag.initial_dt, sim_time = solver.sim_time, world_time=end_world_time, wall_time=end_wall_time, iteration=solver.iteration)
    post.merge_process_files('checkpoint',cleanup=True)

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()