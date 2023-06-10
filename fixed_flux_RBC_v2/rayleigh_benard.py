"""
Dedalus script for 2D Rayleigh-Benard convection.

This script uses a Fourier basis in the x direction with periodic boundary
conditions.  The equations are scaled in units of the buoyancy time (Fr = 1).

This script can be ran serially or in parallel, and uses the built-in analysis
framework to save data snapshots in HDF5 files.  The `merge_procs` command can
be used to merge distributed analysis sets from parallel runs, and the
`plot_slices.py` script can be used to plot the snapshots.

To run, merge, and plot using 4 processes, for instance, you could use:
    $ mpiexec -n 4 python3 rayleigh_benard.py
    $ mpiexec -n 4 python3 -m dedalus merge_procs snapshots
    $ mpiexec -n 4 python3 plot_slices.py snapshots/*.h5

This script can restart the simulation from the last save of the original
output to extend the integration.  This requires that the output files from
the original simulation are merged, and the last is symlinked or copied to
`restart.h5`.

To run the original example and the restart, you could use:
    $ mpiexec -n 4 python3 rayleigh_benard.py
    $ mpiexec -n 4 python3 -m dedalus merge_procs snapshots
    $ ln -s snapshots/snapshots_s2.h5 restart.h5
    $ mpiexec -n 4 python3 rayleigh_benard.py

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
flag.Lx, flag.Lz = (4., 1.)
flag.Prandtl = 1.
flag.Rayleigh = 1e8
flag.Nx=512
flag.Nz=128
flag.A_noise=1e-3
flag.initial_dt=0.125
flag.stop_sim_time=1000
flag.post_store_dt=10


# Create bases and domain
x_basis = de.Fourier('x', flag.Nx, interval=(0, flag.Lx), dealias=3/2)
z_basis = de.Chebyshev('z', flag.Nz, interval=(-flag.Lz/2, flag.Lz/2), dealias=3/2)
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)

# 2D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['p','b','u','w','bz','uz','wz'])
problem.parameters['P'] = (flag.Rayleigh * flag.Prandtl)**(-1/2)
problem.parameters['R'] = (flag.Rayleigh / flag.Prandtl)**(-1/2)
problem.parameters['F'] = F = 1
#problem.parameters['Lx'] = Lx
#problem.parameters['Lz'] = Lz
problem.add_equation("dx(u) + wz = 0")
problem.add_equation("dt(b) - P*(dx(dx(b)) + dz(bz))             = -(u*dx(b) + w*bz)")
problem.add_equation("dt(u) - R*(dx(dx(u)) + dz(uz)) + dx(p)     = -(u*dx(u) + w*uz)")
problem.add_equation("dt(w) - R*(dx(dx(w)) + dz(wz)) + dz(p) - b = -(u*dx(w) + w*wz)")
problem.add_equation("bz - dz(b) = 0")
problem.add_equation("uz - dz(u) = 0")
problem.add_equation("wz - dz(w) = 0")
problem.add_bc("bz(z='left') = -1")
problem.add_bc("bz(z='right') = -1")
problem.add_bc("u(z='left') = 0")
problem.add_bc("u(z='right') = 0")
problem.add_bc("w(z='left') = 0")
problem.add_bc("w(z='right') = 0", condition="(nx != 0)")
problem.add_bc("integ(p) = 0", condition="(nx == 0)")

# Build solver
solver = problem.build_solver(de.timesteppers.RK222)
logger.info('Solver built')

# Initial conditions or restart
if not pathlib.Path('restart.h5').exists():

    # Initial conditions
    x, z = domain.all_grids()
    b = solver.state['b']
    bz = solver.state['bz']

    # Random perturbations, initialized globally for same results in parallel
    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand = np.random.RandomState(seed=42)
    noise = rand.standard_normal(gshape)[slices]

    # Linear background + perturbations damped at walls
    zb, zt = z_basis.interval
    pert =  flag.A_noise * noise * (zt - z) * (z - zb)
    b['g'] = F * pert
    b.differentiate('z', out=bz)

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
#analysis.add_system(solver.state)
analysis.add_task("b",layout='g',name='b')
analysis.add_task("w",layout='g',name='w')
analysis.add_task("bz",layout='g',name='bz')

# CFL
CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=0.5,
                     max_change=1.5, min_change=0.5, max_dt=0.125, threshold=0.05)
CFL.add_velocities(('u', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("sqrt(u*u + w*w) / R", name='Re')
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
            logger.info('Max Re = %f' %flow.max('Re'))
            #dy_T_mean_q=flow_out.volume_average('wb')-1                
            #logger.info('dy_T_mean_q: {}'.format(dy_T_mean_q))
            #logger.info('Nu: {}'.format(-1/dy_T_mean_q))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()