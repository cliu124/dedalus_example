"""
Dedalus script for 2D Rayleigh-Benard convection.

This script uses a Fourier basis in the x direction with periodic boundary
conditions.  The equations are scaled in units of the buoyancy time (Fr = 1).

This script can be ran serially or in parallel, and uses the built-in analysis
framework
folder can be used to merge distributed analysis sets from parallel runs,
and the `plot_2d_series.py` script can be used to plot the snapshots.

To run, merge, and plot using 4 processes, for instance, you could use:
    $ mpiexec -n 4 python3 rayleigh_benard.py
    $ mpiexec -n 4 python3 merge.py snapshots
    $ mpiexec -n 4 python3 plot_2d_series.py snapshots/*.h5


"""

import numpy as np
from mpi4py import MPI
import time

from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.core import operators

import logging
logger = logging.getLogger(__name__)


rand = np.random.RandomState(seed=42)

# Define a function to get back the time-step needed to rescale white noise
def forcingx(deltaT):
    gshape = domain.dist.grid_layout.global_shape(scales=3/2)
    slices = domain.dist.grid_layout.slices(scales=3/2)
    noise = rand.standard_normal(gshape)[slices]
    return noise/np.sqrt(deltaT)

def forcingy(deltaT):
    gshape = domain.dist.grid_layout.global_shape(scales=3/2)
    slices = domain.dist.grid_layout.slices(scales=3/2)
    noise = rand.standard_normal(gshape)[slices]
    return noise/np.sqrt(deltaT)


# Parameters
Lx, Lz = (2*np.pi, 2*np.pi)
Re = 10

# Create bases and domain
x_basis = de.Fourier('x', 128, interval=(0, Lx), dealias=3/2)
z_basis = de.Fourier('z', 128, interval=(0, Lz), dealias=3/2)
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)
x = domain.grid(0)
z = domain.grid(1)

# Define the internal heat forcing function (a constant usually)
#forcing_func = domain.new_field(name='forcing_func')
#forcing_func['g'] = 1.

forcing_func_x = operators.GeneralFunction(domain,'g',forcingx,args=[])
forcing_func_y = operators.GeneralFunction(domain,'g',forcingy,args=[])

# 2D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['u','v','p'])
problem.parameters['R'] = Re
problem.parameters['forcing_func_x'] = forcing_func_x
problem.parameters['forcing_func_y'] = forcing_func_y

problem.add_equation("dx(u) + dz(v) = 0", condition="(nx!=0) or (nz!=0)")
problem.add_equation("dt(u) + dx(p) - R**(-1)*(dx(dx(u)) + dz(dz(u)))  = - u*dx(u) - v*dz(u) + forcing_func_x")
problem.add_equation("dt(v) + dx(p) - R**(-1)*(dx(dx(v)) + dz(dz(v)))  = - u*dx(v) - v*dz(v) + forcing_func_y")
problem.add_equation("p=0",condition = "(nx==0) and (nz==0)")

# Build solver
solver = problem.build_solver(de.timesteppers.RK111)
logger.info('Solver built')
logger.info('dt')

#forcing_func.args = [solver.dt]
forcing_func_x.original_args = [0.0001]
forcing_func_y.original_args = [0.0001]
# Initial conditions
x = domain.grid(0)
z = domain.grid(1)

# Random perturbations, initialized globally for same results in parallel
gshape = domain.dist.grid_layout.global_shape(scales=1)
slices = domain.dist.grid_layout.slices(scales=1)
rand = np.random.RandomState(seed=42)
noise = rand.standard_normal(gshape)[slices]

# Linear background + perturbations damped at walls
zb, zt = z_basis.interval
pert =  1e-3 * noise * (zt - z) * (z - zb)

# Integration parameters
solver.stop_sim_time = 1000
solver.stop_wall_time = 10 * 60.
solver.stop_iteration = np.inf

# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots_4e4', sim_dt=0.1, max_writes=50)
snapshots.add_system(solver.state)

# Scalar Data
analysis1 = solver.evaluator.add_file_handler("scalar_data_4e4", iter=10)
analysis1.add_task("integ(0.5*(u*u+v*v))", name="Ek")
analysis1.add_task("integ(0.5*(u*u))", name="Ekx")
analysis1.add_task("integ(0.5*(v*v))", name="Ekz")

analysis1.add_task("z", name="z")     #try to add z for Tz profile graph


analysis1.add_task("R")       #try to add Ra in graph


# horizontally averaged profiles
#analysis2 = solver.evaluator.add_file_handler("profile_data_4e4", iter=100)
#analysis2.add_task("integ(u,'x')", name="u_profile")
#analysis2.add_task("integ(w,'x')", name="w_profile")
#analysis2.add_task("integ(T,'x')", name="T_profile")

#analysis2.add_task("z", name="z")     #try to add z for Tz profile graph

# CFL
CFL = flow_tools.CFL(solver, initial_dt=0.0001, cadence=10, safety=1,
                     max_change=1.5, min_change=0.5, max_dt=0.1)
CFL.add_velocities(('u', 'v'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("sqrt(u*u + v*v) / R", name='Re_rms')

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.ok:
#    for i in range(10):
        dt = CFL.compute_dt()
        forcing_func_x.args = [dt]
        forcing_func_y.args = [dt]
        solver.step(dt)
        if (solver.iteration-1) % 2 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('Max Re = %f' %flow.max('Re_rms'))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
