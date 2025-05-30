"""
Dedalus script for 2D Rayleigh-Benard convection.

This script uses a Fourier basis in the x direction with periodic boundary
conditions.  The equations are scaled in units of the buoyancy time (Fr = 1).

This script can be ran serially or in parallel, and uses the built-in analysis
framework
folder can be used to merge distributed analysis sets from parallel runs,
and the `plot_2d_series.py` script can be used to plot the snapshots.

To run, merge, and plot using 4 processes, for instance, you could use:
    $ mpiexec -n 4 python3 stochastic_forcing_2D_v2.py
    $ mpiexec -n 4 python3 merge.py snapshots
    $ mpiexec -n 4 python3 plot_slices.py snapshots/*.h5


"""

import numpy as np
from scipy.ndimage import gaussian_filter
from mpi4py import MPI
import time
import matplotlib.pyplot as plt
from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.core import operators

import logging
logger = logging.getLogger(__name__)

rand = np.random.RandomState(seed=42)

# Parameters
delta = 1.5 # aspect ratio
Lx, Ly = (2*np.pi*delta, 2*np.pi)
nu  = 0.02
eps = 1.0
Nx = 192
Ny = 128

k1 = 7.5
k2 = 8.5

# Create bases and domain
x_basis = de.Fourier('x', Nx, interval=(0, Lx), dealias=3/2)
y_basis = de.Fourier('y', Ny, interval=(0, Ly), dealias=3/2)
domain = de.Domain([x_basis, y_basis], grid_dtype=np.float64)
x = domain.grid(0)
y = domain.grid(1)

tmp_grid=domain.new_field()
tmp_grid.set_scales(3/2)
tmp_grid_filter=domain.new_field()
tmp_grid_filter.set_scales(3/2)
tmp_grid_filter['g']=0

# Define a function to get back the time-step needed to rescale white noise
kx = domain.elements(0); 
ky = domain.elements(1);
print('kx=',kx,'ky=',ky)

def forcingx(deltaT):
    phase = 2*np.pi*np.random.rand(len(kx),len(ky[0]))
    force_range = (kx**2+ky[0]**2>k1**2)*(kx**2+ky[0]**2<k2**2)
    #######################################
    # check graphically  that forcing range is correct
    # KX,KY = np.meshgrid(ky,kx); plt.clf()
    # plt.contourf(KX,KY,force_range); plt.show()
    ######################################
    A = np.cos(phase); B = np.sin(phase)
    for j in range(1,int(len(ky[0])/2)):
       B[0,-j] = - B[0,j]
    tmp_grid_filter['c'][force_range]   =  A[force_range] + 1.0j * B[force_range]
    tmp_grid_filter['c'][1-force_range] =  0.0j
    noise_filter = tmp_grid_filter['g']
    tmpx         = 2*np.mean(noise_filter**2)
    noise_filter = noise_filter*np.sqrt(2*eps)/np.sqrt(tmpx)/np.sqrt(deltaT)
    return noise_filter

def forcingy(deltaT):
    phase = 2*np.pi*np.random.rand(len(kx),len(ky[0]))
    A = np.cos(phase); B = np.sin(phase)
    for j in range(1,int(len(ky[0])/2)):
       B[0,-j] = - B[0,j]    
    force_range = (kx**2+ky[0]**2>k1**2)*(kx**2+ky[0]**2<k2**2)
    tmp_grid_filter['c'][force_range] = A[force_range] + 1.0j * B[force_range]
    tmp_grid_filter['c'][1-force_range] =  0.0j
    noise_filter   = tmp_grid_filter['g']
    tmpy           = 2*np.mean(noise_filter**2)
    noise_filter   = noise_filter*np.sqrt(2*eps)/np.sqrt(tmpy)/np.sqrt(deltaT)
    return noise_filter

forcing_func_x = operators.GeneralFunction(domain,'g',forcingx,args=[])
forcing_func_y = operators.GeneralFunction(domain,'g',forcingy,args=[])

# 2D Boussinesq hydrodynamics
problem = de.IVP(domain, variables=['u','v','p','forcing_var_x','forcing_var_y'])
problem.parameters['nu'] = nu
problem.parameters['forcing_func_x'] = forcing_func_x
problem.parameters['forcing_func_y'] = forcing_func_y

problem.add_equation("dx(u) + dy(v) = 0", condition="(nx!=0) or (ny!=0)")
problem.add_equation("dt(u) + dx(p) - nu*(dx(dx(u)) + dy(dy(u)))  = - u*dx(u) - v*dy(u) + forcing_func_x")
problem.add_equation("dt(v) + dy(p) - nu*(dx(dx(v)) + dy(dy(v)))  = - u*dx(v) - v*dy(v) + forcing_func_y")
problem.add_equation("p=0",condition = "(nx==0) and (ny==0)")
problem.add_equation("forcing_var_x=forcing_func_x")
problem.add_equation("forcing_var_y=forcing_func_y")

# Build solver
solver = problem.build_solver(de.timesteppers.RK111)
logger.info('Solver built')
logger.info('dt')

#forcing_func.args = [solver.dt]
forcing_func_x.original_args = [0.0001]
forcing_func_y.original_args = [0.0001]
# Initial conditions
#x = domain.grid(0)
#z = domain.grid(1)

# Random perturbations, initialized globally for same results in parallel
#gshape = domain.dist.grid_layout.global_shape(scales=1)
#slices = domain.dist.grid_layout.slices(scales=1)
#rand = np.random.RandomState(seed=42)
#noise = rand.standard_normal(gshape)[slices]

# Integration parameters
solver.stop_sim_time = 100
solver.stop_wall_time = 1000 * 60.
solver.stop_iteration = np.inf

# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=0.1)
snapshots.add_system(solver.state)

# Scalar Data
analysis1 = solver.evaluator.add_file_handler("scalar_data", sim_dt=0.01)
analysis1.add_task("integ(0.5*(u*u+v*v))", name="Ek")
analysis1.add_task("integ(0.5*(u*u))", name="Ekx")
analysis1.add_task("integ(0.5*(v*v))", name="Eky")
#analysis1.add_task("forcing_var_x",name="forcing_var_x")
#analysis1.add_task("forcing_var_y",name="forcing_var_y")
#analysis1.add_task("z", name="z")     #try to add z for Tz profile graph
#analysis1.add_task("R")       #try to add Ra in graph


# CFL
CFL = flow_tools.CFL(solver, initial_dt=0.0001, cadence=10, safety=0.5,
                     max_change=1.1, min_change=0.9, max_dt=0.1)
CFL.add_velocities(('u', 'v'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("integ(0.5*(u*u + v*v))", name='E_kin')
flow.add_property("integ(sqrt(u*u + v*v))", name='U_rms')
flow.add_property("(integ(v*v)-integ(u*u))/integ(u*u+v*v)", name='m')

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        dt = CFL.compute_dt()
        forcing_func_x.args = [dt]
        forcing_func_y.args = [dt]
        ###############################
        #CHECK FORCING NORMALIZATION
        ###############################
        # print(np.mean(forcing_func_x['g']**2 + forcing_func_y['g']**2)*(dt/(2*eps)))
        # print(forcing_func_x['g'][0,1])
        solver.step(dt)
        if (solver.iteration-1) % 100 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('E_kin = %f' %flow.max('E_kin'))
            logger.info('U_rms = %f' %flow.max('U_rms'))
            logger.info('m = %f' %flow.max('m'))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
