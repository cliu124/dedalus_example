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

from dedalus import public as de
from dedalus.extras import flow_tools
from dedalus.core import operators

import logging
logger = logging.getLogger(__name__)

rand = np.random.RandomState(seed=42)

# Parameters
Lx, Ly = (2*np.pi, 2*np.pi)
nu  = 0.01
eps = 1.0
Nx=128
Ny=128

#rand = np.random.RandomState(seed=42)

k1=7
k2=9

# Create bases and domain
x_basis = de.Fourier('x', Nx, interval=(0, Lx), dealias=3/2)
y_basis = de.Fourier('y', Ny, interval=(0, Ly), dealias=3/2)
domain = de.Domain([x_basis, y_basis], grid_dtype=np.float64)
x = domain.grid(0)
y = domain.grid(1)

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
            
# print(mask)
tmp_grid=domain.new_field()
tmp_grid.set_scales(3/2)
tmp_grid_filter=domain.new_field()
tmp_grid_filter.set_scales(3/2)
tmp_grid_filter['g']=0
#mask=bool(mask)
# Define a function to get back the time-step needed to rescale white noise
kx=domain.elements(0)
ky=domain.elements(1)

def forcingx(deltaT):
    gshape = domain.dist.grid_layout.global_shape(scales=3/2)
    slices = domain.dist.grid_layout.slices(scales=3/2)
    noise = rand.standard_normal(gshape)[slices]
    #print(deltaT)
    #print(noise)
    #print(slices)
    #print(len(noise[0,:]))
    #print(len(noise[:,0]))
    tmp_grid['g']=noise
    #tmp_grid_copy = ((tmp_grid + 1e-16) - 1e-16).evaluate()
    #mask_slices=mask[slices]
    #tmp_grid_copy['c'][mask_slices]=0j
    #print((kx**2+ky**2>k1**2)*(kx**2+ky**2<k2**2))
    tmp_grid_filter['c'][(kx**2+ky**2>k1**2)*(kx**2+ky**2<k2**2)]=tmp_grid['c'][(kx**2+ky**2>k1**2)*(kx**2+ky**2<k2**2)]
    noise_filter=tmp_grid_filter['g']
    tmpx  = 2*np.mean(noise_filter**2)
    noise_filter_normalized = noise_filter*np.sqrt(2*eps)/np.sqrt(tmpx)/np.sqrt(deltaT)
    #noise = gaussian_filter(noise, sigma=1)
    return noise_filter_normalized

def forcingy(deltaT):
    gshape = domain.dist.grid_layout.global_shape(scales=3/2)
    slices = domain.dist.grid_layout.slices(scales=3/2)
    noise = rand.standard_normal(gshape)[slices]
    #print(deltaT)
    #print(noise)
    #print(len(noise[0,:]))
    #print(len(noise[:,0]))
    tmp_grid['g']=noise
    #tmp_grid_copy = ((tmp_grid + 1e-16) - 1e-16).evaluate()
    #mask_slices=mask[slices]
    #tmp_grid_copy['c'][mask_slices]=0j
    tmp_grid_filter['c'][(kx**2+ky**2>k1**2)*(kx**2+ky**2<k2**2)]=tmp_grid['c'][(kx**2+ky**2>k1**2)*(kx**2+ky**2<k2**2)]
    noise_filter=tmp_grid_filter['g']
    tmpx  = 2*np.mean(noise_filter**2)
    noise_filter_normalized = noise_filter*np.sqrt(2*eps)/np.sqrt(tmpx)/np.sqrt(deltaT)
    #noise = gaussian_filter(noise, sigma=1)
    return noise_filter_normalized

# Define the internal heat forcing function (a constant usually)
#forcing_func = domain.new_field(name='forcing_func')
#forcing_func['g'] = 1.

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
solver = problem.build_solver(de.timesteppers.RK222)
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
solver.stop_sim_time = 10
solver.stop_wall_time = 10 * 60.
solver.stop_iteration = np.inf

# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=0.1)
snapshots.add_system(solver.state)
snapshots.add_task("forcing_var_x",layout='c',name='forcing_var_x_coeff')
snapshots.add_task("forcing_var_y",layout='c',name='forcing_var_y_coeff')
snapshots.add_task("u",layout='c',name='u_coeff')
snapshots.add_task("v",layout='c',name='v_coeff')

# Scalar Data
analysis1 = solver.evaluator.add_file_handler("scalar_data", sim_dt=0.1)
analysis1.add_task("integ(0.5*(u*u+v*v))", name="Ek")
analysis1.add_task("integ(0.5*(u*u))", name="Ekx")
analysis1.add_task("integ(0.5*(v*v))", name="Eky")
#analysis1.add_task("forcing_var_x",name="forcing_var_x")
#analysis1.add_task("forcing_var_y",name="forcing_var_y")
#analysis1.add_task("z", name="z")     #try to add z for Tz profile graph
#analysis1.add_task("R")       #try to add Ra in graph


# CFL
CFL = flow_tools.CFL(solver, initial_dt=0.0001, cadence=10, safety=1,
                     max_change=1.5, min_change=0.5, max_dt=0.1)
CFL.add_velocities(('u', 'v'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("integ(sqrt(u*u + v*v))", name='U_rms')
flow.add_property("(integ( v*v)-integ(u*u))/integ(u*u+v*v)", name='m')

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
        if (solver.iteration-1) % 100 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
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
