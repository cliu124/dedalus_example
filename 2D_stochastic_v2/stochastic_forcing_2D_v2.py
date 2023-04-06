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

class flag:
    pass
# Parameters
flag=flag()
flag.Lx, flag.Ly = (2*np.pi, 2*np.pi)
flag.nu  = 0.1
flag.eps = 0.1
flag.Nx=128
flag.Ny=128
flag.k1=7
flag.k2=9
flag.stop_sim_time=10
flag.initial_dt=0.0001

# Create bases and domain
x_basis = de.Fourier('x', flag.Nx, interval=(0, flag.Lx), dealias=3/2)
y_basis = de.Fourier('y', flag.Ny, interval=(0, flag.Ly), dealias=3/2)
domain = de.Domain([x_basis, y_basis], grid_dtype=np.float64)
x = domain.grid(0)
y = domain.grid(1)

tmp_grid_x=domain.new_field()
tmp_grid_x.set_scales(3/2)
#tmp_grid_x_filter=domain.new_field()
#tmp_grid_x_filter.set_scales(3/2)

tmp_grid_y=domain.new_field()
tmp_grid_y.set_scales(3/2)
#tmp_grid_y_filter=domain.new_field()
#tmp_grid_y_filter.set_scales(3/2)


kx=domain.elements(0)
ky=domain.elements(1)
rand = np.random.RandomState(seed=42)
mask=np.invert(kx**2+ky**2>=flag.k1**2)*(kx**2+ky**2<=flag.k2**2)
   

def forcingx(deltaT):
    gshape = domain.dist.grid_layout.global_shape(scales=3/2)
    slices = domain.dist.grid_layout.slices(scales=3/2)
 
    noise = rand.standard_normal(gshape)[slices]
    tmp_grid_x['g']=noise
    tmp_grid_x['c'][mask]=0j
    noise_filter=tmp_grid_x['g']
    tmpx  = 2*np.mean(noise_filter**2)
    noise_filter_normalized = noise_filter*np.sqrt(2*flag.eps)/np.sqrt(tmpx)/np.sqrt(deltaT)
    return noise_filter_normalized

def forcingy(deltaT):
    gshape = domain.dist.grid_layout.global_shape(scales=3/2)
    slices = domain.dist.grid_layout.slices(scales=3/2)
    noise = rand.standard_normal(gshape)[slices]
    tmp_grid_y['g']=noise
    tmp_grid_y['c'][mask]=0j
    noise_filter=tmp_grid_y['g']
    tmpy  = 2*np.mean(noise_filter**2)
    noise_filter_normalized = noise_filter*np.sqrt(2*flag.eps)/np.sqrt(tmpy)/np.sqrt(deltaT)
    return noise_filter_normalized

forcing_func_x = operators.GeneralFunction(domain,'g',forcingx,args=[])
forcing_func_y = operators.GeneralFunction(domain,'g',forcingy,args=[])

#governing equations
problem = de.IVP(domain, variables=['u','v','p','forcing_var_x','forcing_var_y'])
problem.parameters['nu'] = flag.nu
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

#forcing_func.args = [solver.dt]
forcing_func_x.original_args = flag.initial_dt
forcing_func_y.original_args = flag.initial_dt

# Integration parameters
solver.stop_sim_time = flag.stop_sim_time
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf

# Analysis
analysis = solver.evaluator.add_file_handler('analysis', sim_dt=0.1)
analysis.add_system(solver.state)
analysis.add_task("forcing_var_x",layout='c',name='forcing_var_x_coeff')
analysis.add_task("forcing_var_y",layout='c',name='forcing_var_y_coeff')
analysis.add_task("u",layout='c',name='u_coeff')
analysis.add_task("v",layout='c',name='v_coeff')

# Scalar Data
scalar_data = solver.evaluator.add_file_handler("scalar_data", sim_dt=0.1)
scalar_data.add_task("integ(0.5*(u*u+v*v))", name="Ek")
scalar_data.add_task("integ(0.5*(u*u))", name="Ekx")
scalar_data.add_task("integ(0.5*(v*v))", name="Eky")
scalar_data.add_task("(integ(v*v)-integ(u*u))/integ(u*u+v*v)",name="m")

# CFL
CFL = flow_tools.CFL(solver, initial_dt=flag.initial_dt, cadence=10, safety=1,
                     max_change=1.5, min_change=0.5, max_dt=0.1)
CFL.add_velocities(('u', 'v'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("integ(sqrt(u*u + v*v))", name='U_rms')
flow.add_property("(integ( v*v)-integ(u*u))/integ(u*u+v*v)", name='m')

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
            