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
flag.Lz = (1.) #domain size
flag.Rayleigh = 1.3*10**7 #Rayleigh number
flag.Prandtl = 1
flag.Nz=128 #grid point number in z
flag.kx=20*2*np.pi/5
flag.ky=0
flag.Q=1e6
flag.A_noise=0.1

#a parameter determine the boundary condition, kappa=0 is Dirichlet, and kappa=1 for Neumann
#The top layer boundary condition reads as (1-kappa)*T(z=1)+kappa dT/dz(z=1)=0

#parameter to control simulation and storage time
flag.initial_dt=0.001 #the initial time step
flag.stop_sim_time=100 #The simulation time to stop
flag.post_store_dt=1 #The time step to store the data

# Create bases and domain
z_basis = de.Chebyshev('z', flag.Nz, interval=(0, flag.Lz), dealias=3/2)
domain = de.Domain([z_basis], grid_dtype=np.complex128)

# 2D Boussinesq hydrodynamics
conj = lambda A: np.conj(A)
problem = de.IVP(domain, variables=['u','v','w','p','Jx','Jy','Jz','phi','T','T0','U0','V0', \
                                    'uz','vz','wz','Tz','U0z','V0z','T0z','phi_z'])
    
problem.parameters['Ra'] = flag.Rayleigh
problem.parameters['kappa'] = (flag.Rayleigh * flag.Prandtl)**(-1/2)
problem.parameters['nu'] = (flag.Rayleigh / flag.Prandtl)**(-1/2)
problem.parameters['zi'] = 1j
problem.parameters['kx'] = flag.kx
problem.parameters['ky'] = flag.ky
problem.parameters['Q'] = flag.Q

problem.add_equation("uz - dz(u) = 0")
problem.add_equation("vz - dz(v) = 0")
problem.add_equation("wz - dz(w) = 0")
problem.add_equation("Tz - dz(T) = 0")
problem.add_equation("U0z - dz(U0) = 0")
problem.add_equation("V0z - dz(V0) = 0")
problem.add_equation("T0z - dz(T0) = 0")
problem.add_equation("phi_z-dz(phi) = 0")

problem.add_equation("dt(u)+zi*kx*p-nu*(dz(uz)-kx*kx*u-ky*ky*u)-Q*nu*Jy=-zi*kx*U0*u-zi*ky*V0*u-w*U0z")
problem.add_equation("dt(v)+zi*ky*p-nu*(dz(vz)-kx*kx*v-ky*ky*v)+Q*nu*Jx=-zi*kx*U0*v-zi*ky*V0*v-w*V0z")
problem.add_equation("dt(w)+dz(p)-nu*(dz(wz)-kx*kx*w-ky*ky*w)-T=-zi*kx*U0*w-zi*ky*V0*w")
problem.add_equation("zi*kx*u+zi*ky*v+wz=0")

problem.add_equation("Jx+zi*kx*phi-v=0")
problem.add_equation("Jy+zi*ky*phi+u=0")
problem.add_equation("Jz+phi_z=0")
problem.add_equation("zi*kx*Jx+zi*ky*Jy+dz(Jz)=0")

problem.add_equation("dt(T)-w-kappa*(dz(Tz)-kx*kx*T-ky*ky*T)=-zi*kx*U0*T-zi*ky*V0*T-w*T0z")
problem.add_equation("dt(T0)-kappa*dz(T0z)=-dz(conj(w)*T+conj(T)*w)")
problem.add_equation("dt(U0)-nu*dz(U0z)=-dz(conj(w)*u+conj(u)*w)")
problem.add_equation("dt(V0)-nu*dz(V0z)=-dz(conj(w)*v+conj(v)*w)")

#add B.C.
problem.add_bc("w(z='left')=0")# No penetration
problem.add_bc("w(z='right')=0")
problem.add_bc("uz(z='left')=0")# Stress free
problem.add_bc("uz(z='right')=0")
problem.add_bc("vz(z='left')=0")# Stress free
problem.add_bc("vz(z='right')=0")

problem.add_bc("U0z(z='left')=0")#stress free
problem.add_bc("U0z(z='right')=0")
problem.add_bc("V0z(z='left')=0")
problem.add_bc("V0z(z='right')=0")

problem.add_bc("phi_z(z='left')=0")
problem.add_bc("phi_z(z='right')=0")

#B.C. for temperature
problem.add_bc("T(z='left')=0")
problem.add_bc("T(z='right')=0")
problem.add_bc("T0(z='left')=0")
problem.add_bc("T0(z='right')=0")

# Build solver
solver = problem.build_solver(de.timesteppers.RK222)
logger.info('Solver built')

if not pathlib.Path('restart.h5').exists():

    print('Set up initial condition!')
    # Initial conditions
    z = domain.grid(0)
    
    T = solver.state['T']
    Tz = solver.state['Tz']
    w = solver.state['w']
    wz = solver.state['wz']
    u = solver.state['u']
    T0z = solver.state['T0z']
    # Random perturbations, initialized globally for same results in parallel
    gshape = domain.dist.grid_layout.global_shape(scales=1)
    slices = domain.dist.grid_layout.slices(scales=1)
    rand = np.random.RandomState(seed=42)
    noise = rand.standard_normal(gshape)[slices]

    # Linear background + perturbations damped at walls
    #zb, zt = z_basis.interval
    pert = flag.A_noise * noise * z * (1 - z)
    
    T['g'] = np.sin(np.pi*z) + pert
    
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
solver.stop_sim_time = flag.stop_sim_time

# Analysis
analysis = solver.evaluator.add_file_handler('analysis', sim_dt=flag.post_store_dt)
analysis.add_system(solver.state)

# # CFL
# CFL = flow_tools.CFL(solver, initial_dt=flag.initial_dt, cadence=10, safety=0.5,
#                      max_change=1.5, min_change=0.5, max_dt=0.125, threshold=0.05)
# CFL.add_velocities(('u', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
flow.add_property("integ(sqrt(u*conj(u) +v*conj(v)+ w*conj(w)))/2", name='TKE')           

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
        # dt = CFL.compute_dt()
        dt = solver.step(dt)
        if (solver.iteration-1) % 1000 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e, TKE = %f, Nu-1=%f' %(solver.iteration, solver.sim_time, dt, flow.max('TKE'),-T0z['g'][0]))
            # logger.info('TKE = %f' %flow.max('TKE'))

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