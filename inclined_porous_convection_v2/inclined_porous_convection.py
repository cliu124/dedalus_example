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
flag.Lx, flag.Lz = (15., 1.) #domain size
flag.phi = 35/180*np.pi #inclination angle
flag.Rayleigh = 1e2 #Rayleigh number
flag.Nx=384 #grid point number in x
flag.Nz=64 #grid point number in z

#a parameter determine the boundary condition, kappa=0 is Dirichlet, and kappa=1 for Neumann
#The top layer boundary condition reads as (1-kappa)*T(z=1)+kappa dT/dz(z=1)=0
flag.kappa=0

#parameter to control simulation and storage time
flag.initial_dt=0.001 #the initial time step
flag.stop_sim_time=100 #The simulation time to stop
flag.post_store_dt=10 #The time step to store the data

#paramter for the initial guess
flag.A_noise=0 #random noise magnitude in the initial condition
flag.A_LS=0 #The magnitude of initial localized structure guess
flag.modulation='sin'# The modulation function shape, either 'sin' or 'gaussian'
flag.initial_phase=np.pi/2
flag.gaussian_sigma=1 #The sigma parameter in the Gaussian modulation
flag.restart_t0=1 #if 1, the simulation time will start from zero. Otherwise, will continue the previous one 

#collision index
#flag.collision1=flag.collision2=0 will give the same results as before
#flag.collision1=[1,2,3,4,5,6] will specify the pulse number in the left domain. Note that 6 corresponds to conduction state
#flag.collision2=[1,2,3,4,5,6] will specify the pulse number in the right domain. Note that 6 corresponds to conduction state
flag.collision1=0
flag.collision2=0

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
#problem.add_equation("dt(T) - (dx(dx(T)) + dz(Tz))+w*(kappa-1)+Ra*sin_phi*((kappa-1)*z-(kappa-1)/2)*dx(T) = -((u)*dx(T) + w*Tz)")
problem.add_equation("dt(T) - (dx(dx(T)) + dz(Tz))-w+Ra*sin_phi*(1/2-z)*dx(T) = -((u)*dx(T) + w*Tz)")
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
if flag.collision1!=0 and flag.collision2!=0:
    #half horizontal domain, just read the data
    x_basis_half1 = de.Fourier('x', flag.Nx/2, interval=(0, flag.Lx/2), dealias=1)
    
    #ignore below, just repeat building solvers.
    z_basis_half1 = de.Chebyshev('z', flag.Nz, interval=(0, flag.Lz), dealias=1)
    domain_half1 = de.Domain([x_basis_half1, z_basis_half1], grid_dtype=np.float64)
    problem_half1 = de.IVP(domain_half1, variables=['p','T','u','w','Tz','wz'])
    problem_half1.parameters['Ra'] = flag.Rayleigh
    problem_half1.parameters['sin_phi'] = np.sin(flag.phi)
    problem_half1.parameters['cos_phi'] = np.cos(flag.phi)
    problem_half1.parameters['kappa'] = flag.kappa
    problem_half1.add_equation("dx(u) + wz = 0")
    problem_half1.add_equation("dt(T) - (dx(dx(T)) + dz(Tz))-w+Ra*sin_phi*(1/2-z)*dx(T) = -((u)*dx(T) + w*Tz)")
    problem_half1.add_equation(" u + dx(p) - Ra*sin_phi*T = 0")
    problem_half1.add_equation(" w + dz(p) - Ra*cos_phi*T = 0")
    problem_half1.add_equation("Tz - dz(T) = 0")
    problem_half1.add_equation("wz - dz(w) = 0")
    problem_half1.add_bc("T(z='left') = 0")
    problem_half1.add_bc("(1-kappa)*T(z='right')+kappa*Tz(z='right') = 0")
    problem_half1.add_bc("w(z='left') = 0")
    problem_half1.add_bc("w(z='right') = 0", condition="(nx != 0)")
    problem_half1.add_bc("integ(p) = 0", condition="(nx == 0)")
    #ignore above, nothing special

    solver_half1 = problem_half1.build_solver(de.timesteppers.RK222)
    write, last_dt_half1 = solver_half1.load_state('X'+str(flag.collision1)+'_checkpoint_s1.h5',-1) #edit this path if you are reading data from a different path
    
    x_basis2 = de.Fourier('x', flag.Nx/2, interval=(0, flag.Lx/2), dealias=1)
    # #ignore below, just repeat building solvers.
    z_basis2 = de.Chebyshev('z', flag.Nz, interval=(0, flag.Lz), dealias=1)
    domain2 = de.Domain([x_basis2, z_basis2], grid_dtype=np.float64)
    problem2 = de.IVP(domain2, variables=['p','T','u','w','Tz','wz'])
    problem2.parameters['Ra'] = flag.Rayleigh
    problem2.parameters['sin_phi'] = np.sin(flag.phi)
    problem2.parameters['cos_phi'] = np.cos(flag.phi)
    problem2.parameters['kappa'] = flag.kappa
    problem2.add_equation("dx(u) + wz = 0")
    problem2.add_equation("dt(T) - (dx(dx(T)) + dz(Tz))-w+Ra*sin_phi*(1/2-z)*dx(T) = -((u)*dx(T) + w*Tz)")
    problem2.add_equation(" u + dx(p) - Ra*sin_phi*T = 0")
    problem2.add_equation(" w + dz(p) - Ra*cos_phi*T = 0")
    problem2.add_equation("Tz - dz(T) = 0")
    problem2.add_equation("wz - dz(w) = 0")
    problem2.add_bc("T(z='left') = 0")
    problem2.add_bc("(1-kappa)*T(z='right')+kappa*Tz(z='right') = 0")
    problem2.add_bc("w(z='left') = 0")
    problem2.add_bc("w(z='right') = 0", condition="(nx != 0)")
    problem2.add_bc("integ(p) = 0", condition="(nx == 0)")
    # #ignore above, nothing special

    solver_half2 = problem2.build_solver(de.timesteppers.RK222)
    write, last_dt_half2 = solver_half2.load_state('X'+str(flag.collision2)+'_checkpoint_s1.h5',-1) #Edit this path if you are reading data from a different path
    
    solver.state['T']['g']=np.vstack((solver_half1.state['T']['g'],solver_half2.state['T']['g']))
    solver.state['Tz']['g']=np.vstack((solver_half1.state['Tz']['g'],solver_half2.state['Tz']['g']))
    solver.state['w']['g']=np.vstack((solver_half1.state['w']['g'],solver_half2.state['w']['g']))
    solver.state['wz']['g']=np.vstack((solver_half1.state['wz']['g'],solver_half2.state['wz']['g']))
    solver.state['u']['g']=np.vstack((solver_half1.state['u']['g'],solver_half2.state['u']['g']))
    solver.state['p']['g']=np.vstack((solver_half1.state['p']['g'],solver_half2.state['p']['g']))
    
    dt = np.min([last_dt_half1,last_dt_half2])
    stop_sim_time = flag.stop_sim_time
    if flag.restart_t0:
        solver.sim_time=0
        fh_mode='overwrite'
    else: 
        solver.sim_time=np.max([solver_half1.sim_time,solver_half2.sim_time])
        fh_mode = 'append'

   
elif flag.collision1==0 and flag.collision2!=0:
    #the second is zero, so it is not active. This is only for flip the direction of collision 1 state
    solver.load_state('X'+str(flag.collision2)+'_checkpoint_s1.h5',-1)
    # if flag.collision2<0:
    #     #create a temporary variable
        
    #     #flip the T variable and change sign
    #     tmp1=domain.new_field()
    #     tmp1['c']=solver.state['T']['c']
    #     #tmp1.require_layout(domain.dist.layouts[0])
    #     #tmp1.data[1:,:]=tmp1.data[:0:-1,:]
    #     tmp1.require_layout(domain.dist.layouts[1])
    #     tmp1.data[:,:]=np.flipud(tmp1.data[:,::-1])
    #     #tmp1.data=np.flipud(tmp1.data)
    #     tmp1.set_scales(1.5)
    #     solver.state['T']['g']=-tmp1['g']
        
    #     tmp2=domain.new_field()
    #     tmp2['c']=solver.state['Tz']['c']
    #     #tmp2.require_layout(domain.dist.layouts[0])
    #     #tmp2.data[1:,:]=tmp2.data[:0:-1,:]
    #     tmp2.require_layout(domain.dist.layouts[1])
    #     tmp2.data[:,:]=np.flipud(tmp2.data[:,::-1])
    #     ##tmp2.data=np.flipud(tmp2.data)
    #     tmp2.set_scales(1.5)
    #     solver.state['Tz']['g']=tmp2['g']
        
    #     tmp3=domain.new_field()
    #     tmp3['c']=solver.state['w']['c']
    #     #tmp3.require_layout(domain.dist.layouts[0])
    #     #tmp3.data[1:,:]=tmp3.data[:0:-1,:]
    #     tmp3.require_layout(domain.dist.layouts[1])
    #     tmp3.data[:,:]=np.flipud(tmp3.data[:,::-1])
    #     #tmp3.data=np.flipud(tmp3.data)
    #     tmp3.set_scales(1.5)
    #     solver.state['w']['g']=-tmp3['g']
        
    #     tmp4=domain.new_field()
    #     tmp4['c']=solver.state['wz']['c']
    #     #tmp4.require_layout(domain.dist.layouts[0])
    #     #tmp4.data[1:,:]=tmp4.data[:0:-1,:]
    #     tmp4.require_layout(domain.dist.layouts[1])
    #     tmp4.data[:,:]=np.flipud(tmp4.data[:,::-1])
    #     #tmp4.data=np.flipud(tmp4.data)
    #     tmp4.set_scales(1.5)
    #     solver.state['wz']['g']=tmp4['g']
        
    #     tmp5=domain.new_field()
    #     tmp5['c']=solver.state['u']['c']
    #     #tmp5.require_layout(domain.dist.layouts[0])
    #     #tmp5.data[1:,:]=tmp5.data[:0:-1,:]
    #     tmp5.require_layout(domain.dist.layouts[1])
    #     tmp5.data[:,:]=np.flipud(tmp5.data[:,::-1])
    #     #tmp5.data=np.flipud(tmp5.data)
    #     tmp5.set_scales(1.5)
    #     solver.state['u']['g']=-tmp5['g']
        
    #     tmp6=domain.new_field()
    #     tmp6['c']=solver.state['p']['c']
    #     #tmp6.require_layout(domain.dist.layouts[0])
    #     #tmp6.data[1:,:]=tmp6.data[:0:-1,:]
    #     tmp6.require_layout(domain.dist.layouts[1])
    #     tmp6.data[:,:]=np.flipud(tmp6.data[:,::-1])
    #     #tmp6.data=np.flipud(tmp6.data)
    #     tmp6.set_scales(1.5)
    #     solver.state['p']['g']=tmp6['g']
    if flag.restart_t0:
        solver.sim_time=0
        fh_mode='overwrite'
        
elif flag.collision1==0 and flag.collision2==0:
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
            T['g'] = flag.A_LS*np.sin(2*np.pi/(2*flag.Lx)*x+flag.initial_phase)*np.sin(2*np.pi/2*x)*(1-z)*z +pert
            #pert
            T.differentiate('z', out=Tz)
            w['g'] = flag.A_LS*np.sin(2*np.pi/(2*flag.Lx)*x+flag.initial_phase)*np.sin(2*np.pi/2*x)*(1-z)*z +pert
            u['g'] = flag.A_LS*np.sin(2*np.pi/(2*flag.Lx)*x+flag.initial_phase)*np.sin(2*np.pi/2*x)*(1-z)*z +pert
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
solver.stop_sim_time = flag.stop_sim_time

# Analysis
analysis = solver.evaluator.add_file_handler('analysis', sim_dt=flag.post_store_dt)
analysis.add_system(solver.state)

# CFL
CFL = flow_tools.CFL(solver, initial_dt=flag.initial_dt, cadence=10, safety=0.5,
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