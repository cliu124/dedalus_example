import numpy as np
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)

#### Parameters ###
Lx, Ly, Lz = (0.6*np.pi, 2.0, 0.18*np.pi)
#Lx, Ly, Lz = (4*np.pi, 2.0, 2*np.pi)
Re = 690 # U_b*H/nu , 16200
#Retau = 180 # = u_tau*H/nu
dtype = np.float64
stop_sim_time = 50

#timestepper, max time step, initial time step, and stop time. 
#timestepper = d3.RK222
timestepper = d3.SBDF4
max_timestep = 1e-5  # 0.125 to 0.01
initial_dt = 1e-5
stop_sim_time = 1000

# Create bases and domain
#nx, ny, nz = 192, 129, 160 #54, 129, 52
nx, ny, nz = 54, 129, 52 #54, 129, 52

geometry='yz' #xy (only streamwise and wall-normal) yz (only wall-normal and spanwise) or xyz(3D)
wavy_wall='spanwise' #'streamwise': streamwise wavy wall; 'spanwise': spanwise wavy wall; 'streamwise_spanwise': 3D wavy wall varying in both streamwise and spanwise
k_inv_scheme='LHS' #RHS: put k_inv term on the RHS of momentum equation, and LHS: put k_inv term on the LHS of the momentum equations
noise_amp_IC=1e-6
solution_method='NLBVP' # or 'NLBVP'
ncc_cutoff=1e-3

if geometry=='xy':
    #coordinates
    coords = d3.CartesianCoordinates('x', 'y')
    dist = d3.Distributor(coords, dtype=np.float64)
    xbasis = d3.RealFourier(coords['x'], size=nx, bounds=(0, Lx), dealias=3/2)
    ybasis = d3.Chebyshev(coords['y'], size=ny, bounds=(-Ly/2, Ly/2), dealias=3/2)
    
    # Fields
    p = dist.Field(name='p', bases=(xbasis,ybasis))
    u = dist.VectorField(coords, name='u', bases=(xbasis,ybasis))
    tau_u1 = dist.VectorField(coords, name='tau_u1', bases=(xbasis))
    tau_u2 = dist.VectorField(coords, name='tau_u2', bases=(xbasis))
    tau_p = dist.Field(name='tau_p')
    
    # Substitutions
    x, y= dist.local_grids(xbasis, ybasis)
    ex, ey = coords.unit_vector_fields(dist)
    lift_basis = ybasis.derivative_basis(1) # Chebyshev U basis
    lift = lambda A: d3.Lift(A, lift_basis, -1) # Shortcut for multiplying by U_{N-1}(y)
    grad_u = d3.grad(u) - ey*lift(tau_u1) # Operator representing G
    x_average = lambda A: d3.Average(A,'x')
    xz_average = lambda A: d3.Average(A, 'x')

elif geometry =='yz':
    #coordinates
    coords = d3.CartesianCoordinates('y','z')
    dist = d3.Distributor(coords, dtype=np.float64)
    ybasis = d3.Chebyshev(coords['y'], size=ny, bounds=(-Ly/2, Ly/2), dealias=3/2)
    zbasis = d3.RealFourier(coords['z'], size=nz, bounds=(0, Lz), dealias=3/2)
    
    # Fields
    #p = dist.Field(name='p', bases=(ybasis,zbasis))
    #for yz plane, we only need to solve U(y,z) and thus just define a velocity field instead of a vector
    u = dist.Field(name='u', bases=(ybasis,zbasis))
    tau_u1 = dist.Field(name='tau_u1', bases=(zbasis))
    tau_u2 = dist.Field(name='tau_u2', bases=(zbasis))
    #tau_p = dist.Field(name='tau_p')
    
    # Substitutions
    y, z = dist.local_grids(ybasis, zbasis)
    ey, ez = coords.unit_vector_fields(dist)
    lift_basis = ybasis.derivative_basis(1) # Chebyshev U basis
    lift = lambda A: d3.Lift(A, lift_basis, -1) # Shortcut for multiplying by U_{N-1}(y)
    grad_u = d3.grad(u) - ey*lift(tau_u1) # Operator representing G
    z_average = lambda A: d3.Average(A,'z')
    xz_average = lambda A: d3.Average(A, 'z')    

elif geometry =='xyz':
    #coordinates
    coords = d3.CartesianCoordinates('x', 'y','z')
    dist = d3.Distributor(coords, dtype=np.float64)
    xbasis = d3.RealFourier(coords['x'], size=nx, bounds=(0, Lx), dealias=3/2)
    ybasis = d3.Chebyshev(coords['y'], size=ny, bounds=(-Ly/2, Ly/2), dealias=3/2)
    zbasis = d3.RealFourier(coords['z'], size=nz, bounds=(0, Lz), dealias=3/2)
    
    # Fields
    p = dist.Field(name='p', bases=(xbasis,ybasis,zbasis))
    u = dist.VectorField(coords, name='u', bases=(xbasis,ybasis,zbasis))
    tau_u1 = dist.VectorField(coords, name='tau_u1', bases=(xbasis,zbasis))
    tau_u2 = dist.VectorField(coords, name='tau_u2', bases=(xbasis,zbasis))
    tau_p = dist.Field(name='tau_p')
    
    # Substitutions
    x, y, z = dist.local_grids(xbasis, ybasis, zbasis)
    ex, ey, ez = coords.unit_vector_fields(dist)
    lift_basis = ybasis.derivative_basis(1) # Chebyshev U basis
    lift = lambda A: d3.Lift(A, lift_basis, -1) # Shortcut for multiplying by U_{N-1}(y)
    grad_u = d3.grad(u) - ey*lift(tau_u1) # Operator representing G
    x_average = lambda A: d3.Average(A,'x')
    xz_average = lambda A: d3.Average(d3.Average(A, 'x'), 'z')

dPdx = -2/Re #Pressure gradient

#parameters for volume penalty
K_inv=100 #690
y0=0.8 #average distance from the channel centerline of the wavy wall. 
A=0.1 #amplitude of wavy wall. 
sharpness=30 #sharpness for the solid fluid boundary, a parameter in tanh function

#
if wavy_wall =='streamwise':
    mask = dist.Field(name='mask', bases=(xbasis, ybasis)) # mask function, only as a function in x and y
    mask['g']= np.tanh(sharpness*(y-(y0+A*np.sin(2*np.pi/Lx*x))))+1-np.tanh(sharpness*(y+y0+A*np.sin(2*np.pi/Lx*x)))

elif wavy_wall =='spanwise':
    mask = dist.Field(name='mask', bases=(ybasis, zbasis)) # mask function, only as a function in y and z
    mask['g']= np.tanh(sharpness*(y-(y0+A*np.sin(2*np.pi/Lz*z))))+1-np.tanh(sharpness*(y+y0+A*np.sin(2*np.pi/Lz*z)))

elif wavy_wall =='streamwise_spanwise':
    mask = dist.Field(name='mask', bases=(xbasis, ybasis, zbasis)) # mask function, is a function of x, y, and z
    mask['g']= np.tanh(sharpness*(y-(y0+A*np.sin(2*np.pi/Lx*x)*np.sin(2*np.pi/Lz*z))))+1-np.tanh(sharpness*(y+y0+A*np.sin(2*np.pi/Lx*x)*np.sin(2*np.pi/Lz*z)))

# Problem
if wavy_wall=='spanwise' and geometry=='yz':
    #This is a scalar equation for U(y,z) in streamwise momentum equation. The nonlinear term, pressure gradient disappear
    #continuity is automatically satiafied and does not need to add. 
    if solution_method=='IVP':
        problem = d3.IVP([u, tau_u1, tau_u2], namespace=locals())
        if k_inv_scheme=='RHS':
            #print('RHS')
            problem.add_equation("dt(u) - 1/Re*div(grad_u) + lift(tau_u2) =-dPdx -K_inv*mask*u")
        elif k_inv_scheme == 'LHS':
            #print('LHS')
            problem.add_equation("dt(u) - 1/Re*div(grad_u) + lift(tau_u2)+K_inv*mask*u =-dPdx")
        #B.C.
        problem.add_equation("u(y=-1) = 0") 
        problem.add_equation("u(y=+1) = 0")
        
        mask.set_scales(1.5)
        # initial condition: Laminar solution + perturbations damped at walls
        np.random.seed(0)
        u['g'] = 0
        #(1-y**2) + np.random.randn(*u['g'].shape) * noise_amp_IC*np.sin(np.pi*(y+1)*0.5) # Laminar solution (plane Poiseuille)+  random perturbation
        
        #In this case, u is a scalar not vector and does not support CFL condition. 
        solver = problem.build_solver(timestepper)
        solver.stop_sim_time = stop_sim_time
        
        # u,tau_u1,tau_u2 = (solver.state[name] for name in problem.variables)
        # for field in [u, tau_u1, tau_u2]: 
        #     field.set_scales(3/2)
        #     field['g'] = 0
        
    elif solution_method=='NLBVP':
        problem = d3.NLBVP([u, tau_u1, tau_u2], namespace=locals())
        if k_inv_scheme=='RHS':
            #print('RHS')
            problem.add_equation("- 1/Re*div(grad_u) + lift(tau_u2) =-dPdx -K_inv*mask*u")
        elif k_inv_scheme == 'LHS':
            #print('LHS')
            problem.add_equation("- 1/Re*div(grad_u) + lift(tau_u2)+K_inv*mask*u =-dPdx")
        #B.C.
        problem.add_equation("u(y=-1) = 0") 
        problem.add_equation("u(y=+1) = 0")
        solver = problem.build_solver(ncc_cutoff=ncc_cutoff)

    # Flow properties
    flow = d3.GlobalFlowProperty(solver, cadence=20) # changed cadence from 10 to 50
    flow.add_property(np.sqrt(u**2)/2, name='TKE')

else:
    problem = d3.IVP([p, u, tau_p, tau_u1, tau_u2], namespace=locals())
    if k_inv_scheme=='RHS':
        problem.add_equation("dt(u) - 1/Re*div(grad_u) + grad(p) + lift(tau_u2) =-dPdx*ex -dot(u,grad(u))-K_inv*mask*u")
    elif k_inv_scheme=='LHS':
        problem.add_equation("dt(u) - 1/Re*div(grad_u) + grad(p) + lift(tau_u2)+K_inv*mask*u =-dPdx*ex -dot(u,grad(u))")
    
    #mass conservation and pressure gauge condition
    problem.add_equation("trace(grad_u) + tau_p = 0")
    problem.add_equation("integ(p) = 0")
    
    #B.C.
    problem.add_equation("u(y=-1) = 0") 
    problem.add_equation("u(y=+1) = 0")
    
    # initial condition: Laminar solution + perturbations damped at walls
    np.random.seed(0)
    u['g'][0] = (1-y**2) + np.random.randn(*u['g'][0].shape) * noise_amp_IC*np.sin(np.pi*(y+1)*0.5) # Laminar solution (plane Poiseuille)+  random perturbation
    
    solver = problem.build_solver(timestepper)
    solver.stop_sim_time = stop_sim_time
    
    CFL = d3.CFL(solver, initial_dt=initial_dt, cadence=5, safety=0.5, threshold=0.05,
                 max_change=1.5, min_change=0.5, max_dt=max_timestep)
    CFL.add_velocity(u) # changed threshold from 0.05 to 0.01
    
    # Flow properties
    flow = d3.GlobalFlowProperty(solver, cadence=20) # changed cadence from 10 to 50
    flow.add_property(np.sqrt(u@u)/2, name='TKE')

# Build Solver
fh_mode = 'overwrite'


snapshots = solver.evaluator.add_file_handler('snapshots_channel', sim_dt=1e-4, max_writes=400)

snapshots.add_task(u, name='velocity')
snapshots.add_task(K_inv*mask, name='stiffness',layout='g')
snapshots.add_task(mask,name="mask",layout="g")

if geometry=='xyz':
    #only output this stress when computing 3D problem.
    snapshots_stress = solver.evaluator.add_file_handler('snapshots_channel_stress', sim_dt=1, max_writes=400)
    snapshots_stress.add_task(xz_average(u),name = 'ubar')
    snapshots_stress.add_task(xz_average(((u-xz_average(u))@ex)**2),name = 'u_prime_u_prime')
    snapshots_stress.add_task(xz_average(((u-xz_average(u))@ey)**2),name = 'v_prime_v_prime')
    snapshots_stress.add_task(xz_average(((u-xz_average(u))@ez)**2),name = 'w_prime_w_prime')
    snapshots_stress.add_task(xz_average(((u-xz_average(u))@ex)*(u-xz_average(u))@ey),name = 'u_prime_v_prime')

# CFL

# Main loop
startup_iter = 10

if wavy_wall=='spanwise' and geometry=='yz':

    if solution_method=='IVP':
        #Fixed time stepper, without using CFL
        try:
            logger.info('Starting main loop')
            while solver.proceed:
                timestep = initial_dt
                solver.step(timestep)
                if (solver.iteration-1) % 10 == 0:
                    max_TKE = flow.max('TKE')
                    logger.info('Iteration=%i, Time=%e, dt=%e, max(TKE)=%f' %(solver.iteration, solver.sim_time, timestep, max_TKE))
        except:
            logger.error('Exception raised, triggering end of main loop.')
            raise
        finally:
            solver.log_stats()
    elif solution_method=='NLBVP':
        pert_norm = np.inf
        #u.change_scales(3/2) #dealising factor
        #steps = [u['g'].ravel().copy()]
        tolerance = 1e-10
        while pert_norm > tolerance:
            solver.newton_iteration()
            pert_norm = sum(pert.allreduce_data_norm('c', 2) for pert in solver.perturbations)
            logger.info(f'Perturbation norm: {pert_norm:.3e}')
            #max_TKE = flow.max('TKE')
            #logger.info('Iteration=%i, max(TKE)=%f' %(solver.iteration, max_TKE))
            logger.info('Iteration=%i' %(solver.iteration))

        
else: 
    #Using CFL condition to update the time stepper.
    try:
        logger.info('Starting main loop')
        while solver.proceed:
            timestep = CFL.compute_timestep()
            solver.step(timestep)
            if (solver.iteration-1) % 10 == 0:
                max_TKE = flow.max('TKE')
                logger.info('Iteration=%i, Time=%e, dt=%e, max(TKE)=%f' %(solver.iteration, solver.sim_time, timestep, max_TKE))
    except:
        logger.error('Exception raised, triggering end of main loop.')
        raise
    finally:
        solver.log_stats()
