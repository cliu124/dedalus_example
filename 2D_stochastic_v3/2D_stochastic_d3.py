"""
Dedalus script simulating a 2D periodic incompressible shear flow with a passive
tracer field for visualization. This script demonstrates solving a 2D periodic
initial value problem. It can be ran serially or in parallel, and uses the
built-in analysis framework to save data snapshots to HDF5 files. The
`plot_snapshots.py` script can be used to produce plots from the saved data.
The simulation should take about 10 cpu-minutes to run.

The initial flow is in the x-direction and depends only on z. The problem is
non-dimensionalized usign the shear-layer spacing and velocity jump, so the
resulting viscosity and tracer diffusivity are related to the Reynolds and
Schmidt numbers as:

    nu = 1 / Reynolds
    D = nu / Schmidt

To run and plot using e.g. 4 processes:
    $ mpiexec -n 4 python3 shear_flow.py
    $ mpiexec -n 4 python3 plot_snapshots.py snapshots/*.h5
"""

import numpy as np
import dedalus.public as d3
import logging
import matplotlib.pyplot as plt
logger = logging.getLogger(__name__)


# Parameters
Lx, Ly = 2*np.pi,2*np.pi
Nx, Ny = 64, 64
nu = 0.0 #0.01
eps = 1.0
u0  = 0.01
dealias = 3/2
stop_sim_time = 20
timestepper = d3.RK443
max_timestep = 1e-2
dtype = np.float64

# Bases
coords = d3.CartesianCoordinates('x', 'y')
dist = d3.Distributor(coords, dtype=dtype)
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
ybasis = d3.RealFourier(coords['y'], size=Ny, bounds=(0, Ly), dealias=dealias)

# Fields
p = dist.Field(name='p', bases=(xbasis,ybasis))
u = dist.VectorField(coords, name='u', bases=(xbasis,ybasis))
tau_p = dist.Field(name='tau_p')
f = dist.VectorField(coords, name='g', bases=(xbasis,ybasis))

# Substitutions
x, y = dist.local_grids(xbasis, ybasis)
ex, ey = coords.unit_vector_fields(dist)

# Problem
problem = d3.IVP([u, p, tau_p], namespace=locals())
problem.add_equation("dt(u) + grad(p) - nu*lap(u) = - u@grad(u) + f")
problem.add_equation("div(u) + tau_p = 0")
problem.add_equation("integ(p) = 0") # Pressure gauge


# Solver
solver = problem.build_solver(timestepper)
#snapshots.add_system(solver.state)
solver.stop_sim_time = stop_sim_time

# Initial conditions
u['g'][0] = u0*(1-2*np.random.rand(len(u['g'][0][:,0]),len(u['g'][0][0,:])))

# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=0.1, max_writes=10)
#snapshots.add_system(solver.state)
snapshots.add_task(p, name='pressure')
snapshots.add_task(-d3.div(d3.skew(u)), name='vorticity')
snapshots.add_task(d3.ave(0.5*u@u), name='E_kin')

# Scalar Data
analysis1 = solver.evaluator.add_file_handler("scalar_data", sim_dt=0.01)
#analysis1.add_task(solver.sim_time,name = "time",dtype = dtype)
analysis1.add_task(d3.ave(0.5*(u@u)), name="Ek")
analysis1.add_task((d3.ave((u@ey)**2) - (d3.ave((u@ex)**2)))/(d3.ave(u@u)), name='m')

# CFL
CFL = d3.CFL(solver, initial_dt=max_timestep, cadence=10, safety=0.2, threshold=0.25, max_change=1.5, min_change=0.5, max_dt=max_timestep)
CFL.add_velocity(u)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=10)
flow.add_property(d3.integ(0.5*(u@u)), name='E_kin_tot')
flow.add_property(d3.ave(0.5*(u@u)), name='E_kin_avg')
flow.add_property(d3.ave(0.5*((u@ex)**2)), name='E_kinx')
flow.add_property(d3.ave(0.5*((u@ey)**2)), name='E_kiny')
flow.add_property((d3.ave((u@ey)**2) - (d3.ave((u@ex)**2)))/d3.ave(u@u), name='m')
flow.add_property(d3.ave(f@f), name = 'fvar')

#Define wavenumbers (order: cos(0*x), -sin(0*x), cos(kmin*x), -sin(kmin*x), cos(2*kmin*x), -sin(2*kmin*x), cos(3*kmin*x), - sin(3*kmin*x), ...; kmin= 2*pi/L)
kx = xbasis.wavenumbers; kx_long = np.array([i//2 for i in range(Nx)])
ky = ybasis.wavenumbers; ky_long = np.array([j//2 for j in range(Ny)])
#print('kx=',kx,'ky=',ky)

k1 = 10
k2 = 12
KX, KY = np.meshgrid(kx_long,ky_long)
KA2 = KX**2 + KY**2
force_range = (KA2 >= k1**2)*(KA2 <=k2**2)
###########################################
# Graphically check forcing range is correct
#plt.contourf(force_range); plt.show()
###########################################
# Main loop
try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        #update random forcing
        f['c'][0], f['c'][1] = 0, 0
        random_numbers1 = np.random.normal(0,2,[len(f['c'][0][:,0]),len(f['c'][0][0,:])])
        random_numbers2 = np.random.normal(0,2,[len(f['c'][1][:,0]),len(f['c'][1][0,:])])
        f['c'][0][force_range] = random_numbers1[force_range]; f['g'][0] /= np.sqrt(2*np.var(f['g'][0]));  
        f['c'][1][force_range] = random_numbers2[force_range]; f['g'][1] /= np.sqrt(2*np.var(f['g'][1]));
        f['g'] *= np.sqrt(2*eps/timestep)
        ########################################

        # check forcing normalization
        #print('FORCE AMP'+str((np.var(f['g'][0])+np.var(f['g'][1]))*timestep/(2*eps)))
        solver.step(timestep)
        if (solver.iteration-1) % 10 == 0:
            # max_w = np.sqrt(flow.max('w2'))
            E_kin_avg = flow.max('E_kin_avg'); 
            E_kin_avg_rate = flow.max('E_kin_avg')/solver.sim_time;
            #E_kinx= flow.max('E_kinx');
            #E_kiny= flow.max('E_kiny');
            m     = flow.max('m')
            fvar  = flow.max('fvar')*timestep/(2*eps)
            logger.info('Iteration=%i, Time=%e, dt=%e, Ekin_avg=%f, Ekin_avg_rate=%f, m=%f, fvar = %f' %(solver.iteration, solver.sim_time, timestep, E_kin_avg, E_kin_avg_rate,  m, fvar))
            
except:     
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()
