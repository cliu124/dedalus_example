"""
Dedalus script simulating forced 2D periodic incompressible turbulence for visualization. 
It can be ran serially or in parallel, and uses the built-in analysis framework to save data snapshots to HDF5 files.
The `plot_snapshots.py` script can be used to produce plots from the saved data.
The simulation should take about 10 cpu-minutes to run.

The initial flow is in the x-direction and depends only on z. The problem is
non-dimensionalized usign the shear-layer spacing and velocity jump, so the
resulting viscosity and tracer diffusivity are related to the Reynolds and
Schmidt numbers as:

    nu  = kinematic viscosity
    eps = energy injection rate

To run and plot using e.g. 4 processes:
    $ mpiexec -n 4 python3 2D_stochastic_d3_COMPLEX_FOURIER.py
    $ mpiexec -n 4 python3 plot_snapshots.py snapshots/*.h5
"""

import numpy as np
import dedalus.core.operators as ops
import dedalus.public as d3
from mpi4py import MPI
import logging
import dedalus.tools.logging as mpi_logging
import dedalus.extras.flow_tools as flow_tools
import dedalus.tools.parallel as parallel
import matplotlib.pyplot as plt
from dedalus.extras.plot_tools import plot_bot_2d
figkw = {'figsize':(6,4), 'dpi':200}

logger = logging.getLogger(__name__)


# Parameters
Lx, Ly        = 2*np.pi,2*np.pi # Size of domain in x and y directions
Nx, Ny        = 64,64           # Resolution in x and y directions
nu            = 0.005           # Viscosity
u0            = 0.1             # Amplitude of initial condition
dealias       = 3/2             # Dealiasing parameter (real-space grid is larger by factor dealias compared to Fourier grid
stop_sim_time = 100             # Integration time after which simulation terminates
timestepper   = d3.RK443        # Time integration scheme, usually: d3.RK443 = 4th-order Runge-Kutta
max_timestep  = 1e-2            # Largest permitted timestep
dtype         = np.complex128   # Complex data type required for complex Fourier transform
forcing_type  = 0               # 0: Random forcing, 1: Kolmogorov forcing
eps           = 1.0             # Forcing amplitude (injection rate by random forcing / prefactor of Kolmogorov forcing)
 
# Bases
coords = d3.CartesianCoordinates('x', 'y')
dist   = d3.Distributor(coords, dtype=dtype)
xbasis = d3.ComplexFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
ybasis = d3.ComplexFourier(coords['y'], size=Ny, bounds=(0, Ly), dealias=dealias)

# Fields
p = dist.Field(name='p', bases=(xbasis,ybasis))
psi = dist.Field(name='psi', bases=(xbasis,ybasis))
u = dist.VectorField(coords, name='u', bases=(xbasis,ybasis))
tau_p = dist.Field(name='tau_p')
tau_psi = dist.Field(name='tau_psi')
f = dist.VectorField(coords, name='f', bases=(xbasis,ybasis))

# Substitutions
x = dist.local_grid(xbasis,scale=dealias)
y = dist.local_grid(ybasis,scale=dealias)
ex, ey = coords.unit_vector_fields(dist)

# Problem
problem = d3.IVP([u, p, psi, tau_p, tau_psi], namespace=locals())
problem.add_equation("dt(u) + grad(p) - nu*lap(u) = - u@grad(u) + f")
problem.add_equation("div(u) + tau_p = 0")
problem.add_equation("integ(p) = 0") # Pressure gauge
problem.add_equation("div(grad(psi)) + d3.div(d3.skew(u))+tau_psi = 0") #definition of streamfunction as inverse laplacian of vorticity
problem.add_equation("integ(psi) = 0") # Pressure gauge

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

#GlobalArrayReducer
reducer   = flow_tools.GlobalArrayReducer(comm=MPI.COMM_WORLD)

# Initial conditions
x = dist.local_grid(xbasis,scale=len(f['g'][0])/Nx)
y = dist.local_grid(ybasis,scale=len(f['g'][1])/Ny)
kini      = 10*np.pi/Lx             
u['g'][0] =  np.cos(kini * y + 0 * x)
u['g'][1] = np.sin(kini * x + 0 * y)
norm      =  reducer.global_mean(u['g']**2)
u['g']    *= u0/norm

# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=0.1, max_writes=100)
snapshots.add_task(p, name='pressure')
snapshots.add_task(-d3.div(d3.skew(u)), name='vorticity')
snapshots.add_task(psi, name='psi')
#snapshots.add_task(d3.ave(0.5*u@u), name='E_kin')

# Scalar Data
analysis1 = solver.evaluator.add_file_handler("scalar_data", sim_dt=0.01)
analysis1.add_task(d3.ave(0.5*(u@u)), name="Ek")
analysis1.add_task((d3.integ((u@ey)**2) - (d3.integ((u@ex)**2)))/(d3.integ(u@u)), name='m')

# CFL
CFL = d3.CFL(solver, initial_dt=max_timestep, cadence=5, safety=0.4, threshold=0.5, max_change=1e3, min_change=0.1, max_dt=max_timestep)
CFL.add_velocity(u)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=10)
flow.add_property(d3.integ(0.5*(u@u)), name='E_kin_tot')
flow.add_property(d3.ave(0.5*(u@u)), name='E_kin_avg')
flow.add_property(d3.ave(0.5*((u@ex)**2)), name='E_kinx')
flow.add_property(d3.ave(0.5*((u@ey)**2)), name='E_kiny')
flow.add_property((d3.ave((u@ey)*(u@ey)) - (d3.ave((u@ex)*(u@ex))))/d3.ave(u@u), name='m')

flow2 = d3.GlobalFlowProperty(solver, cadence=1)
flow2.add_property(d3.ave(f@f), name = 'fvar')

comm = MPI.COMM_WORLD

#Define wavenumbers
kx = dist.local_modes(xbasis); #print('kx='+str(kx)) 
ky = dist.local_modes(ybasis); #print('ky='+str(ky))
#print('kx=',kx,'ky=',ky)

#DEFINE FORCING BAND
k1 = 10
k2 = 12
#KX, KY = np.meshgrid(kx,ky)
#KA2 = KX**2 + KY**2
#force_range = (KA2 >= k1**2)*(KA2 <=k2**2)
#force_range[KA2 > k2**2] = 0.0
###########################################
# Graphically check forcing range is correct
#plt.contourf(KX,KY,force_range); plt.colorbar(); plt.show()
###########################################
x = dist.local_grid(xbasis,scale=dealias)
y = dist.local_grid(ybasis,scale=dealias)

#iteration index
cnt = 0

# Main loop
try:
    logger.info('Starting main loop')
    while solver.proceed:
        cnt += 1
        np.random.seed(123+cnt)
        #print('DBG len(x),len(f[g][0] = '+str(len(x)),str(len(f['g'][0])))
        timestep = CFL.compute_timestep()
       
        x = dist.local_grid(xbasis,scale=len(f['g'][0])/Nx)
        y = dist.local_grid(ybasis,scale=len(f['g'][1])/Ny)
        
        #update random forcing
        f['g'][0] = np.zeros((len(x),len(y)))
        f['g'][1] = np.zeros((len(x),len(y)))
       
        #f['g'][0] =     np.zeros_like(f['g'][0])
        #f['g'][1] =     np.zeros_like(f['g'][1])        
        
        #RANDOM FORCING
        if forcing_type == 0 :
           #update random forcing
           #f['g'][0] =     np.zeros_like(f['g'][0])
           #f['g'][1] =     np.zeros_like(f['g'][1])
        
           for kx_i in kx[(kx>=0)&(kx<k2)]:
              for ky_j in ky[(ky>=0)&(ky<k2)]:
                 if (kx_i**2 + ky_j**2 >= k1**2) and (kx_i**2 + ky_j**2 <= k2**2):
                    X, Y = np.meshgrid(x,y)
                    #print('inside forcing range, kx,ky = '+str(kx_i)+', '+str(ky_j))
                    phase_x    = np.random.rand(1)
                    #print(np.shape(f['g'][0]))
                    #print(np.shape(np.exp(1.0j*kx_i*x + 1.0j*ky_j*y)))
                    f['g'][0]  += np.real( np.exp(2*np.pi*1.0j*phase_x) * np.exp(1.0j*kx_i*x + 1.0j*ky_j*y) ) 
                    phase_y    = np.random.rand(1)
                    f['g'][1]  += np.real( np.exp(2*np.pi*1.0j*phase_y) * np.exp(1.0j*kx_i*x + 1.0j*ky_j*y) )

           comm.Bcast(f['g'],0)
           
           #plot_bot_2d(f[0], figkw=figkw, title="f['g']"); plt.show()
           #PLOT FORCING
           #if mpi_logging.MPI_RANK == 0:  plt.contourf(f['g'][0]); plt.pause(1)
        
           # PRINT MY RANK
           #print('my rank is'+str(mpi_logging.MPI_RANK))

        
           ### RESCALING FORCING ###
           data = f['g']**2
           tmp = reducer.global_mean(data) 
           #print('tmp='+str(tmp))
        
           f['g'][0] = f['g'][0]/np.sqrt(tmp)
           f['g'][1] = f['g'][1]/np.sqrt(tmp)
                
           f['g'][0] = f['g'][0] * np.sqrt(2*eps/timestep)
           f['g'][1] = f['g'][1] * np.sqrt(2*eps/timestep)
        
        # KOLMOGOROV FORCING: f(x,y) = cos(kf * x) * ey + 0 * ex
        if forcing_type == 1:
           X,Y = np.meshgrid(x,y)
           kf_x = 2*np.pi/Lx
           #print(kf_x)
           f['g'][0] = np.zeros_like(f['g'][0])
           f['g'][1] = eps*np.cos(kf_x * x + 0 *y)

           #if mpi_logging.MPI_RANK == 1:
           #plt.contourf(X,Y,np.transpose(f['g'][1])); plt.pause(1)
 		
        solver.step(timestep)

        if (solver.iteration-1) % 10 == 0:
            # max_w = np.sqrt(flow.max('w2'))
            E_kin_avg = flow.max('E_kin_avg'); 
            E_kin_avg_rate = flow.max('E_kin_avg')/solver.sim_time;
            E_kinx= flow.max('E_kinx');
            #E_kiny= flow.max('E_kiny');
            m     = flow.max('m')
            #fvar  = flow.max('fvar')*timestep/(2*eps)
            #favg  = flow.max('favg')
            max_v_frac = 0.0
            if forcing_type == 1:  
               max_v_frac = reducer.global_max(abs(u['g'][1]))/(eps/(kf_x**2*nu))  # CHECK THAT velocity amplitude in laminar flow state is correct
               logger.info('Iteration=%i, Time=%e, dt=%e, Ekin_avg=%f,Ekin_avg/t_sim=%f, m=%f, max_v_frac=%f' %(solver.iteration, solver.sim_time, timestep, E_kin_avg, E_kin_avg_rate,  m, max_v_frac))
            elif forcing_type == 0: 
               logger.info('Iteration=%i, Time=%e, dt=%e, Ekin_avg=%f,Ekin_avg/t_sim=%f, m=%f' %(solver.iteration, solver.sim_time, timestep, E_kin_avg, E_kin_avg_rate,  m))                   
            
except:     
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    1
    #solver.log_stats()
