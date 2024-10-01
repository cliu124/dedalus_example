"""
Dedalus script computing the eigenmodes of waves on a clamped string.
This script demonstrates solving a 1D eigenvalue problem and produces
plots of the first few eigenmodes and the relative error of the eigenvalues.
It should take just a few seconds to run (serial only).

We use a Legendre basis to solve the EVP:
    s*u + dx(dx(u)) = 0
    u(x=0) = 0
    u(x=Lx) = 0
where s is the eigenvalue.

For the second derivative on a closed interval, we need two tau terms.
Here we choose to use a first-order formulation, putting one tau term
on an auxiliary first-order variable and another in the PDE, and lifting
both to the first derivative basis.

To run and plot:
    $ python3 waves_on_a_string.py
"""

import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)


# Parameters
Lz = 1
Nz = 128
dtype = np.complex128
Rayleigh = 1.5*10**7
Prandtl = 1
Q=1e6
kx=2*np.pi/5*21
ky=0
timestepper = d3.SBDF2
timestep = 2e-3
stop_sim_time=1000

zi=1j

# Bases
zcoord = d3.Coordinate('z')
dist = d3.Distributor(zcoord, dtype=dtype)
zbasis = d3.Chebyshev(zcoord, size=Nz, bounds=(0, Lz))

# Define derivative and lift operators. 
dz = lambda A: d3.Differentiate(A, zcoord)
lift_basis = zbasis.derivative_basis(1)
lift = lambda A: d3.Lift(A, lift_basis, -1)
conj = lambda A: np.conj(A)

# Fields: all of these are harmonics with second order derivative. 
u = dist.Field(name='u', bases=zbasis)
tau_u_1 = dist.Field(name='tau_u_1')
tau_u_2 = dist.Field(name='tau_u_2')
uz = dz(u) + lift(tau_u_1) # First-order reduction
uzz = dz(uz) + lift(tau_u_2)

v = dist.Field(name='v', bases=zbasis)
tau_v_1 = dist.Field(name='tau_v_1')
tau_v_2 = dist.Field(name='tau_v_2')
vz = dz(v) + lift(tau_v_1) # First-order reduction
vzz = dz(vz) + lift(tau_v_2)

w = dist.Field(name='w', bases=zbasis)
tau_w_1 = dist.Field(name='tau_w_1')
tau_w_2 = dist.Field(name='tau_w_2')
wz = dz(w) + lift(tau_w_1) # First-order reduction
wzz = dz(wz) + lift(tau_w_2)

T = dist.Field(name='T', bases=zbasis)
tau_T_1 = dist.Field(name='tau_T_1')
tau_T_2 = dist.Field(name='tau_T_2')
Tz = dz(T) + lift(tau_T_1) # First-order reduction
Tzz = dz(Tz) + lift(tau_T_2)

# Fields: without second order derivative so we do not need tau there
p = dist.Field(name='p', bases=zbasis)

phi = dist.Field(name='phi',bases=zbasis)

tau_phi_1= dist.Field(name='tau_phi_1')
tau_phi_2 = dist.Field(name='tau_phi_2')
phi_z = dz(phi) + lift(tau_phi_1) # First-order reduction
phi_zz = dz(phi_z) + lift(tau_phi_2) # First-order reduction

Jx = dist.Field(name='Jx',bases=zbasis)
Jy = dist.Field(name='Jy',bases=zbasis)
Jz = dist.Field(name='Jz',bases=zbasis)

#horizontal average mode
U0 = dist.Field(name='U0', bases=zbasis)
tau_U0_1 = dist.Field(name='tau_U0_1')
tau_U0_2 = dist.Field(name='tau_U0_2')
U0z = dz(U0) + lift(tau_U0_1) # First-order reduction
U0zz = dz(U0z) + lift(tau_U0_2)

V0 = dist.Field(name='V0', bases=zbasis)
tau_V0_1 = dist.Field(name='tau_V0_1')
tau_V0_2 = dist.Field(name='tau_V0_2')
V0z = dz(V0) + lift(tau_V0_1) # First-order reduction
V0zz = dz(V0z) + lift(tau_V0_2)

T0 = dist.Field(name='T0', bases=zbasis)
tau_T0_1 = dist.Field(name='tau_T0_1')
tau_T0_2 = dist.Field(name='tau_T0_2')
T0z = dz(T0) + lift(tau_T0_1) # First-order reduction
T0zz = dz(T0z) + lift(tau_T0_2)

#pressure and phi gauge variable
tau_p = dist.Field(name='tau_p')
tau_phi = dist.Field(name='tau_phi')

kappa = (Rayleigh * Prandtl)**(-1/2)
nu = (Rayleigh / Prandtl)**(-1/2)
# Problem

problem = d3.IVP([u, v, w, p, Jx, Jy, Jz, phi, T, U0, V0, T0, \
                  tau_u_1, tau_u_2, tau_v_1, tau_v_2, tau_w_1, tau_w_2, tau_T_1, tau_T_2, \
                  tau_U0_1, tau_U0_2, tau_V0_1, tau_V0_2, tau_T0_1, tau_T0_2, \
                      tau_phi_1, tau_phi_2,\
                          tau_p, tau_phi], namespace=locals())

problem.add_equation("dt(u)+zi*kx*p-nu*(uzz-kx*kx*u-ky*ky*u)-Q*nu*Jy=-zi*kx*U0*u-zi*ky*V0*u-w*U0z")
problem.add_equation("dt(v)+zi*ky*p-nu*(vzz-kx*kx*v-ky*ky*v)+Q*nu*Jx=-zi*kx*U0*v-zi*ky*V0*v-w*V0z")
problem.add_equation("dt(w)+dz(p)-nu*(wzz-kx*kx*w-ky*ky*w)-T=-zi*kx*U0*w-zi*ky*V0*w")
problem.add_equation("zi*kx*u+zi*ky*v+wz+tau_p=0")
problem.add_equation("integ(p)=0")

problem.add_equation("Jx+zi*kx*phi-v=0")
problem.add_equation("Jy+zi*ky*phi+u=0")
problem.add_equation("Jz+phi_z=0")
problem.add_equation("phi_zz-kx*kx*phi-ky*ky*phi-zi*kx*v+zi*ky*u+tau_phi=0")
problem.add_equation("integ(phi)=0")

problem.add_equation("dt(T)-w-kappa*(Tzz-kx*kx*T-ky*ky*T)=-zi*kx*U0*T-zi*ky*V0*T-w*T0z")
problem.add_equation("dt(U0)-nu*U0zz=-dz(conj(w)*u+w*conj(u))")
problem.add_equation("dt(V0)-nu*V0zz=-dz(conj(w)*v+w*conj(v))")
problem.add_equation("dt(T0)-kappa*T0zz=-dz(conj(w)*T+w*conj(T))")

#add B.C.
problem.add_equation("w(z=0)=0")# No penetration
problem.add_equation("w(z=Lz)=0")
problem.add_equation("uz(z=0)=0")# Stress free
problem.add_equation("uz(z=Lz)=0")
problem.add_equation("vz(z=0)=0")# Stress free
problem.add_equation("vz(z=Lz)=0")

problem.add_equation("U0z(z=0)=0")
problem.add_equation("U0z(z=Lz)=0")
problem.add_equation("V0z(z=0)=0")
problem.add_equation("V0z(z=Lz)=0")

problem.add_equation("phi_z(z=0)=0")
problem.add_equation("phi_z(z=Lz)=0")

#B.C. for temperature
problem.add_equation("T(z=0)=0")
problem.add_equation("T(z=Lz)=0")
problem.add_equation("T0(z=0)=0")
problem.add_equation("T0(z=Lz)=0")

delta=1
z = dist.local_grid(zbasis)
w['g']=delta*np.sin(np.pi*z)
T['g']=delta*np.sin(np.pi*z)

solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Main loop
u.change_scales(1)
u_list = [np.copy(u['g'])]
w_list = [np.copy(w['g'])]
T_list = [np.copy(T['g'])]
T0_list = [np.copy(T0['g'])]
t_list = [solver.sim_time]
Nu_list = [np.copy(T0z['g'][1])]
while solver.proceed:
    solver.step(timestep)
    if solver.iteration % 100 == 0:
        logger.info('Iteration=%i, Time=%e, dt=%e, Nu-1=%e' %(solver.iteration, solver.sim_time, timestep, T0z['g'][1]))
        
    if solver.iteration % 25 == 0:
        u.change_scales(1)
        u_list.append(np.copy(u['g']))
        w_list.append(np.copy(w['g']))
        T_list.append(np.copy(T['g']))
        T0_list.append(np.copy(T0['g']))
        Nu_list.append(np.copy(T0z['g'][1]))
        t_list.append(solver.sim_time)
        
        
        

# Plot
plt.figure(figsize=(6, 4))
plt.pcolormesh(z.ravel(), np.array(t_list), np.real(np.array(u_list)), cmap='RdBu_r', shading='gouraud', rasterized=True, clim=(-0.8, 0.8))
plt.xlim(0, Lz)
plt.ylim(0, stop_sim_time)
plt.xlabel('z')
plt.ylabel('t')
plt.title('Single-mode of magnetoconvection')
plt.tight_layout()
plt.savefig('u.png', dpi=200)


plt.figure(figsize=(6, 4))
plt.pcolormesh(z.ravel(), np.array(t_list), np.real(np.array(w_list)), cmap='RdBu_r', shading='gouraud', rasterized=True, clim=(-0.8, 0.8))
plt.xlim(0, Lz)
plt.ylim(0, stop_sim_time)
plt.xlabel('z')
plt.ylabel('t')
plt.title('Single-mode of magnetoconvection')
plt.tight_layout()
plt.savefig('w.png', dpi=200)


plt.figure(figsize=(6, 4))
plt.pcolormesh(z.ravel(), np.array(t_list), np.real(np.array(T_list)), cmap='RdBu_r', shading='gouraud', rasterized=True, clim=(-0.8, 0.8))
plt.xlim(0, Lz)
plt.ylim(0, stop_sim_time)
plt.xlabel('z')
plt.ylabel('t')
plt.title('Single-mode of magnetoconvection')
plt.tight_layout()
plt.savefig('T.png', dpi=200)