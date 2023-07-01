

import numpy as np
import dedalus.public as de


# Bases and domain
xbasis = de.Fourier('x', 64, (-1, 1))
zbasis = de.Chebyshev('z', 64, (-1, 1))
domain = de.Domain([xbasis, zbasis], np.float64)
x, z = domain.grids(1)

# Initial field
u = domain.new_field()
u['g'] = np.sin(4*np.pi*z + 2*np.pi*x)

# Flip analytically
v = domain.new_field()
v['g'] = np.sin(4*np.pi*(-z) + 2*np.pi*x)

# Flip numerically
# First copy coefficients
w = domain.new_field()
w['c'] = u['c']
# Make z dimension grid space, but local
w.require_layout(domain.dist.layouts[1])
# Reverse z data
w.data[:,:] = w.data[:,::-1]

# Check result
print(np.allclose(w['g'], v['g']))

