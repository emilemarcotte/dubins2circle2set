import numpy as np
np.set_printoptions(precision=3)
from math import sin, cos, tan, pi, atan, atan2, acos, asin, sqrt
import matplotlib.pyplot as plt
from time import perf_counter

"""
Will need to add a way to make sure the trajectory does not cross C.
"""

# Known values
# rho = 1
rho = 2
# thet0 = 260*pi/180 # given between [0,2pi]
rT = 1
xT = 0
yT = 1
rc = 2
xc = 1
yc = -2

# Find known distances
d_TC = ( (xc - xT)**2 + (yc - yT)**2 )**.5
d_mT = abs(rT-rho) # tangent case
d_mc = rho + rc

### Verify triangle inequality, tangent case
assert d_TC <= abs(rT - rho) + rho + rc, "Condition 3 not met for this type of solution."
assert d_TC >= abs( abs(rT-rho) - rho - rc ), "Condition 4 not met for this type of solution."

# Find d_rp
d_rp = (d_mT**2 + d_TC**2 - d_mc**2)/(2*d_TC)

# Find d_mp
d_mp = (d_mT**2 - d_rp**2)**.5

# Find point P, version 1
# thet_TC = atan2(yc-yT, xc-xT)
# px = xT + d_rp*cos(thet_TC)
# py = yT + d_rp*sin(thet_TC)

# Find point P, version 2
v_parx = (xc - xT)
v_pary = (yc - yT)
px = xT + (d_rp/d_TC)*v_parx
py = yT + (d_rp/d_TC)*v_pary

# Find v_parx
v_parx_u = (xc - xT)*(1/d_TC)
v_pary_u = (yc - yT)*(1/d_TC)

# Find c1x, 2 possibilities!!!!!!!!!!!!!!!
c1x_1 = px + d_mp*(-v_pary_u)
c1x_2 = px - d_mp*(-v_pary_u)
# Find c1y, 2 possibilities!!!!!!!!!!!!!!!
c1y_1 = py + d_mp*v_parx_u
c1y_2 = py - d_mp*v_parx_u

### Plot the results
fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
arrow_length = 1

# Intermediate circle C
plt.plot(xc, yc, "ko")
circleC = plt.Circle((xc,yc), rc, fill=False, linestyle="--", color="k")
ax.add_patch(circleC)

# Terminal zone T
plt.plot(xT, yT, "ko")
circleT = plt.Circle((xT,yT), rT, fill=False, linestyle="--", color="k")
ax.add_patch(circleT)

# circle c1, 2 possibilities (NOT 4! c_x1 and c_y1 are NOT independent)
plt.plot(c1x_1, c1y_1, "go")
# assert False, "Intentional code break!"
c11 = plt.Circle((c1x_1, c1y_1), rho, fill=False, linestyle="--", color="g")
ax.add_patch(c11)
plt.plot(c1x_2, c1y_2, "mo")
c22 = plt.Circle((c1x_2, c1y_2), rho, fill=False, linestyle="--", color="m")
ax.add_patch(c22)

# point P
plt.plot(px, py, "ro")

plt.grid(False)
plt.axis("equal")
plt.show()
