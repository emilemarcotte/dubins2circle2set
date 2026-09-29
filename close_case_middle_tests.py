import numpy as np
np.set_printoptions(precision=3)
from math import sin, cos, tan, pi, atan, atan2, acos, asin, sqrt
import matplotlib.pyplot as plt
from time import perf_counter

# Known values
rho = 1.5
# thet0 = 260*pi/180 # given between [0,2pi]
rT = 1
xT = 0
yT = 1
rc = 2
xc = 1
yc = -2

condition = ((xc-xT)**2 + (yc-yT)**2)**.5 - rc
print(condition)
assert condition <= 2*rho, "Condition not met for this type of solution."

# Find known distances
d_TC = ( (xc - xT)**2 + (yc - yT)**2 )**.5
d_mT = rho
d_mc = rho + rc

# Find d_rp
d_rp = (d_mT**2 + d_TC**2 - d_mc**2)/(2*d_TC)

# Find d_mp
d_mp = (d_mT**2 + d_rp**2)**.5

# Find v_parx
v_parx = (xc - xT)*(1/d_TC)
v_pary = (yc - yT)*(1/d_TC)

# Find point p
px = xT + (d_rp/d_TC)*v_parx
py = yT + (d_rp/d_TC)*v_pary

# Find c1x, 2 possibilities!!!!!!!!!!!!!!!
c1x_1 = px + d_mp*(-v_pary)
c1x_2 = px - d_mp*(-v_pary)
# Find c1y, 2 possibilities!!!!!!!!!!!!!!!
c1y_1 = py + d_mp*v_parx
c1y_2 = py - d_mp*v_parx


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
c11 = plt.Circle((c1x_1, c1y_1), rho, fill=False, linestyle="--", color="g")
ax.add_patch(c11)
# plt.plot(c1x_1, c1y_2, "bo")
# c12 = plt.Circle((c1x_1, c1y_2), rho, fill=False, linestyle="--", color="b")
# ax.add_patch(c12)
# plt.plot(c1x_2, c1y_1, "ro")
# c21 = plt.Circle((c1x_2, c1y_1), rho, fill=False, linestyle="--", color="r")
# ax.add_patch(c21)
plt.plot(c1x_2, c1y_2, "mo")
c22 = plt.Circle((c1x_2, c1y_2), rho, fill=False, linestyle="--", color="m")
ax.add_patch(c22)

plt.grid(False)
plt.axis("equal")
plt.show()
