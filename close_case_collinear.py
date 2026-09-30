import numpy as np
np.set_printoptions(precision=3)
from math import sin, cos, tan, pi, atan, atan2, acos, asin, sqrt
import matplotlib.pyplot as plt
from time import perf_counter

"""
Will need to add a way to make sure the trajectory does not cross C.
"""

# Known values
rho = 0.75
# rho = 3
# thet0 = 260*pi/180 # given between [0,2pi]
rT = 1
xT = 0
yT = 1
rc = 2
xc = 1
yc = -2

# Find known distances
d_TC = ( (xc - xT)**2 + (yc - yT)**2 )**.5
d_mc = rho + rc

### Verify triangle inequality, collinear case
assert d_TC <= 2*rho + rc, "Condition 5 not met for this type of solution."
assert d_TC >= 2*rho - rT + rc, "Condition 6 not met for this type of solution."

# Angle of d_TC
thet_TC = atan2(yT-yc, xT-xc)

# Find the middle of circle c1
c1x = xc + d_mc*cos(thet_TC)
c1y = yc + d_mc*sin(thet_TC)

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
plt.plot(c1x, c1y, "go")
# assert False, "Intentional code break!"
c1 = plt.Circle((c1x, c1y), rho, fill=False, linestyle="--", color="g")
ax.add_patch(c1)

plt.grid(False)
plt.axis("equal")
plt.show()
