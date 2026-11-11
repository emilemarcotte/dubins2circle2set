import numpy as np
np.set_printoptions(precision=3)
from math import sin, cos, tan, pi, atan, atan2, acos, asin, sqrt, hypot
import matplotlib.pyplot as plt
from time import perf_counter

"""
Will need to add a way to make sure the trajectory does not cross C.
"""

# Known values
# rho = 1
rho = 1.9
# thet0 = 260*pi/180 # given between [0,2pi]
rT = 1
xT = 0
yT = 1
rc = 2
xc = 1
yc = -2

# Initial conditions
x0 = -6
y0 = 1
thet0 = 3*pi/2

# Get center of "turn radius circle" for the first L, NEEDS TO BE GENERALIZED
c0x = x0 - rho*sin(thet0)
c0y = y0 + rho*cos(thet0)

def find_internal_tangent_point(x1, y1, r1, x2, y2, r2):
    """
    Finds the single intersection point of two internally tangent circles 
    (one circle is inside the other and they touch at one point).

    Parameters:
        x1, y1, r1: Center (x, y) and radius of the first circle
        x2, y2, r2: Center (x, y) and radius of the second circle
    Returns:
        A tuple (x, y) representing the intersection point, 
        or None if they are not internally tangent.
    """
    dx = x2 - x1
    dy = y2 - y1
    d = hypot(dx, dy)

    # Check if internally tangent and not concentric
    if d != abs(r1 - r2) or d == 0:
        return None

    ux = dx / d
    uy = dy / d

    if r1 > r2:
        x = x1 + r1 * ux
        y = y1 + r1 * uy
    else:
        x = x1 - r1 * ux
        y = y1 - r1 * uy

    return (x, y)

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

# NEED TO REMOVE THIS LATER
c1x = c1x_2
c1y = c1y_2

# Angle of S segment joining C0 and C1, NEEDS TO BE GENERALIZED
psi = (atan2(c1y-c0y, c1x-c0x))%(2*pi)

### Plot the results
fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
arrow_length = 0.5
ax.set_axis_off() # remove the axes and their numbers

# Initial configuration
plt.plot(x0,y0, "ro")
plt.arrow(x0, y0, arrow_length*cos(thet0), arrow_length*sin(thet0), head_width=0.1)
circleC0 = plt.Circle((c0x,c0y), rho, fill=False, linestyle="--", color="g")
ax.add_patch(circleC0)

# point that starts the first S
xt0 = c0x + rho*cos(psi-pi/2)
yt0 = c0y + rho*sin(psi-pi/2)
# plt.plot(xt0,yt0,"bo")

# Intermediate circle C
plt.plot(xc, yc, "ko")
circleC = plt.Circle((xc,yc), rc, fill=False, linestyle="-", color="b")
ax.add_patch(circleC)
ax.text(3, -1, r'$\mathcal{C}$', fontsize=15)

# Terminal zone T
plt.plot(xT, yT, "ko")
circleT = plt.Circle((xT,yT), rT, fill=False, linestyle="-", color="k")
ax.add_patch(circleT)
ax.text(0,2.1, r'$\mathcal{T}$', fontsize=15)

# point P
# plt.plot(px, py, "ro")

# circle c1, 2 possibilities (NOT 4! c_x1 and c_y1 are NOT independent)
# plt.plot(c1x_1, c1y_1, "mo")
# c11 = plt.Circle((c1x_1, c1y_1), rho, fill=False, linestyle="--", color="m")
# ax.add_patch(c11)
plt.plot(c1x_2, c1y_2, "go")
c22 = plt.Circle((c1x_2, c1y_2), rho, fill=False, linestyle="--", color="g")
ax.add_patch(c22)

# Get the points for the first curved segment, and make sure angle range is ok, NEEDS TO BE GENERALIZED
# if psi < thet0: # changed for L
    # angles0 = np.linspace(thet0-pi, psi+pi, 20) # changed for L
#     print("angles01")
# else:
angles0 = np.linspace(thet0-pi, psi+pi, 20) # changed for L
# print(angles0)
x_L0 = []
y_L0 = []
for angle in angles0:
    x, y = c0x - rho*sin(angle), c0y + rho*cos(angle) # changed for L
    x_L0.append(x)
    y_L0.append(y)
plt.plot(x_L0, y_L0, "g-")
# print(x_L0, y_L0)

# point that ends the first S, and starts the 2nd L, NEEDS TO BE GENERALIZED
xt1 = c1x + rho*cos(psi-pi/2)
yt1 = c1y + rho*sin(psi-pi/2)
# plt.plot(xt1, yt1, "bo")

# Straight line, NEEDS TO BE GENERALIZED
pointsx = np.linspace(xt0, xt1, 20)
pointsy = np.linspace(yt0, yt1, 20)
plt.plot(pointsx, pointsy, "g-")

# find the points where circle c1 intersects with T, NEEDS TO BE GENERALIZED
intersection_points = find_internal_tangent_point(c1x, c1y, rho, xT, yT, rT)
print(intersection_points)
x1,y1 = intersection_points
# plt.plot(x1,y1,"mo")
plt.plot(x1,y1,"ro")
thetf = atan2(y1-c1y,x1-c1x)%(2*pi)
# print(thetf*180/pi)

# Draw the last L, NEEDS TO BE GENERALIZED
angles1 = np.linspace(psi+pi, thetf+pi/2-pi, 20)
print(angles1)
x_L1 = []
y_L1 = []
for angle in angles1:
    x, y = c1x - rho*sin(angle), c1y + rho*cos(angle) # changed for L
    x_L1.append(x)
    y_L1.append(y)
plt.plot(x_L1, y_L1, "g-")
# print(x_L1, y_L1)

plt.grid(False)
plt.axis("equal")
plt.tight_layout() # to remove white margins
plt.show()
