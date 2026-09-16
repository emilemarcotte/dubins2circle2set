import numpy as np
np.set_printoptions(precision=3)
from math import sin, cos, tan, pi, atan, atan2, acos, asin
import matplotlib.pyplot as plt
from time import perf_counter

rc = 1.5
rho = 1
thet0 = 290*pi/180
x0 = 2
y0 = 0
xT = 12
yT = 0
xc = 8
yc = 0

def get_circle_center(x, y, thet, rho, turn):
    """
    Get the coordinates of the circle center, getting those using the coordinates
    of a point on its perimeter, and the tangent heading angle.
    x,y: coordinates of the reference circle
    thet: heading angle
    rho: circle radius
    turn: "L" or "R"
    """
    if turn == "L":
        c_x = x - rho*sin(thet)
        c_y = y + rho*cos(thet)
    elif turn == "R":
        c_x = x + rho*sin(thet)
        c_y = y - rho*cos(thet)
    else:
        raise ValueError("""The "turn" argument is invalid. It should be either "L" or "R".""")
    return c_x, c_y

# Subvariables to simplify
x0c = x0 - xc
y0c = y0 - yc
xTc = xT - xc
yTc = yT - yc

# MANAGE THIS
l0 = 1 # left
# l0 = -1 # right

# MANAGE THIS
# delta = 1 # counterclockwise
delta = -1 # clockwise

# Get angle theta_1
mu = x0c**2 + y0c**2 + 2*l0*rho*(y0c*cos(thet0) - x0c*sin(thet0))
thet1_num1 = (y0c + l0*rho*cos(thet0))*(l0*rho - delta*rc)
thet1_num2 = (x0c - l0*rho*sin(thet0))*(mu - rc**2 + 2*delta*rc*l0*rho)**.5
thet1_denom = mu + rho**2

# Get the angle theta1
# MANAGE the plus/minus that gives us thet1_1 and thet1_2, add ValueError for acos
thet1_1 = acos(thet1_num1/thet1_denom + thet1_num2/thet1_denom)
thet1_2 = acos(thet1_num1/thet1_denom - thet1_num2/thet1_denom)

# MANAGE solutions of theta1 = arccos or theta1 = -arccos


# Get the final theta1
# thet1 = thet1_1
thet1 = thet1_2

thet1 = thet1%(2*pi)

# Get P1x and P1y
P1x = x0 - l0*rho*sin(thet0) + l0*rho*sin(thet1)
P1y = y0 + l0*rho*cos(thet0) - l0*rho*cos(thet1)

# Get the angle theta2
# MANAGE the plus/minus that gives us thet2_1 and thet2_2, add ValueError for acos
thet2_num_1 = -delta*rc*yTc + (xTc**4 - xTc*rc**2 + (xTc**2)*(yTc**2))**.5
thet2_num_2 = -delta*rc*yTc - (xTc**4 - xTc*rc**2 + (xTc**2)*(yTc**2))**.5
thet2_1 = acos(thet2_num_1/(xTc**2 + yTc**2))
thet2_2 = acos(thet2_num_2/(xTc**2 + yTc**2))

# MANAGE solutions of theta2 = arccos or theta2 = -arccos


# Get the final theta2
thet2 = thet2_1
# thet2 = thet2_2

thet2 = thet2%(2*pi)


### Plot the results
fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
arrow_length = 1

c_x0, c_y0 = get_circle_center(x0, y0, thet0, rho, "L")

# Initial point and circle
plt.plot(x0, y0, "ro")
plt.arrow(x0, y0, arrow_length*cos(thet0), arrow_length*sin(thet0), head_width=0.2)
circle0 = plt.Circle((c_x0,c_y0), rho, fill=False, linestyle="--")
ax.add_patch(circle0)

# Mid circle
circle_c = plt.Circle((xc,yc), rc, fill=False, linestyle="--", color="b")
ax.add_patch(circle_c)

# End point and circle
plt.plot(xT, yT, "ro")
plt.arrow(xT, yT, arrow_length*cos(thet2), arrow_length*sin(thet2), head_width=0.2)

# Plot P1
plt.plot(P1x, P1y, "ro")

# Plot first S
S1_x = np.linspace(P1x, P1x+4, 20)
b = P1y - tan(thet1)*P1x
plt.plot(S1_x, tan(thet1)*S1_x + b, "g-")

# Plot first curve
angles1 = np.linspace(thet0-pi, thet1+pi, 20) # changed for L
x_L1 = []
y_L1 = []
for angle in angles1:
    x, y = get_circle_center(c_x0, c_y0, angle, rho, "L") # changed for L
    x_L1.append(x)
    y_L1.append(y)
plt.plot(x_L1, y_L1, "g-")

plt.grid(False)
plt.axis("equal")
plt.show()
