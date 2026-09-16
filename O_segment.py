import numpy as np
np.set_printoptions(precision=3)
from math import sin, cos, tan, pi, atan, atan2, acos, asin
import matplotlib.pyplot as plt
from time import perf_counter

rc = 1.5
rho = 1
thet0 = 260*pi/180
x0 = 2
y0 = 0
xT = 3
yT = -6
xc = 8
yc = -3

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
# MANAGE the plus/minus that gives us thet1_1 and thet1_2
thet1_1 = acos(thet1_num1/thet1_denom + thet1_num2/thet1_denom)
thet1_2 = acos(thet1_num1/thet1_denom - thet1_num2/thet1_denom)

# MANAGE solutions of theta1 = arccos or theta1 = -arccos


# Get the final theta1
thet1 = thet1_1
# thet1 = thet1_2

# Get P1x and P1y
P1x = x0 - l0*rho*sin(thet0) + l0*rho*sin(thet1)
P1y = y0 + l0*rho*cos(thet0) - l0*rho*cos(thet1)

# Get the angle theta2
# MANAGE the plus/minus that gives us thet2_1 and thet2_2
thet2_num_1 = -delta*rc*yTc + (xTc**4 - xTc*rc**2 + (xTc**2)*(yTc**2))**.5
thet2_num_2 = -delta*rc*yTc - (xTc**4 - xTc*rc**2 + (xTc**2)*(yTc**2))**.5
thet2_1 = acos(thet2_num_1/(xTc**2 + yTc**2))
thet2_2 = acos(thet2_num_2/(xTc**2 + yTc**2))

# MANAGE solutions of theta2 = arccos or theta2 = -arccos


# Get the final theta2
thet2 = thet2_1
# thet2 = thet2_2

# Plot the results

