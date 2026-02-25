import numpy as np
np.set_printoptions(precision=3)
from math import sin, cos, tan, pi, atan
import matplotlib.pyplot as plt
plt.close()

r_c = 2
rho = 1
thet0 = 260*pi/180
x0 = -4
y0 = 4
thetf = 105*pi/180
xf = 1
yf = 4
xc = 1
yc = -1

"""
LSLSL case
"""

def get_thet1_LSLSL(r_c, rho, thet0, x0, y0, thetf, xf, yf, xc, yc):
    x_c0 = xc - x0
    y_c0 = yc - y0
    x_fc = xf - xc
    y_fc = yf - yc
    c_c0 = x_c0**2 + y_c0**2 +2*x_c0*sin(thet0) -2*y_c0*cos(thet0) + r_c**2 + 2*rho**2 + 2*r_c*rho
    c_fc = x_fc**2 + y_fc**2 -2*x_fc*rho*sin(thetf) +2*y_fc*cos(thetf) + r_c**2 + 2*rho**2 + 2*r_c*rho
    # Terms of degree 0
    T_0_1 = c_fc*(x_c0 + rho*sin(thet0))**2 -2*(r_c+rho)*(y_fc + rho*cos(thetf))*(x_c0 + rho*sin(thet0))**2
    T_0_2 = c_c0*(x_fc - rho*sin(thetf))**2 +2*(r_c+rho)*(y_c0 - rho*cos(thet0))*(x_fc - rho*sin(thetf))**2
    T_0 = T_0_1 - T_0_2
    # Terms of degree 1
    T_1_1 = (
        4*(r_c+rho)*(x_fc - rho*sin(thetf))*(x_c0 + rho*sin(thet0))**2 + 4*c_fc*(y_c0-rho*cos(thet0))*(x_c0 + rho*sin(thet0))
        -8*(r_c+rho)*(y_fc + rho*cos(thetf))*(y_c0 - rho*cos(thet0))*(x_c0 + rho*sin(thet0))
    )
    T_1_2 = (
        -4*(r_c+rho)*(x_c0 + rho*sin(thet0))*(x_fc - rho*sin(thetf))**2 + 4*c_c0*(y_fc - rho*cos(thetf))*(x_fc - rho*sin(thetf))
        +8*(r_c+rho)*(y_c0 - rho*cos(thet0))*(y_fc - rho*cos(thetf))*(x_fc - rho*sin(thetf))
    )
    T_1 = T_1_1 - T_1_2
    # Terms of degree 2
    T_2_1 = (
        -c_fc*(x_c0 + rho*sin(thet0))**2 + 6*(r_c+rho)*(y_fc + rho*cos(thetf))*(x_c0 + rho*sin(thet0))**2
        + 16*(r_c+rho)*(x_fc - rho*sin(thetf))*(y_c0 - rho*cos(thet0))*(x_c0 + rho*sin(thet0))
        + 8*(r_c+rho)*(y_fc + rho*cos(thetf))*(y_c0 - rho*cos(thet0))*(x_c0 + rho*sin(thet0))
        + 4*c_fc*(y_c0 - rho*cos(thet0))**2 -8*(r_c+rho)*(y_fc + cos(thetf))*(y_c0 - rho*cos(thet0))**2
    )
    T_2_2 = (
        -c_c0*(x_fc - rho*sin(thetf))**2 - 6*(r_c+rho)*(y_c0 - rho*cos(thet0))*(x_fc - rho*sin(thetf))**2
        -16*(r_c+rho)*(x_c0 + rho*sin(thet0))*(y_fc - rho*cos(thetf))*(x_fc - rho*sin(thetf)) + 4*c_c0*(y_fc - rho*cos(thetf))**2
        +8*(r_c+rho)*(y_c0 - rho*cos(thet0))*(y_fc - rho*cos(thetf))**2
    )
    T_2 = T_2_1 - T_2_2
    # Terms of degree 3
    T_3_1 = (
        -8*(r_c+rho)*(x_fc - rho*sin(thetf))*(x_c0 + rho*sin(thet0))**2
        +8*(r_c+rho)*(y_fc + rho*cos(thetf))*(y_c0 - rho*cos(thet0))*(x_c0 + rho*sin(thet0))
        +16*(r_c+ rho)*(x_fc - rho*sin(thetf))*(y_c0 - rho*cos(thet0))**2
    )
    T_3_2 = (
        8*(r_c+rho)*(x_c0 + rho*sin(thet0))*(x_fc - rho*sin(thetf))**2
        -16*(r_c+rho)*(y_c0 - rho*cos(thet0))*(y_fc - rho*cos(thetf))*(x_fc - rho*sin(thetf))
        -16*(r_c+rho)*(x_c0 + rho*sin(thet0))*(y_fc - rho*cos(thetf))**2
    )
    T_3 = T_3_1 - T_3_2
    # Terms of degree 4
    T_4_1 = (
        -c_fc*(x_c0 + rho*sin(thet0))**2 -6*(r_c+rho)*(y_fc + rho*cos(thetf))*(x_c0 + rho*sin(thet0))**2
        -16*(r_c+rho)*(x_fc - rho*sin(thetf))*(y_c0 - rho*cos(thet0))*(x_c0 + rho*sin(thet0)) + 4*c_fc*(y_c0 - rho*cos(thet0))**2
        +8*(r_c+rho)*(y_fc + rho*cos(thetf))*(y_c0 - rho*cos(thet0))**2
    )
    T_4_2 = (
        -c_c0*(x_fc - rho*sin(thetf))**2 + 6*(r_c+rho)*(y_c0 -rho*cos(thet0))*(x_fc - rho*sin(thetf))**2
        +16*(r_c+rho)*(x_c0 + rho*sin(thet0))*(y_fc - rho*cos(thetf))*(x_fc - rho*sin(thetf)) + 4*c_c0*(y_fc - rho*cos(thetf))**2
        -8*(r_c+rho)*(y_c0 - rho*cos(thet0))*(y_fc - rho*cos(thetf))**2
    )
    T_4 = T_4_1 - T_4_2
    # Terms of degree 5
    T_5_1 = (
        4*(r_c+rho)*(x_fc - rho*sin(thetf))*(x_c0 + rho*sin(thet0))**2 - 4*c_fc*(y_c0 - rho*cos(thet0))*(x_c0 + rho*sin(thet0))
        -8*(r_c+rho)*(y_fc + rho*cos(thetf))*(y_c0 - rho*cos(thet0))*(x_c0 + rho*sin(thet0))
    )
    T_5_2 = (
        -4*(r_c+rho)*(x_c0 + rho*sin(thet0))*(x_fc - rho*sin(thetf))**2 - 4*c_c0*(y_fc - rho*cos(thetf))*(x_fc - rho*sin(thetf))
        +8*(r_c+rho)*(y_c0 - rho*cos(thet0))*(y_fc - rho*cos(thetf))*(x_fc - rho*sin(thetf))
    )
    T_5 = T_5_1 - T_5_2
    # Terms of degree 6
    T_6_1 = (
        c_fc*(x_c0 + rho*sin(thet0))**2 + 2*(r_c+rho)*(y_fc + rho*cos(thetf))*(x_c0 + rho*sin(thet0))**2
    )
    T_6_2 = (
        c_c0*(x_fc - rho*sin(thetf))**2 - 2*(r_c+rho)*(y_c0 - rho*cos(thet0))*(x_fc - rho*sin(thetf))**2
    )
    T_6 = T_6_1 - T_6_2
    # Build the companion matrix
    coeffs = np.array([[T_0/T_6], [T_1/T_6], [T_2/T_6], [T_3/T_6], [T_4/T_6], [T_5/T_6]])
    n = 6 # order of the equation that we need to find the roots of
    sub_mat = np.eye(n-1)
    upp_mat = np.zeros((1,n-1))
    left_mat = np.vstack((upp_mat, sub_mat))
    companion = np.hstack((left_mat, coeffs))
    # Get the roots by finding the eigenvalues of the comp matrix
    # print(companion)
    roots, _ = np.linalg.eig(companion)
    # print(roots)
    real_roots = [np.real(root) for root in roots if np.imag(root) == 0]
    # print(real_roots)
    return real_roots

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

roots = get_thet1_LSLSL(r_c, rho, thet0, x0, y0, thetf, xf, yf, xc, yc)
# print(roots)
thet_1_candidates = [2*atan(half_angle)*180/pi for half_angle in roots]
# print(angles)
c_x0, c_y0 = get_circle_center(x0, y0, thet0, rho, "L")
c_xf, c_yf = get_circle_center(xf, yf, thetf, rho, "L")
c1_candidates = []
for thet1 in thet_1_candidates:
    c_x1, c_y1 = get_circle_center(xc, yc, thet1, rho+r_c, "L")
    c1_candidates.append((c_x1,c_y1))
    # Account for the fact that adding 180 degrees can also give the same tangent
    c_x1, c_y1 = get_circle_center(xc, yc, thet1+pi, rho+r_c, "L")
    c1_candidates.append((c_x1,c_y1))
print(c1_candidates)

fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
arrow_length = 1
# Initial point and circle
plt.plot(x0, y0, "ro")
plt.arrow(x0, y0, arrow_length*cos(thet0), arrow_length*sin(thet0), head_width=0.2)
circle0 = plt.Circle((c_x0,c_y0), rho, fill=False, linestyle="--")
ax.add_patch(circle0)
# Mid circle
circle1 = plt.Circle((xc,yc), r_c, fill=False, linestyle="--", color="b")
ax.add_patch(circle1)
# End point and circle
plt.plot(xf, yf, "ro")
plt.arrow(xf, yf, arrow_length*cos(thetf), arrow_length*sin(thetf), head_width=0.2)
circlef = plt.Circle((c_xf,c_yf), rho, fill=False, linestyle="--")
ax.add_patch(circlef)
# Candidate circles tangent to the target
for c1_xy in c1_candidates:
    circle1 = plt.Circle(c1_xy, rho, fill=False, linestyle="--", color="g")
    ax.add_patch(circle1)


plt.grid(True)
plt.axis("equal")
plt.show()
