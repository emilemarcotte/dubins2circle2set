import numpy as np
np.set_printoptions(precision=3)
from math import sin, cos, tan, pi, atan, atan2, acos, asin
import matplotlib.pyplot as plt
from time import perf_counter

r_c = 2
rho = 1
thet0 = 260*pi/180
x0 = -4
y0 = 4
# thetf = 105*pi/180
xf = 1
yf = 4
# xf = 5
# yf = -4
xc = 1
yc = -2# + 3 or 1

"""
RSLS case
"""

def get_thet1_RSLS(r_c, rho, theta0, x_0, y_0, x_f, y_f, x_c, y_c):

    T0 = 
    T1 = 
    T2 = 
    T3 = 
    T4 = 
    T5 = 
    T6 = 
    T7 = 
    T8 = 
    T9 = 
    T10 = 
    T11 = 
    T12 = 
    T13 = 
    T14 = 
    T15 = 
    T16 = 
    T17 = 
    T18 = 
    T19 = 
    T20 = 
    T21 = 
    T22 = 
    T23 = 
    T24 = 
    T25 = 
    T26 = 
    T27 = 
    T28 = 
    T29 = 
    T30 = 
    T31 = 
    T32 = 
    T33 = 
    T34 = 
    T35 = 
    T36 = 
    T37 = 
    T38 = 
    T39 = 
    T40 = 
    T41 = 
    T42 = 
    T43 = 
    T44 = 
    T45 = 
    T46 = 
    T47 = 
    T48 = 
    T49 = 
    T50 = 
    T51 = 
    T52 = 
    T53 = 
    T54 = 
    T55 = 
    T56 = 
    T57 = 
    T58 = 
    T59 = 
    T60 = 
    T61 = 
    T62 = 
    T63 = 
    T64 = 
    T65 = 
    T66 = 
    T67 = 
    T68 = 
    T69 = 
    T70 = 
    T71 = 
    T72 = 
    T73 = 
    T74 = 
    T75 = 
    T76 = 
    T77 = 
    T78 = 

    # Build the companion matrix, the highest order coeff is 78
    coeffs = np.array([[T_0/T_78], [T_1/T_78], [T_2/T_78], [T_3/T_78], [T_4/T_78],
                       [T_5/T_78], [T_6/T_78], [T_7/T_78], [T_8/T_78], [T_9/T_78],
                       [T_10/T_78], [T_11/T_78], [T_12/T_78], [T_13/T_78], [T_14/T_78],
                       [T_15/T_78], [T_16/T_78], [T_17/T_78], [T_18/T_78], [T_19/T_78],
                       [T_20/T_78], [T_21/T_78], [T_22/T_78], [T_23/T_78], [T_24/T_78],
                       [T_25/T_78], [T_26/T_78], [T_27/T_78], [T_28/T_78], [T_29/T_78]
                       [T_30/T_78], [T_31/T_78], [T_32/T_78], [T_33/T_78], [T_34/T_78],
                       [T_35/T_78], [T_36/T_78], [T_37/T_78], [T_38/T_78], [T_39/T_78],
                       [T_40/T_78], [T_41/T_78], [T_42/T_78], [T_43/T_78], [T_44/T_78],
                       [T_45/T_78], [T_46/T_78], [T_47/T_78], [T_48/T_78], [T_49/T_78],
                       [T_50/T_78], [T_51/T_78], [T_52/T_78], [T_53/T_78], [T_54/T_78],
                       [T_55/T_78], [T_56/T_78], [T_57/T_78], [T_58/T_78], [T_59/T_78],
                       [T_60/T_78], [T_61/T_78], [T_62/T_78], [T_63/T_78], [T_64/T_78],
                       [T_65/T_78], [T_66/T_78], [T_67/T_78], [T_68/T_78], [T_69/T_78],
                       [T_70/T_78], [T_71/T_78], [T_72/T_78], [T_73/T_78], [T_74/T_78],
                       [T_75/T_78], [T_76/T_78], [T_77/T_78]])
    n = 78 # order of the equation that we need to find the roots of
    sub_mat = np.eye(n-1)
    upp_mat = np.zeros((1,n-1))
    left_mat = np.vstack((upp_mat, sub_mat))
    companion = np.hstack((left_mat, coeffs))
    # Get the roots by finding the eigenvalues of the comp matrix
    roots, _ = np.linalg.eig(companion)
    real_roots = [np.real(root) for root in roots if np.imag(root) == 0]
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

def dist_points(x1, y1, x2, y2):
    """
    Returns the distance between 2 points in a 2D plane.
    """
    return ((x1 - x2)**2 + (y1 - y2)**2)**.5

# Calculate the total execution time
start = perf_counter()

# # Find the candidate values for the angle thet1
# roots = get_thet1_RSRS(r_c, rho, thet0, x0, y0, xf, yf, xc, yc)
# print("Real roots:", roots)
# # end = perf_counter() # To know how quick the code is
# # print("time =", end-start)

# thet1_candidates = [2*atan(half_angle) for half_angle in roots]
# # Account for the fact that adding 180 degrees can also give the same tangent
# thet1_candidates_180 = []
# for angle in thet1_candidates:
#     thet1_candidates_180.append(angle+pi)
# thet1_candidates += thet1_candidates_180
# c1_candidates = []
# for thet1 in thet1_candidates:
#     c_x1_, c_y1_ = get_circle_center(xc, yc, thet1, rho+r_c, "R")
#     c1_candidates.append((c_x1_,c_y1_))
# # Keep only the candidate value that gives the shortest overall path length
# c_x0, c_y0 = get_circle_center(x0, y0, thet0, rho, "R")
# distance = -1
# min_dist_idx = None
# c_x1 = None
# c_y1 = None
# x_t0 = None
# y_t0 = None
# x_t1 = None
# y_t1 = None
# thet1_t = None
# phi = None
# theta_f = None
# for i, c1 in enumerate(c1_candidates):
#     _c_x1, _c_y1 = c1
#     _phi = atan2(c_y0-_c_y1, c_x0-_c_x1) # angle of segment joining the centers of c0 and c1 (same as S segment angle)
#     _x_t0, _y_t0 = get_circle_center(c_x0, c_y0, _phi, rho, "R") # point that starts the first S.
#     _x_t1, _y_t1 = get_circle_center(_c_x1, _c_y1, _phi, rho, "R") # point that ends the first S, and starts the 2nd L.
#     # arc1_length = abs(phi - (thet0 - pi))*rho # arc = theta*rho
#     straight1_length = ( (_x_t0 - _x_t1)**2 + (_y_t0 - _y_t1)**2 )**.5
#     arc1_length = abs(_phi - thet0)*rho # arc = theta*rho
#     # Get the final angle
#     _theta1_t = thet1_candidates[i]
#     L = ( (_c_x1 - xf)**2 + (_c_y1 - yf)**2 )**.5
#     # D = ( L**2 - rho**2 )**.5
#     # print("rho/L", rho/L)
#     _theta = acos(rho/L)
#     # _theta = asin(D/L)
#     # _theta = atan2(D, rho)
#     _delta = acos((xf - xc - (r_c + rho)*sin(_theta1_t))/L) # changed for R
#     _theta_f = _delta + _theta - pi/2 # changed for R
#     arc2_length = abs(_theta_f - _phi)*rho
#     _x_tf, _y_tf = get_circle_center(_c_x1, _c_y1, _theta_f, rho, "R") # point that starts the last S.
#     straight2_length = ( (_x_tf - xf)**2 + (_y_tf - yf)**2 )**.5
#     total_path_length = arc1_length + straight1_length + arc2_length + straight2_length
#     if (i == 0) or (total_path_length < distance):
#     # if (i == 1): # for testing purposes
#         # Only keep the parameters for the shortest path
#         distance = total_path_length
#         min_dist_idx = i
#         c_x1 = _c_x1
#         c_y1 = _c_y1
#         x_t0 = _x_t0
#         y_t0 = _y_t0
#         x_t1 = _x_t1
#         y_t1 = _y_t1
#         thet1_t = _theta1_t
#         phi = _phi
#         theta_f = _theta_f
#     # print("total_path_length:", total_path_length)

end = perf_counter() # To know how quick the code is
print("time =", end-start)

# # Plot the results
# fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
# arrow_length = 1
# # Initial point and circle
# plt.plot(x0, y0, "ro")
# plt.arrow(x0, y0, arrow_length*cos(thet0), arrow_length*sin(thet0), head_width=0.2)
# circle0 = plt.Circle((c_x0,c_y0), rho, fill=False, linestyle="--")
# ax.add_patch(circle0)
# # Mid circle
# circle_c = plt.Circle((xc,yc), r_c, fill=False, linestyle="--", color="b")
# ax.add_patch(circle_c)
# # End point and circle
# plt.plot(xf, yf, "ro")
# plt.arrow(xf, yf, arrow_length*cos(theta_f), arrow_length*sin(theta_f), head_width=0.2)
# # # Candidate circles tangent to the target
# # for c1_xy in c1_candidates:
# #     circle1 = plt.Circle(c1_xy, rho, fill=False, linestyle="--", color="g")
# #     ax.add_patch(circle1)
# circle1 = plt.Circle((c_x1,c_y1), rho, fill=False, linestyle="--", color="k")
# ax.add_patch(circle1)
# # Get and plot the tangent point
# p_xt, p_yt = get_circle_center(xc, yc, thet1_t, r_c, "R")
# plt.plot(p_xt, p_yt, "ro")
# plt.arrow(p_xt, p_yt, arrow_length*cos(thet1_t), arrow_length*sin(thet1_t), head_width=0.2)
# # Plot the point that starts the first S segment
# plt.plot(x_t0, y_t0, "go")
# plt.plot(x_t1, y_t1, "go")
# # Get the points for the first curved segment
# angles0 = np.linspace(thet0+pi, phi, 10) # changed for R
# x_L0 = []
# y_L0 = []
# for angle in angles0:
#     x, y = get_circle_center(c_x0, c_y0, angle, rho, "R")
#     x_L0.append(x)
#     y_L0.append(y)
# plt.plot(x_L0, y_L0, "g-")
# # Get the points for the first straight segment
# x_s01 = np.linspace(x_t0, x_t1, 10)
# y_s01 = np.linspace(y_t0, y_t1, 10)
# plt.plot(x_s01, y_s01, "g-")
# # Get the points for the second curved segment
# angles1 = np.linspace(phi, theta_f-pi, 10) # changed for R
# x_L1 = []
# y_L1 = []
# for angle in angles1:
#     x, y = get_circle_center(c_x1, c_y1, angle, rho, "R")
#     x_L1.append(x)
#     y_L1.append(y)
# plt.plot(x_L1, y_L1, "g-")
# plt.plot(x_L1[-1], y_L1[-1], "go")
# # Get the points for the last straight segment
# x_S_f = np.linspace(x_L1[-1], xf, 10)
# y_S_f = np.linspace(y_L1[-1], yf, 10)
# plt.plot(x_S_f, y_S_f, "g-")

# all_x_trajectory = x_L0 + x_s01.tolist()[1:-1] + x_L1 + x_S_f.tolist()[1:]
# all_y_trajectory = y_L0 + y_s01.tolist()[1:-1] + y_L1 + y_S_f.tolist()[1:]
# plt.plot(all_x_trajectory, all_y_trajectory, "c-")

# plt.grid(False)
# plt.axis("equal")
# plt.show()
