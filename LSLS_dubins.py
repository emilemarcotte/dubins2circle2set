import numpy as np
np.set_printoptions(precision=3)
from math import sin, cos, tan, pi, atan, atan2
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
xc = 1
yc = -1

"""
LSLSL case
"""

def get_thet1_LSLS(r_c, rho, thet0, x0, y0, xf, yf, xc, yc):
    R = r_c + rho
    x_c0 = xc - x0
    y_c0 = yc - y0
    x_fc = xf - xc
    y_fc = yf - yc
    c_c0 = x_c0**2 + y_c0**2 +2*x_c0*sin(thet0) -2*y_c0*cos(thet0) + r_c**2 + 2*rho**2 + 2*r_c*rho

    c_v0 = 
    c_v1 = 
    c_v2 = 
    c_v3 = 
    c_v4 = 
    c_v5 = 
    c_v6 = 
    c_v7 = 
    c_v8 = 
    c_v9 = 
    c_v10 = 
    c_v11 = 
    c_v12 = 
    c_v13 = 
    c_v14 = 
    c_v15 = 
    c_v16 = 
    c_v17 = 
    c_v18 = 
    c_v19 = 
    c_v20 = 
    c_v21 = 
    c_v22 = 
    c_v23 = 
    c_v24 = 
    c_v25 = 
    c_v26 = 
    c_v27 = 
    c_v28 = 

    # Eqs 317 and 318
    c_TT0 = 
    c_TT1 = 
    c_TT2 = 
    c_TT3 = 
    c_TT4 = 
    c_TT5 = 
    c_TT6 = 
    c_TT7 = 
    c_TT8 = 
    c_TT9 = 
    c_TT10 = 
    c_TT11 = 
    c_TT12 = 

    # Eqs 203 and 204
    c_q0 = 
    c_q1 = 
    c_q2 = 
    c_q3 = 
    c_q4 = 
    c_q5 = 
    c_q6 = 
    c_q7 = 
    c_q8 = 

    # Eqs 211 and 212
    c_lambda0 = 
    c_lambda1 = 
    c_lambda2 = 
    c_lambda3 = 
    c_lambda4 = 
    c_lambda5 = 
    c_lambda6 = 
    c_lambda7 = 
    c_lambda8 = 
    c_lambda9 = 
    c_lambda10 = 
    c_lambda11 = 
    c_lambda12 = 

    # Eqs 333 and 334
    c_nd0 = 
    c_nd1 = 
    c_nd2 = 
    c_nd3 = 
    c_nd4 = 
    c_nd5 = 
    c_nd6 = 
    c_nd7 = 
    c_nd8 = 
    c_nd9 = 
    c_nd10 = 
    c_nd11 = 
    c_nd12 = 
    c_nd13 = 
    c_nd14 = 

    # Eqs 222 and 224:
    c_T0 = 
    c_T1 = 
    c_T2 = 
    c_T3 = 
    c_T4 = 
    c_T5 = 
    c_T6 = 
    c_T7 = 
    c_T8 = 
    c_T9 = 
    c_T10 = 
    c_T11 = 
    c_T12 = 

    # Eqs 40 and 41
    c__1 = - 2*x0*R - 2*rho*R*sin(thet0)
    c__2 = 2*y_c0*R -2*rho*R*cos(thet0)

    # Eqs 326 and 327:
    c__3 = c_c0 + c__2
    c__4 = 2*c__1
    c__5 = c_c0 - c__2

    # Eqs 357 and 358:
    c_u0 = c_T0-c_nd0
    c_u1 = c_T1-c_nd1
    c_u2 = c_T2+ c_T0-c_nd2
    c_u3 = c_T3+c_T1-c_nd3
    c_u4 = c_T4+c_T2-c_nd4
    c_u5 = c_T5+c_T3-c_nd5
    c_u6 = c_T6+c_T4-c_nd6
    c_u7 = c_T7+c_T5-c_nd7
    c_u8 = c_T8+c_T6-c_nd8
    c_u9 = c_T9+c_T7-c_nd9
    c_10 = c_T10+c_T8-c_nd10
    c_u11 = c_T11+c_T9-c_nd11
    c_u12 = c_T12+c_T10-c_nd12
    c_u13 = c_T11-c_nd13
    c_u14 = c_T12-c_nd14

    # Eqs 76 and 77
    c_D20 = (x_fc**2 + y_fc**2 + R**2 - 2*y_fc*R - rho**2)
    c_D21 = 4*x_fc*R
    c_D22 = (x_fc**2 + y_fc**2 + R**2 + 2*y_fc*R- rho**2)

    # Eqs 335 and 336:
    c_dd0 = c__3*c_TD0
    c_dd1 = c__3*c_TD1 + c__4*c_TD0
    c_dd2 = c__3*c_TD2 + c__4*c_TD1 + c__5*c_TD0
    c_dd3 = c__3*c_TD3 + c__4*c_TD2 + c__5*c_TD1
    c_dd4 = c__3*c_TD4 + c__4*c_TD3 + c__5*c_TD2
    c_dd5 = c__3*c_TD5 + c__4*c_TD4 + c__5*c_TD3
    c_dd6 = c__3*c_TD6 + c__4*c_TD5 + c__5*c_TD4
    c_dd7 = c__3*c_TD7 + c__4*c_TD6 + c__5*c_TD5
    c_dd8 = c__3*c_TD8 + c__4*c_TD7 + c__5*c_TD6
    c_dd9 = c__3*c_TD9 + c__4*c_TD8 + c__5*c_TD7
    c_dd10 = c__3*c_TD10 + c__4*c_TD9 + c__5*c_TD8
    c_dd11 = c__3*c_TD11 + c__4*c_TD10 + c__5*c_TD9
    c_dd12 = c__3*c_TD12 + c__4*c_TD11 + c__5*c_TD10
    c_dd13 = c__4*c_TD12 + c__5*c_TD11
    c_dd14 = c__5*c_TD12

    # Eqs 362 and 363:
    c_v0 = c_u0*c_u0
    c_v1 = c_u0*c_u1 + c_u1*c_u0
    c_v2 = c_u0*c_u2 + c_u1*c_u1 + c_u2*c_u0
    c_v3 = c_u0*c_u3 + c_u1*c_u2 + c_u2*c_u1 + c_u3*c_u0
    c_v4 = c_u0*c_u4 + c_u1*c_u3 + c_u2*c_u2 + c_u3*c_u1 + c_u4*c_u0
    c_v5 = c_u0*c_u5 + c_u1*c_u4 + c_u2*c_u3 + c_u3*c_u2 + c_u4*c_u1 + c_u5*c_u0
    c_v6 = c_u0*c_u6 + c_u1*c_u5 + c_u2*c_u4 + c_u3*c_u3 + c_u4*c_u2 + c_u5*c_u1 + c_u6*c_u0
    c_v7 = c_u0*c_u7 + c_u1*c_u6 + c_u2*c_u5 + c_u3*c_u4 + c_u4*c_u3 + c_u5*c_u2 + c_u6*c_u1 + c_u7*c_u0
    c_v8 = c_u0*c_u8 + c_u1*c_u7 + c_u2*c_u6 + c_u3*c_u5 + c_u4*c_u4 + c_u5*c_u3 + c_u6*c_u2 + c_u7*c_u1 + c_u8*c_u0
    c_v9 = c_u0*c_u9 + c_u1*c_u8 + c_u2*c_u7 + c_u3*c_u6 + c_u4*c_u5 + c_u5*c_u4 + c_u6*c_u3 + c_u7*c_u2 + c_u8*c_u1 + c_u9*c_u0
    c_v10 = c_u0*c_10 + c_u1*c_u9 + c_u2*c_u8 + c_u3*c_u7 + c_u4*c_u6 + c_u5*c_u5 + c_u6*c_u4 + c_u7*c_u3 + c_u8*c_u2 + c_u9*c_u1 + c_10*c_u0
    c_v11 = c_u0*c_u11 + c_u1*c_10 + c_u2*c_u9 + c_u3*c_u8 + c_u4*c_u7 + c_u5*c_u6 + c_u6*c_u5 + c_u7*c_u4 + c_u8*c_u3 + c_u9*c_u2 + c_10*c_u1 + c_u11*c_u0
    c_v12 = c_u0*c_u12 + c_u1*c_u11 + c_u2*c_10 + c_u3*c_u9 + c_u4*c_u8 + c_u5*c_u7 + c_u6*c_u6 + c_u7*c_u5 + c_u8*c_u4 + c_u9*c_u3 + c_10*c_u2 + c_u11*c_u1 + c_u12*c_u0
    c_v13 = c_u0*c_u13 + c_u1*c_u12 + c_u2*c_u11 + c_u3*c_10 + c_u4*c_u9 + c_u5*c_u8 + c_u6*c_u7 + c_u7*c_u6 + c_u8*c_u5 + c_u9*c_u4 + c_10*c_u3 + c_u11*c_u2 + c_u12*c_u1 + c_u13*c_u0
    c_v14 = c_u0*c_u14 + c_u1*c_u13 + c_u2*c_u12 + c_u3*c_u11 + c_u4*c_10 + c_u5*c_u9 + c_u6*c_u8 + c_u7*c_u7 + c_u8*c_u6 + c_u9*c_u5 + c_10*c_u4 + c_u11*c_u3 + c_u12*c_u2 + c_u13*c_u1 + c_u14*c_u0
    c_v15 = c_u1*c_u14 + c_u2*c_u13 + c_u3*c_u12 + c_u4*c_u11 + c_u5*c_10 + c_u6*c_u9 + c_u7*c_u8 + c_u8*c_u7 + c_u9*c_u6 + c_10*c_u5 + c_u11*c_u4 + c_u12*c_u3 + c_u13*c_u2 + c_u14*c_u1
    c_v16 = c_u2*c_u14 + c_u3*c_u13 + c_u4*c_u12 + c_u5*c_u11 + c_u6*c_10 + c_u7*c_u9 + c_u8*c_u8 + c_u9*c_u7 + c_10*c_u6 + c_u11*c_u5 + c_u12*c_u4 + c_u13*c_u3 + c_u14*c_u2
    c_v17 = c_u3*c_u14 + c_u4*c_u13 + c_u5*c_u12 + c_u6*c_u11 + c_u7*c_10 + c_u8*c_u9 + c_u9*c_u8 + c_10*c_u7 + c_u11*c_u6 + c_u12*c_u5 + c_u13*c_u4 + c_u14*c_u3
    c_v18 = c_u4*c_u14 + c_u5*c_u13 + c_u6*c_u12 + c_u7*c_u11 + c_u8*c_10 + c_u9*c_u9 + c_10*c_u8 + c_u11*c_u7 + c_u12*c_u6 + c_u13*c_u5 + c_u14*c_u4
    c_v19 = c_u5*c_u14 + c_u6*c_u13 + c_u7*c_u12 + c_u8*c_u11 + c_u9*c_10 + c_10*c_u9 + c_u11*c_u8 + c_u12*c_u7 + c_u13*c_u6 + c_u14*c_u5
    c_v20 = c_u6*c_u14 + c_u7*c_u13 + c_u8*c_u12 + c_u9*c_u11 + c_10*c_10 + c_u11*c_u9 + c_u12*c_u8 + c_u13*c_u7 + c_u14*c_u6
    c_v21 = c_u7*c_u14 + c_u8*c_u13 + c_u9*c_u12 + c_10*c_u11 + c_u11*c_10 + c_u12*c_u9 + c_u13*c_u8 + c_u14*c_u7
    c_v22 = c_u8*c_u14 + c_u9*c_u13 + c_10*c_u12 + c_u11*c_u11 + c_u12*c_10 + c_u13*c_u9 + c_u14*c_u8
    c_v23 = c_u9*c_u14 + c_10*c_u13 + c_u11*c_u12 + c_u12*c_u11 + c_u13*c_10 + c_u14*c_u9
    c_v24 = c_10*c_u14 + c_u11*c_u13 + c_u12*c_u12 + c_u13*c_u11 + c_u14*c_10
    c_v25 = c_u11*c_u14 + c_u12*c_u13 + c_u13*c_u12 + c_u14*c_u11
    c_v26 = c_u12*c_u14 + c_u13*c_u13 + c_u14*c_u12
    c_v27 = c_u13*c_u14 + c_u14*c_u13
    c_v28 = c_u14*c_u14

    # Eqs 366 and 267:
    c_w0 = c_dd0*c_dd0
    c_w1 = c_dd0*c_dd1 + c_dd1*c_dd0
    c_w2 = c_dd0*c_dd2 + c_dd1*c_dd1 + c_dd2*c_dd0
    c_w3 = c_dd0*c_dd3 + c_dd1*c_dd2 + c_dd2*c_dd1 + c_dd3*c_dd0
    c_w4 = c_dd0*c_dd4 + c_dd1*c_dd3 + c_dd2*c_dd2 + c_dd3*c_dd1 + c_dd4*c_dd0
    c_w5 = c_dd0*c_dd5 + c_dd1*c_dd4 + c_dd2*c_dd3 + c_dd3*c_dd2 + c_dd4*c_dd1 + c_dd5*c_dd0
    c_w6 = c_dd0*c_dd6 + c_dd1*c_dd5 + c_dd2*c_dd4 + c_dd3*c_dd3 + c_dd4*c_dd2 + c_dd5*c_dd1 + c_dd6*c_dd0
    c_w7 = c_dd0*c_dd7 + c_dd1*c_dd6 + c_dd2*c_dd5 + c_dd3*c_dd4 + c_dd4*c_dd3 + c_dd5*c_dd2 + c_dd6*c_dd1 + c_dd7*c_dd0
    c_w8 = c_dd0*c_dd8 + c_dd1*c_dd7 + c_dd2*c_dd6 + c_dd3*c_dd5 + c_dd4*c_dd4 + c_dd5*c_dd3 + c_dd6*c_dd2 + c_dd7*c_dd1 + c_dd8*c_dd0
    c_w9 = c_dd0*c_dd9 + c_dd1*c_dd8 + c_dd2*c_dd7 + c_dd3*c_dd6 + c_dd4*c_dd5 + c_dd5*c_dd4 + c_dd6*c_dd3 + c_dd7*c_dd2 + c_dd8*c_dd1 + c_dd9*c_dd0
    c_w10 = c_dd0*c_dd10 + c_dd1*c_dd9 + c_dd2*c_dd8 + c_dd3*c_dd7 + c_dd4*c_dd6 + c_dd5*c_dd5 + c_dd6*c_dd4 + c_dd7*c_dd3 + c_dd8*c_dd2 + c_dd9*c_dd1 + c_dd10*c_dd0
    c_w11 = c_dd0*c_dd11 + c_dd1*c_dd10 + c_dd2*c_dd9 + c_dd3*c_dd8 + c_dd4*c_dd7 + c_dd5*c_dd6 + c_dd6*c_dd5 + c_dd7*c_dd4 + c_dd8*c_dd3 + c_dd9*c_dd2 + c_dd10*c_dd1 + c_dd11*c_dd0
    c_w12 = c_dd0*c_dd12 + c_dd1*c_dd11 + c_dd2*c_dd10 + c_dd3*c_dd9 + c_dd4*c_dd8 + c_dd5*c_dd7 + c_dd6*c_dd6 + c_dd7*c_dd5 + c_dd8*c_dd4 + c_dd9*c_dd3 + c_dd10*c_dd2 + c_dd11*c_dd1 + c_dd12*c_dd0
    c_w13 = c_dd0*c_dd13 + c_dd1*c_dd12 + c_dd2*c_dd11 + c_dd3*c_dd10 + c_dd4*c_dd9 +  c_dd5*c_dd8 + c_dd6*c_dd7 + c_dd7*c_dd6 + c_dd8*c_dd5 + c_dd9*c_dd4 + c_dd10*c_dd3 + c_dd11*c_dd2 + c_dd12*c_dd1 + c_dd13*c_dd0
    c_w14 = c_dd0*c_dd14 + c_dd1*c_dd13 + c_dd2*c_dd12 + c_dd3*c_dd11 + c_dd4*c_dd10 + c_dd5*c_dd9 + c_dd6*c_dd8 + c_dd7*c_dd7 + c_dd8*c_dd6 + c_dd9*c_dd5 + c_dd10*c_dd4 + c_dd11*c_dd3 + c_dd12*c_dd2 + c_dd13*c_dd1 + c_dd14*c_dd0
    c_w15 = c_dd1*c_dd14 + c_dd2*c_dd13 + c_dd3*c_dd12 + c_dd4*c_dd11 + c_dd5*c_dd10 + c_dd6*c_dd9 + c_dd7*c_dd8 + c_dd8*c_dd7 + c_dd9*c_dd6 + c_dd10*c_dd5 + c_dd11*c_dd4 + c_dd12*c_dd3 + c_dd13*c_dd2 + c_dd14*c_dd1
    c_w16 = c_dd2*c_dd14 + c_dd3*c_dd13 + c_dd4*c_dd12 + c_dd5*c_dd11 + c_dd6*c_dd10 + c_dd7*c_dd9 + c_dd8*c_dd8 + c_dd9*c_dd7 + c_dd10*c_dd6 + c_dd11*c_dd5 + c_dd12*c_dd4 + c_dd13*c_dd3 + c_dd14*c_dd2
    c_w17 = c_dd3*c_dd14 + c_dd4*c_dd13 + c_dd5*c_dd12 + c_dd6*c_dd11 + c_dd7*c_dd10 + c_dd8*c_dd9 + c_dd9*c_dd8 + c_dd10*c_dd7 + c_dd11*c_dd6 + c_dd12*c_dd5 + c_dd13*c_dd4 + c_dd14*c_dd3
    c_w18 = c_dd4*c_dd14 + c_dd5*c_dd13 + c_dd6*c_dd12 + c_dd7*c_dd11 + c_dd8*c_dd10 + c_dd9*c_dd9 + c_dd10*c_dd8 + c_dd11*c_dd7 + c_dd12*c_dd6 + c_dd13*c_dd5 + c_dd14*c_dd4
    c_w19 = c_dd5*c_dd14 + c_dd6*c_dd13 + c_dd7*c_dd12 + c_dd8*c_dd11 + c_dd9*c_dd10 + c_dd10*c_dd9 + c_dd11*c_dd8 + c_dd12*c_dd7 + c_dd13*c_dd6 + c_dd14*c_dd5
    c_w20 = c_dd6*c_dd14 + c_dd7*c_dd13 + c_dd8*c_dd12 + c_dd9*c_dd11 + c_dd10*c_dd10 + c_dd11*c_dd9 + c_dd12*c_dd8 + c_dd13*c_dd7 + c_dd14*c_dd6
    c_w21 = c_dd7*c_dd14 + c_dd8*c_dd13 + c_dd9*c_dd12 + c_dd10*c_dd11 + c_dd11*c_dd10 + c_dd12*c_dd9 + c_dd13*c_dd8 + c_dd14*c_dd7
    c_w22 = c_dd8*c_dd14 + c_dd9*c_dd13 + c_dd10*c_dd12 + c_dd11*c_dd11 + c_dd12*c_dd10 + c_dd13*c_dd9 + c_dd14*c_dd8
    c_w23 = c_dd9*c_dd14 + c_dd10*c_dd13 + c_dd11*c_dd12 + c_dd12*c_dd11 + c_dd13*c_dd10 + c_dd14*c_dd9
    c_w24 = c_dd10*c_dd14 + c_dd11*c_dd13 + c_dd12*c_dd12 + c_dd13*c_dd11 + c_dd14*c_dd10
    c_w25 = c_dd11*c_dd14 + c_dd12*c_dd13 + c_dd13*c_dd12 + c_dd14*c_dd11
    c_w26 = c_dd12*c_dd14 + c_dd13*c_dd13 + c_dd14*c_dd12
    c_w27 = c_dd13*c_dd14 + c_dd14*c_dd13
    c_w28 = c_dd14*c_dd14

    # Eq 378:
    T_0 = c_v0-c_D20*c_w0
    T_1 = c_v1- c_D20*c_w1 - c_D21*c_w0
    T_2 = c_v2 + c_v0 - c_D20*c_w2 - c_D21*c_w1 - c_D22*c_w0
    T_3 = c_v3 + c_v1 - c_D20*c_w3 - c_D21*c_w2 - c_D22*c_w1
    T_4 = c_v4 + c_v2 - c_D20*c_w4 - c_D21*c_w3 - c_D22*c_w2
    T_5 = c_v5 + c_v3 - c_D20*c_w5 - c_D21*c_w4 - c_D22*c_w3
    T_6 = c_v6 + c_v4 - c_D20*c_w6 - c_D21*c_w5 - c_D22*c_w4
    T_7 = c_v7 + c_v5 - c_D20*c_w7 - c_D21*c_w6 - c_D22*c_w5
    T_8 = c_v8 + c_v6 - c_D20*c_w8 - c_D21*c_w7 - c_D22*c_w6
    T_9 = c_v9 + c_v7 - c_D20*c_w9 - c_D21*c_w8 - c_D22*c_w7
    T_10 = c_v10 + c_v8 - c_D20*c_w10 - c_D21*c_w9 - c_D22*c_w8
    T_11 = c_v11 + c_v9 - c_D20*c_w11 - c_D21*c_w10 - c_D22*c_w9
    T_12 = c_v12 + c_v10 - c_D20*c_w12 - c_D21*c_w11 - c_D22*c_w10
    T_13 = c_v13 + c_v11 - c_D20*c_w13 - c_D21*c_w12 - c_D22*c_w11
    T_14 = c_v14 + c_v12 - c_D20*c_w14 - c_D21*c_w13 - c_D22*c_w12
    T_15 = c_v15 + c_v13 - c_D20*c_w15 - c_D21*c_w14 - c_D22*c_w13
    T_16 = c_v16 + c_v14 - c_D20*c_w16 - c_D21*c_w15 - c_D22*c_w14
    T_17 = c_v17 + c_v15 - c_D20*c_w17 - c_D21*c_w16 - c_D22*c_w15
    T_18 = c_v18 + c_v16 - c_D20*c_w18 - c_D21*c_w17 - c_D22*c_w16
    T_19 = c_v19 + c_v17 - c_D20*c_w19 - c_D21*c_w18 - c_D22*c_w17
    T_20 = c_v20 + c_v18 - c_D20*c_w20 - c_D21*c_w19 - c_D22*c_w18
    T_21 = c_v21 + c_v19 - c_D20*c_w21 - c_D21*c_w20 - c_D22*c_w19
    T_22 = c_v22 + c_v20 - c_D20*c_w22 - c_D21*c_w21 - c_D22*c_w20
    T_23 = c_v23 + c_v21 - c_D20*c_w23 - c_D21*c_w22 - c_D22*c_w21
    T_24 = c_v24 + c_v22 - c_D20*c_w24 - c_D21*c_w23 - c_D22*c_w22
    T_25 = c_v25 + c_v23 - c_D20*c_w25 - c_D21*c_w24 - c_D22*c_w23
    T_26 = c_v26 + c_v24 - c_D20*c_w26 - c_D21*c_w25 - c_D22*c_w24
    T_27 = c_v27 + c_v25 - c_D20*c_w27 - c_D21*c_w26 - c_D22*c_w25
    T_28 = c_v28 + c_v26 - c_D20*c_w28 - c_D21*c_w27 - c_D22*c_w26
    T_29 = c_v27 - c_D21*c_w28 - c_D22*c_w27
    T_30 = c_v28 - c_D22*c_w28
    # Build the companion matrix, the highest order coeff is 1
    coeffs = np.array([[T_0/T_6], [T_1/T_6], [T_2/T_6], [T_3/T_6], [T_4/T_6], [T_5/T_6]])
    n = 6 # order of the equation that we need to find the roots of
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

# Find the candidate values for the angle thet1
start = perf_counter()
roots = get_thet1_LSLSL(r_c, rho, thet0, x0, y0, thetf, xf, yf, xc, yc)
thet1_candidates = [2*atan(half_angle)*180/pi for half_angle in roots]
c_x0, c_y0 = get_circle_center(x0, y0, thet0, rho, "L")
c_xf, c_yf = get_circle_center(xf, yf, thetf, rho, "L")
thet1_candidates_180 = []
for angle in thet1_candidates:
    thet1_candidates_180.append(angle+pi)
# Account for the fact that adding 180 degrees can also give the same tangent
thet1_candidates += thet1_candidates_180
c1_candidates = []
for thet1 in thet1_candidates:
    c_x1, c_y1 = get_circle_center(xc, yc, thet1, rho+r_c, "L")
    c1_candidates.append((c_x1,c_y1))
# Keep only the candidate value that gives the shortest distance between circles c0 and c1
distances = []
for c_x1,c_y1 in c1_candidates:
    dist = dist_points(c_x0, c_y0, c_x1, c_y1)
    distances.append(dist)
min_dist_idx = np.argmin(distances)
c_x1, c_y1 = c1_candidates[min_dist_idx]
end = perf_counter() # To know how quick the code is
print("time =", end-start)

fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
arrow_length = 1
# Initial point and circle
plt.plot(x0, y0, "ro")
plt.arrow(x0, y0, arrow_length*cos(thet0), arrow_length*sin(thet0), head_width=0.2)
circle0 = plt.Circle((c_x0,c_y0), rho, fill=False, linestyle="--")
ax.add_patch(circle0)
# Mid circle
circle_c = plt.Circle((xc,yc), r_c, fill=False, linestyle="--", color="b")
ax.add_patch(circle_c)
# End point and circle
plt.plot(xf, yf, "ro")
plt.arrow(xf, yf, arrow_length*cos(thetf), arrow_length*sin(thetf), head_width=0.2)
circlef = plt.Circle((c_xf,c_yf), rho, fill=False, linestyle="--")
ax.add_patch(circlef)
# Candidate circles tangent to the target
# for c1_xy in c1_candidates:
#     circle1 = plt.Circle(c1_xy, rho, fill=False, linestyle="--", color="g")
#     ax.add_patch(circle1)
circle1 = plt.Circle((c_x1,c_y1), rho, fill=False, linestyle="--", color="k")
ax.add_patch(circle1)
# Get and plot the tangent point
thet1_t = thet1_candidates[min_dist_idx]
p_xt, p_yt = get_circle_center(xc, yc, thet1_t, r_c, "L")
plt.plot(p_xt, p_yt, "ro")
plt.arrow(p_xt, p_yt, arrow_length*cos(thet1_t), arrow_length*sin(thet1_t), head_width=0.2)
# Plot the point that starts the first S segment
# First, angle of segment joining the centers of c0 and c1,
# with an extra pi to get the second tangent point in the circle
phi = pi + atan2(c_y1-c_y0, c_x1-c_x0)
x_t0, y_t0 = get_circle_center(c_x0, c_y0, phi, rho, "L")
plt.plot(x_t0, y_t0, "go")
x_t1, y_t1 = get_circle_center(c_x1, c_y1, phi, rho, "L")
plt.plot(x_t1, y_t1, "go")
# Get the joining segment points
x_s01 = np.linspace(x_t0, x_t1, 10)
slope = tan(phi)
b = y_t0 - slope*x_t0
y_s01 = slope*x_s01 + b
plt.plot(x_s01, y_s01, "g-")
# Do the same thing for the other segment
phi2 = pi + atan2(c_yf-c_y1, c_xf-c_x1)
x_t2, y_t2 = get_circle_center(c_x1, c_y1, phi2, rho, "L")
plt.plot(x_t2, y_t2, "go")
x_tf, y_tf = get_circle_center(c_xf, c_yf, phi2, rho, "L")
plt.plot(x_tf, y_tf, "go")
x_s1f = np.linspace(x_t2, x_tf, 10)
slope2 = tan(phi2)
b2 = y_t2 - slope2*x_t2
y_s1f = slope2*x_s1f + b2
plt.plot(x_s1f, y_s1f, "g-")
# Get the points for the first curved segment
angles0 = np.linspace(thet0-pi, phi, 10)
x_L0 = []
y_L0 = []
for angle in angles0:
    x, y = get_circle_center(c_x0, c_y0, angle, rho, "L")
    x_L0.append(x)
    y_L0.append(y)
plt.plot(x_L0, y_L0, "g-")
# Get the points for the second curved segment
angles1 = np.linspace(phi, phi2, 10)
x_L1 = []
y_L1 = []
for angle in angles1:
    x, y = get_circle_center(c_x1, c_y1, angle, rho, "L")
    x_L1.append(x)
    y_L1.append(y)
plt.plot(x_L1, y_L1, "g-")
# Get the points for the third curved segment
anglesf = np.linspace(phi2, thetf+pi, 10)
x_Lf = []
y_Lf = []
for angle in anglesf:
    x, y = get_circle_center(c_xf, c_yf, angle, rho, "L")
    x_Lf.append(x)
    y_Lf.append(y)
plt.plot(x_Lf, y_Lf, "g-")

all_x_trajectory = x_L0 + x_s01.tolist()[1:-1] + x_L1 + x_s1f.tolist()[1:-1] + x_Lf
all_y_trajectory = y_L0 + y_s01.tolist()[1:-1] + y_L1 + y_s1f.tolist()[1:-1] + y_Lf
# print(all_x_trajectory)
# print(all_y_trajectory)
plt.plot(all_x_trajectory, all_y_trajectory, "c-")

arc1_length = abs(phi - (thet0 - pi))*rho # arc = theta*rho
straight1_length = ( (x_t0 - x_t1)**2 + (y_t0 - y_t1)**2 )**.5
arc2_length = abs(phi2 - phi)*rho
straight2_length = ( (x_t2 - x_tf)**2 + (y_t2 - y_tf)**2 )**.5
arc3_length = abs(thetf+pi - phi2)*rho
total_path_length = arc1_length + straight1_length + arc2_length + straight2_length + arc3_length
print("arc1_length", arc1_length)
print("straight1_length", straight1_length)
print("arc2_length", arc2_length)
print("straight2_length", straight2_length)
print("arc3_length", arc3_length)
print("total_path_length", total_path_length)

plt.grid(False)
plt.axis("equal")
plt.show()
