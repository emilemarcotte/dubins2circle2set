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
# xf = 0
# yf = 0
xc = 1
yc = -2# + 3 or 1

"""
LSLS case
"""

def get_thet1_LSLS(r_c, rho, thet0, x0, y0, xf, yf, xc, yc):
    R = r_c + rho
    x_c0 = xc - x0
    y_c0 = yc - y0
    x_fc = xf - xc
    y_fc = yf - yc

    c_c0 = x_c0**2 + y_c0**2 +2*x_c0*sin(thet0) -2*y_c0*cos(thet0) + r_c**2 + 2*(rho**2) + 2*r_c*rho
    c_sq = x_fc**2 + y_fc**2 + R**2 # Eq 48

    # Eqs 40 and 41:
    c__1 = - 2*x0*R - 2*rho*R*sin(thet0)
    c__2 = 2*y_c0*R -2*rho*R*cos(thet0)

    # Eqs 72 and 73:
    c_L20 = c_sq - 2*y_fc*R
    c_L21 = 4*x_fc*R
    c_L22 = 2*y_fc*R + c_sq

    # Eqs 74 and 75:
    c_L40 = 4*(y_fc**2)*(R**2) - 4*c_sq*y_fc*R + c_sq**2
    c_L41 = - 16*x_fc*y_fc*(R**2) + 8*c_sq*x_fc*R
    c_L42 = 16*(x_fc**2)*(R**2) - 8*(y_fc**2)*(R**2) + 2*(c_sq**2)
    c_L43 = 16*x_fc*y_fc*(R**2) + 8*c_sq*x_fc*R
    c_L44 = 4*(y_fc**2)*(R**2) + 4*c_sq*y_fc*R + c_sq**2

    # Eqs 76 and 77:
    c_D20 = (x_fc**2 + y_fc**2 + R**2 - 2*y_fc*R - rho**2)
    c_D21 = 4*x_fc*R
    c_D22 = (x_fc**2 + y_fc**2 + R**2 + 2*y_fc*R- rho**2)

    # Eq 176:
    c_xm = x_c0*y_c0 - x_c0*rho*cos(thet0) +  rho*y_c0*sin(thet0) - (rho**2)*cos(thet0)*sin(thet0)

    # Between eqs 194 and 195:
    c_d1 = (x_c0+ rho*sin(thet0))**2
    c_d2 = (y_c0- rho*cos(thet0))**2

    # Eqs 203 and 204:
    c_q0 = (c_d1*c_sq -2*c_d1*y_fc*R) + (2*c_xm*c_sq -4*c_xm*y_fc*R )
    c_q1 = 4*c_d1*x_fc*R + 8*c_xm*x_fc*R
    c_q2 = (4*c_d1*y_fc*R) + (4*c_d2*c_sq -8*c_d2*y_fc) + (8*c_xm*c_sq- 8*c_xm*y_fc*R) + 16*y_c0*R*(c_sq-2*y_fc*R)
    c_q3 = -4*c_d1*x_fc*R + 16*c_d2*x_fc*R + 24*c_xm*x_fc*R + 64*x_fc*y_c0*(R**2)
    c_q4 = (-2*c_d1*c_sq ) + (8*c_d2*c_sq ) + (12*c_xm*c_sq) + 64*y_c0*y_fc*(R**2)
    c_q5 = -4*c_d1*x_fc*R + 16*c_d2*x_fc*R + 24*c_xm*x_fc*R -64*x_fc*y_c0*(R**2)
    c_q6 = (-4*c_d1*y_fc*R) + (4*c_d2*c_sq +8*c_d2*y_fc) + (8*c_xm*c_sq + 8*c_xm*y_fc*R) + 16*y_c0*R*(-c_sq - 2*y_fc*R)
    c_q7 = 4*c_d1*x_fc*R + 8*c_xm*x_fc*R
    c_q8 = (c_d1*c_sq +2*c_d1*y_fc*R) + (2*c_xm*c_sq + 4*c_xm*y_fc*R)

    # Eqs 211 and 212
    c_lambda0 = c_L40 + rho**4 -2*(rho**2)*c_L20 + (rho**2)*c_D20
    c_lambda1 = c_L41 -2*(rho**2)*c_L21 + (rho**2)*c_D21
    c_lambda2 = c_L42 + 2*rho**4 -2*(rho**2)*c_L22 -2*(rho**2)*c_L20 + (rho**2)*c_D22 + (rho**2)*c_D20
    c_lambda3 = c_L43 -2*(rho**2)*c_L21 + (rho**2)*c_D21
    c_lambda4 = c_L44 + rho**4 -2*(rho**2)*c_L22 + (rho**2)*c_D22

    # Eqs 222 and 224:
    c_T0 = c_lambda0*c_q0
    c_T1 = c_lambda0*c_q1 + c_lambda1*c_q0
    c_T2 = c_lambda0*c_q2 + c_lambda1*c_q1 + c_lambda2*c_q0
    c_T3 = c_lambda0*c_q3 + c_lambda1*c_q2 + c_lambda2*c_q1 + c_lambda3*c_q0
    c_T4 = c_lambda0*c_q4 + c_lambda1*c_q3 + c_lambda2*c_q2 + c_lambda3*c_q1 + c_lambda4*c_q0
    c_T5 = c_lambda0*c_q5 + c_lambda1*c_q4 + c_lambda2*c_q3 + c_lambda3*c_q2 + c_lambda4*c_q1
    c_T6 = c_lambda0*c_q6 + c_lambda1*c_q5 + c_lambda2*c_q4 + c_lambda3*c_q3 + c_lambda4*c_q2
    c_T7 = c_lambda0*c_q7 + c_lambda1*c_q6 + c_lambda2*c_q5 + c_lambda3*c_q4 + c_lambda4*c_q3
    c_T8 = c_lambda0*c_q8 + c_lambda1*c_q7 + c_lambda2*c_q6 + c_lambda3*c_q5 + c_lambda4*c_q4
    c_T9 = c_lambda1*c_q8 + c_lambda2*c_q7 + c_lambda3*c_q6 + c_lambda4*c_q5
    c_T10 = c_lambda2*c_q8 + c_lambda3*c_q7 + c_lambda4*c_q6
    c_T11 = c_lambda3*c_q8 + c_lambda4*c_q7
    c_T12 = c_lambda4*c_q8

    # Eqs 257 and 259:
    c_a00 = x_fc**2
    c_a01 = 4*x_fc*y_fc
    c_a02 = 2*(2*y_fc**2 -x_fc**2 + x_fc**2)
    c_a03 = 4*x_fc*y_fc
    c_a04 = 2*(-(x_fc**2) + 4*y_fc**2)
    c_a05 = -4*x_fc*y_fc
    c_a06 = 4*y_fc**2
    c_a07 = - 4*x_fc*y_fc
    c_a08 = x_fc**2

    c_a11 = 8*(rho**2)*R*(x_fc)
    c_a12 = 8*(rho**2)*R*(2*y_fc)
    c_a13 = 8*(rho**2)*R*(-x_fc)
    c_a15 = 8*(rho**2)*R*(- x_fc)
    c_a16 = 8*(rho**2)*R*(- 2*y_fc)
    c_a17 = 8*(rho**2)*R*(x_fc)

    c_a22 = 8*(rho**2)*(R**2-R)*(1)
    c_a24 = 8*(rho**2)*(R**2-R)*(-2)
    c_a26 = 8*(rho**2)*(R**2-R)*(1)

    c_a32 = 16*(rho**4)*(R**2)*(1)
    c_a34 = 16*(rho**4)*(R**2)*(-2)
    c_a36 = 16*(rho**4)*(R**2)*(1)

    c_a40 = 2*rho*(x_fc*(y_fc + R))
    c_a42 = 2*rho*(2*x_fc*( - R + 3*y_fc ))
    c_a43 = 2*rho*(8*y_fc*R)
    c_a44 = 2*rho*( 10*x_fc*y_fc)
    c_a45 = 2*rho*(8*y_fc*R)
    c_a46 = 2*rho*(2*x_fc*(R + 2*y_fc))
    c_a48 = 2*rho*(x_fc*(-y_fc -R ))

    c_a51 = 4*rho*( y_fc*R -y_fc**2 + x_fc**2 )
    c_a52 = 4*rho*(2*x_fc*R)
    c_a53 = 4*rho*(-y_fc*R -y_fc**2 + x_fc**2)
    c_a55 = 4*rho*( -y_fc*R +y_fc**2 - x_fc**2 )
    c_a56 = 4*rho*(- 2*x_fc*R)
    c_a57 = 4*rho*(y_fc*R +y_fc**2 - x_fc**2)

    c_a61 = 4*(rho**3)*(-3*y_fc*R + y_fc**2 +2*(R**2) -x_fc**2)
    c_a62 = 4*(rho**3)*(2*x_fc*R)
    c_a63 = 4*(rho**3)*(y_fc**2 -x_fc**2 + 2*(R**2) + 3*y_fc*R)
    c_a65 = 4*(rho**3)*(-y_fc**2 +x_fc**2 + 3*y_fc*R -2*(R**2))
    c_a66 = 4*(rho**3)*(-2*x_fc*R)
    c_a67 = 4*(rho**3)*(-y_fc**2 +x_fc**2 - 3*y_fc*R -2*(R**2))

    c_a70 = (rho**2)*(y_fc**2 -2*y_fc*R+ R**2)
    c_a72 = (rho**2)*(4*(y_fc*R -R**2+ x_fc**2))
    c_a73 = (rho**2)*(16*x_fc*R)
    c_a74 = (rho**2)*(2*(-y_fc**2+ 4*x_fc**2 + 11*R**2))
    c_a75 = (rho**2)*(16*x_fc*R)
    c_a76 = (rho**2)*(4*(-y_fc*R-R**2 + x_fc**2))
    c_a78 = (rho**2)*(y_fc**2 +2*y_fc*R +R**2)

    c_a81 = 4*rho*(x_fc*(-y_fc + R))
    c_a82 = 4*rho*(2*R*(-y_fc + R))
    c_a83 = 4*rho*(x_fc*(-y_fc -R))
    c_a84 = 4*rho*(-4*R**2)
    c_a85 = 4*rho*(x_fc*(y_fc -R))
    c_a86 = 4*rho*(2*R*(y_fc  +R))
    c_a87 = 4*rho*(x_fc*(y_fc + R))

    # Eqs 263 and 264:
    c_Lm0 = c_L20- rho**2
    c_Lm1 = c_L21
    c_Lm2 = c_L22 + c_L20 -2*rho**2
    c_Lm3 = c_L21
    c_Lm4 = c_L22 -(rho**2)

    # Eqs 267 and 268:
    c_Lp0 = c_L40 -2*(rho**2)*c_L20 + rho**4
    c_Lp1 = c_L41 -2*(rho**2)*c_L21
    c_Lp2 = c_L42 -2*(rho**2)*c_L22 -2*(rho**2)*c_L20+ 2*rho**4
    c_Lp3 = c_L43 -2*(rho**2)*c_L21
    c_Lp4 = c_L44 -2*(rho**2)*c_L22 + rho**4

    # Eqs 284 and 285:
    c_b0 = c_Lp0*c_a00
    c_b1 = c_Lp0*c_a01 + c_Lp1*c_a00
    c_b2 = c_Lp0*c_a02 + c_Lp1*c_a01 + c_Lp2*c_a00
    c_b3 = c_Lp0*c_a03 + c_Lp1*c_a02 + c_Lp2*c_a01 + c_Lp3*c_a00
    c_b4 = c_Lp0*c_a04 + c_Lp1*c_a03 + c_Lp2*c_a02 + c_Lp3*c_a01 + c_Lp4*c_a00
    c_b5 = c_Lp0*c_a05 + c_Lp1*c_a04 + c_Lp2*c_a03 + c_Lp3*c_a02 + c_Lp4*c_a01
    c_b6 = c_Lp0*c_a06 + c_Lp1*c_a05 + c_Lp2*c_a04 + c_Lp3*c_a03 + c_Lp4*c_a02
    c_b7 = c_Lp0*c_a07 + c_Lp1*c_a06 + c_Lp2*c_a05 + c_Lp3*c_a04 + c_Lp4*c_a03
    c_b8 = c_Lp0*c_a08 + c_Lp1*c_a07 + c_Lp2*c_a06 + c_Lp3*c_a05 + c_Lp4*c_a04
    c_b9 = c_Lp1*c_a08 + c_Lp2*c_a07 + c_Lp3*c_a06 + c_Lp4*c_a05
    c_b10 = c_Lp2*c_a08 + c_Lp3*c_a07 + c_Lp4*c_a06
    c_b11 = c_Lp3*c_a08 + c_Lp4*c_a07
    c_b12 = c_Lp4*c_a08

    # Eqs 289 and 290:
    c_e1 = c_Lm0*c_a11
    c_e2 = c_Lm0*c_a12 + c_Lm1*c_a11
    c_e3 = c_Lm0*c_a13 +  c_Lm1*c_a12 + c_Lm2*c_a11
    c_e4 = c_Lm1*c_a13 + c_Lm2*c_a12 + c_Lm3*c_a11
    c_e5 = c_Lm0*c_a15 + c_Lm2*c_a13 + c_Lm3*c_a12 + c_Lm4*c_a11
    c_e6 = c_Lm0*c_a16 + c_Lm1*c_a15 + c_Lm3*c_a13 + c_Lm4*c_a12
    c_e7 = c_Lm0*c_a17 + c_Lm1*c_a16 + c_Lm2*c_a15 + c_Lm4*c_a13
    c_e8 = c_Lm1*c_a17 + c_Lm2*c_a16 + c_Lm3*c_a15
    c_e9 = c_Lm2*c_a17 + c_Lm3*c_a16 + c_Lm4*c_a15
    c_e10 = c_Lm3*c_a17 + c_Lm4*c_a16
    c_e11 = c_Lm4*c_a17

    # Eqs 294 and 295:
    c_h2 = c_L20*c_a22
    c_h3 = c_L21*c_a22
    c_h4 = c_L20*c_a24 + (c_L22 +  c_L20)*c_a22
    c_h5 = c_L21*c_a24 + c_L21*c_a22
    c_h6 = c_L20*c_a26 + (c_L22 +  c_L20)*c_a24 + c_L22*c_a22
    c_h7 = c_L21*c_a26 + c_L21*c_a24
    c_h8 = (c_L22 +  c_L20)*c_a26 + c_L22*c_a24
    c_h9 = c_L21*c_a26
    c_h10 = c_L22*c_a26

    # Eqs 298 and 299:
    c_i0 = c_Lm0*c_a40
    c_i1 = c_Lm1*c_a40
    c_i2 = c_Lm0*c_a42 +  c_Lm2*c_a40
    c_i3 = c_Lm0*c_a43 + c_Lm1*c_a42 + c_Lm3*c_a40
    c_i4 = c_Lm0*c_a44 + c_Lm1*c_a43 +  c_Lm2*c_a42 + c_Lm4*c_a40
    c_i5 = c_Lm0*c_a45 + c_Lm1*c_a44 + c_Lm2*c_a43 + c_Lm3*c_a42
    c_i6 = c_Lm0*c_a46 + c_Lm1*c_a45 + c_Lm2*c_a44 + c_Lm3*c_a43 + c_Lm4*c_a42
    c_i7 = c_Lm1*c_a46 + c_Lm2*c_a45 + c_Lm3*c_a44 + c_Lm4*c_a43
    c_i8 = c_Lm0*c_a48 + c_Lm2*c_a46 + c_Lm3*c_a45 + c_Lm4*c_a44
    c_i9 = c_Lm1*c_a48 + c_Lm3*c_a46 + c_Lm4*c_a45
    c_i10 = c_Lm2*c_a48 + c_Lm4*c_a46
    c_i11 = c_Lm3*c_a48
    c_i12 = c_Lm4*c_a48

    # Eqs 303 and 304
    c_j1 = c_L20*c_a51
    c_j2 = c_L20*c_a52 + c_L21*c_a51
    c_j3 = c_L20*c_a53 + c_L21*c_a52 + (c_L22 +  c_L20)*c_a51
    c_j4 = c_L21*c_a53 + (c_L22 + c_L20)*c_a52 + c_L21*c_a51
    c_j5 = c_L20*c_a55 + (c_L22 + c_L20)*c_a53 + c_L21*c_a52 + c_L22*c_a51
    c_j6 = c_L20*c_a56 + c_L21*c_a55 + c_L21*c_a53 + c_L22*c_a52
    c_j7 = c_L20*c_a57 + c_L21*c_a56 + (c_L22 +  c_L20)*c_a55 + c_L22*c_a53
    c_j8 = c_L21*c_a57 + (c_L22 +  c_L20)*c_a56 + c_L21*c_a55
    c_j9 = (c_L22 + c_L20)*c_a57 + c_L21*c_a56 + c_L22*c_a55
    c_j10 = c_L21*c_a57 + c_L22*c_a56
    c_j11 = c_L22*c_a57

    # Eqs 308 and 309:
    c_k0 = c_D20*c_a70
    c_k1 = c_D21*c_a70
    c_k2 = c_D20*c_a72 + (c_D22 + c_D20)*c_a70
    c_k3 = c_D20*c_a73 + c_D21*c_a72 + c_D21*c_a70
    c_k4 = c_D20*c_a74 + c_D21*c_a73 + (c_D22 + c_D20)*c_a72 + c_D22*c_a70
    c_k5 = c_D20*c_a75 + c_D21*c_a74 + (c_D22 + c_D20)*c_a73 + c_D21*c_a72
    c_k6 = c_D20*c_a76 + c_D21*c_a75 + (c_D22 + c_D20)*c_a74 + c_D21*c_a73 + c_D22*c_a72
    c_k7 = c_D21*c_a76 + (c_D22 + c_D20)*c_a75 + c_D21*c_a74 + c_D22*c_a73
    c_k8 = c_D20*c_a78 + (c_D22 + c_D20)*c_a76 + c_D21*c_a75 +  c_D22*c_a74
    c_k9 = c_D21*c_a78 + c_D21*c_a76 + c_D22*c_a75
    c_k10 = (c_D22 + c_D20)*c_a78 + c_D22*c_a76
    c_k11 = c_D21*c_a78
    c_k12 = c_D22*c_a78

    # Eqs 313 and 314:
    c_n1 = c_L20*c_D20*c_a81
    c_n2 = c_L20*c_D20*c_a82 + (c_L20*c_D21 + c_L21*c_D20)*c_a81
    c_n3 = c_L20*c_D20*c_a83 + (c_L20*c_D21 + c_L21*c_D20)*c_a82 + (c_L20*c_D22 + c_L21*c_D21 + c_L22*c_D20)*c_a81
    c_n4 = c_L20*c_D20*c_a84 + (c_L20*c_D21 + c_L21*c_D20)*c_a83 + (c_L20*c_D22 + c_L21*c_D21 + c_L22*c_D20)*c_a82 + (c_L21*c_D22 + c_L22*c_D21)*c_a81
    c_n5 = c_L20*c_D20*c_a85 + (c_L20*c_D21 + c_L21*c_D20)*c_a84 + (c_L20*c_D22 + c_L21*c_D21 + c_L22*c_D20)*c_a83 + (c_L21*c_D22 + c_L22*c_D21)*c_a82 + c_L22*c_D22*c_a81
    c_n6 = c_L20*c_D20*c_a86 + (c_L20*c_D21 + c_L21*c_D20)*c_a85 + (c_L20*c_D22 + c_L21*c_D21 + c_L22*c_D20)*c_a84 + (c_L21*c_D22 + c_L22*c_D21)*c_a83 + c_L22*c_D22*c_a82
    c_n7 = c_L20*c_D20*c_a87 + (c_L20*c_D21 + c_L21*c_D20)*c_a86 + (c_L20*c_D22 + c_L21*c_D21 + c_L22*c_D20)*c_a85 + (c_L21*c_D22 + c_L22*c_D21)*c_a84 + c_L22*c_D22*c_a83
    c_n8 = (c_L20*c_D21 + c_L21*c_D20)*c_a87 + (c_L20*c_D22 + c_L21*c_D21 + c_L22*c_D20)*c_a86 + (c_L21*c_D22 + c_L22*c_D21)*c_a85 + c_L22*c_D22*c_a84
    c_n9 = (c_L20*c_D22 + c_L21*c_D21 + c_L22*c_D20)*c_a87 + (c_L21*c_D22 + c_L22*c_D21)*c_a86 + c_L22*c_D22*c_a85
    c_n10 = (c_L21*c_D22 + c_L22*c_D21)*c_a87 + c_L22*c_D22*c_a86
    c_n11 = c_L22*c_D22*c_a87

    # Eqs 317 and 318:
    c_TT0 = c_b0 + c_k0
    c_TT1 = c_b1 + c_e1 + c_k1 + c_n1
    c_TT2 = c_b2 + c_e2 +  c_h2 + c_a32 + c_k2 + c_n2
    c_TT3 = c_b3 + c_e3 + c_h3 + c_k3 + c_n3
    c_TT4 = c_b4 + c_e4 + c_h4 + (c_a34 + 2*c_a32) + c_k4 + c_n4
    c_TT5 = c_b5 + c_e5 + c_h5 + c_k5 + c_n5
    c_TT6 = c_b6 + c_e6 + c_h6 + (c_a36 + 2*c_a34+ c_a32) + c_k6 +  c_n6
    c_TT7 = c_b7 + c_e7 + c_h7 + c_k7 + c_n7
    c_TT8 = c_b8 + c_e8 + c_h8 + (2*c_a36 + c_a34) + c_k8 + c_n8
    c_TT9 = c_b9 + c_e9 + c_h9 + c_k9 + c_n9
    c_TT10 = c_b10 + c_e10 + c_h10 +  c_a36 + c_k10 + c_n10
    c_TT11 = c_b11 + c_e11 + c_k11 + c_n11
    c_TT12 = c_b12 + c_k12

    # Eqs 319 and 320:
    c_TD0 = c_i0
    c_TD1 = c_i1 + c_j1 + c_a61
    c_TD2 = c_i2 + c_j2 + c_a62
    c_TD3 = c_i3 + c_j3 + (c_a63 + 2*c_a61)
    c_TD4 = c_i4 + c_j4 + 2*c_a62
    c_TD5 = c_i5 + c_j5 + (c_a65 + 2*c_a63 + c_a61)
    c_TD6 = c_i6 + c_j6 + (c_a66 + c_a62)
    c_TD7 = c_i7 + c_j7 + (c_a67 + 2*c_a65 + c_a63)
    c_TD8 = c_i8 + c_j8 + 2*c_a66
    c_TD9 = c_i9 + c_j9 + (2*c_a67 + c_a65)
    c_TD10 = c_i10 + c_j10 + c_a66
    c_TD11 = c_i11 + c_j11 + c_a67
    c_TD12 = c_i12

    # Eqs 326 and 327:
    c__3 = c_c0 + c__2
    c__4 = 2*c__1
    c__5 = c_c0 - c__2


    # Eqs 333 and 334
    c_nd0 = c__3*c_TT0
    c_nd1 = c__3*c_TT1 + c__4*c_TT0
    c_nd2 = c__3*c_TT2 + c__4*c_TT1 + c__5*c_TT0
    c_nd3 = c__3*c_TT3 + c__4*c_TT2 + c__5*c_TT1
    c_nd4 = c__3*c_TT4 + c__4*c_TT3 + c__5*c_TT2
    c_nd5 = c__3*c_TT5 + c__4*c_TT4 + c__5*c_TT3
    c_nd6 = c__3*c_TT6 + c__4*c_TT5 + c__5*c_TT4
    c_nd7 = c__3*c_TT7 + c__4*c_TT6 + c__5*c_TT5
    c_nd8 = c__3*c_TT8 + c__4*c_TT7 + c__5*c_TT6
    c_nd9 = c__3*c_TT9 + c__4*c_TT8 + c__5*c_TT7
    c_nd10 = c__3*c_TT10 + c__4*c_TT9 + c__5*c_TT8
    c_nd11 = c__3*c_TT11 + c__4*c_TT10 + c__5*c_TT9
    c_nd12 = c__3*c_TT12 +c__4*c_TT11 +  c__5*c_TT10
    c_nd13 = c__4*c_TT12 + c__5*c_TT11
    c_nd14 = c__5*c_TT12

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

    # Eqs 366 and 367:
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
    coeffs = np.array([[T_0/T_30], [T_1/T_30], [T_2/T_30], [T_3/T_30], [T_4/T_30],
                       [T_5/T_30], [T_6/T_30], [T_7/T_30], [T_8/T_30], [T_9/T_30],
                       [T_10/T_30], [T_11/T_30], [T_12/T_30], [T_13/T_30], [T_14/T_30],
                       [T_15/T_30], [T_16/T_30], [T_17/T_30], [T_18/T_30], [T_19/T_30],
                       [T_20/T_30], [T_21/T_30], [T_22/T_30], [T_23/T_30], [T_24/T_30],
                       [T_25/T_30], [T_26/T_30], [T_27/T_30], [T_28/T_30], [T_29/T_30]])
    n = 30 # order of the equation that we need to find the roots of
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

# Find the candidate values for the angle thet1
roots = get_thet1_LSLS(r_c, rho, thet0, x0, y0, xf, yf, xc, yc)
print("Real roots:", roots)
thet1_candidates = [2*atan(half_angle) for half_angle in roots]
# Account for the fact that adding 180 degrees can also give the same tangent
thet1_candidates_180 = []
for angle in thet1_candidates:
    thet1_candidates_180.append(angle+pi)
thet1_candidates += thet1_candidates_180
c1_candidates = []
for thet1 in thet1_candidates:
    c_x1_, c_y1_ = get_circle_center(xc, yc, thet1, rho+r_c, "L")
    c1_candidates.append((c_x1_,c_y1_))
# Keep only the candidate value that gives the shortest overall path length
c_x0, c_y0 = get_circle_center(x0, y0, thet0, rho, "L")
distance = np.inf
min_dist_idx = None
c_x1 = None
c_y1 = None
x_t0 = None
y_t0 = None
x_t1 = None
y_t1 = None
thet1_t = None
psi = None
theta_f = None

for i, c1 in enumerate(c1_candidates):
    _c_x1, _c_y1 = c1

    # angle of segment joining the centers of c0 and c1 (same as S segment angle)
    _psi = atan2(_c_y1-c_y0, _c_x1-c_x0)

    # Get the final angle
    _theta1_t = thet1_candidates[i]
    L = ( (_c_x1 - xf)**2 + (_c_y1 - yf)**2 )**.5
    try:
        _theta = acos(rho/L)
        _delta = atan2(yf-_c_y1, xf-_c_x1)
    except (ValueError, ZeroDivisionError):
        # The separation distance L is too short
        continue
    _theta_f = (_delta - _theta + pi/2)%(2*pi)

    # point that starts the first S
    _x_t0 = c_x0 + rho*cos(_psi-pi/2)
    _y_t0 = c_y0 + rho*sin(_psi-pi/2)

    # point that ends the first S, and starts the 2nd L
    _x_t1 = _c_x1 + rho*cos(_psi-pi/2)
    _y_t1 = _c_y1 + rho*sin(_psi-pi/2)

    # Get the length of the first S
    straight1_length = ( (_x_t0 - _x_t1)**2 + (_y_t0 - _y_t1)**2 )**.5

    # Get the length of the first L segment, length = theta*rho
    if _psi < thet0:
        arc1_length = abs((_psi+2*pi - thet0)*rho) # arc = theta*rho
    else:
        arc1_length = abs((_psi - thet0)*rho) # arc = theta*rho

    # Get the length of the second L segment
    # Need to divide it in 2 to account for cases where there is more than 1 rotation
    if _theta1_t < _psi:
        arc2_length_1 = abs((_theta1_t+2*pi - _psi)*rho)
    else:
        arc2_length_1 = abs((_theta1_t - _psi)*rho)
    if _theta_f < _theta1_t:
        arc2_length_2 = abs((_theta_f+2*pi - _theta1_t)*rho)
    else:
        arc2_length_2 = abs((_theta_f - _theta1_t)*rho)

    # Get the length of the last S
    _P2x = _c_x1 + rho*cos(_theta_f-pi/2)
    _P2y = _c_y1 + rho*sin(_theta_f-pi/2)
    straight2_length = ( (_P2x - xf)**2 + (_P2y - yf)**2 )**.5

    # Get the total path length
    total_path_length = arc1_length + straight1_length + arc2_length_1 + arc2_length_2 + straight2_length

    if total_path_length < distance:
    # if (i == 1): # for testing purposes
        # Only keep the parameters for the shortest path
        distance = total_path_length
        min_dist_idx = i
        c_x1 = _c_x1
        c_y1 = _c_y1
        x_t0 = _x_t0
        y_t0 = _y_t0
        x_t1 = _x_t1
        y_t1 = _y_t1
        thet1_t = _theta1_t
        psi = _psi
        theta_f = _theta_f
    # print("total_path_length:", total_path_length)

end = perf_counter() # To know how quick the code is
print("time =", end-start)

# Plot the results
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
plt.arrow(xf, yf, arrow_length*cos(theta_f), arrow_length*sin(theta_f), head_width=0.2)

# Tangent circle giving the shortest path
circle1 = plt.Circle((c_x1,c_y1), rho, fill=False, linestyle="--", color="k")
ax.add_patch(circle1)
# # Candidate circles tangent to the target
# for c1_xy in c1_candidates:
#     circle1 = plt.Circle(c1_xy, rho, fill=False, linestyle="--", color="g")
#     ax.add_patch(circle1)

# Get and plot the tangent point
p_xt, p_yt = get_circle_center(xc, yc, thet1_t, r_c, "L")
plt.plot(p_xt, p_yt, "ro")
plt.arrow(p_xt, p_yt, arrow_length*cos(thet1_t), arrow_length*sin(thet1_t), head_width=0.2)

# Plot the point that starts the first S segment
plt.plot(x_t0, y_t0, "go")

# Plot the point that ends the first S segment
plt.plot(x_t1, y_t1, "go")

# Get the points for the first curved segment, and make sure angle range is ok
if psi < thet0: # changed for L
    angles0 = np.linspace(thet0-pi, psi+pi, 20) # changed for L
    print("angles01")
else:
    angles0 = np.linspace(thet0-pi, psi-pi, 20) # changed for L
    print("angles02")
x_L0 = []
y_L0 = []
for angle in angles0:
    x, y = get_circle_center(c_x0, c_y0, angle, rho, "L") # changed for L
    x_L0.append(x)
    y_L0.append(y)
plt.plot(x_L0, y_L0, "g-")

# Get the points for the first straight segment
x_s01 = np.linspace(x_t0, x_t1, 10)
y_s01 = np.linspace(y_t0, y_t1, 10)
plt.plot(x_s01, y_s01, "g-")

# Get the points for the second curved segment, and make sure angle range is ok
if thet1_t < psi: # changed for L
    angles1_1 = np.linspace(psi-pi, thet1_t+pi, 20) # changed for L
    print("angles1")
else:
    angles1_1 = np.linspace(psi-pi, thet1_t-pi, 20) # changed for L
if theta_f < thet1_t:
    angles1_2 = np.linspace(thet1_t-pi, theta_f+pi, 20) # changed for L
else:
    angles1_2 = np.linspace(thet1_t-pi, theta_f-pi, 20) # changed for L
angles1 = angles1_1.tolist() + angles1_2.tolist()[1:]
x_L1 = []
y_L1 = []
for angle in angles1:
    x, y = get_circle_center(c_x1, c_y1, angle, rho, "L") # changed for L
    x_L1.append(x)
    y_L1.append(y)
plt.plot(x_L1, y_L1, "g-")

# Get the points for the last straight segment
x_S_f = np.linspace(x_L1[-1], xf, 10)
y_S_f = np.linspace(y_L1[-1], yf, 10)
plt.plot(x_S_f, y_S_f, "g-")

all_x_trajectory = x_L0 + x_s01.tolist()[1:-1] + x_L1 + x_S_f.tolist()[1:]
all_y_trajectory = y_L0 + y_s01.tolist()[1:-1] + y_L1 + y_S_f.tolist()[1:]
# plt.plot(all_x_trajectory, all_y_trajectory, "c-")

plt.grid(False)
plt.axis("equal")
plt.show()
