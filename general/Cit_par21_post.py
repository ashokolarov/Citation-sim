# Citation 550 - Linear simulation

from constants import *
import numpy as np

eigen_labels = ['ph', 'sp', 'dr', 'dryd', 'apr', 'spi']


# Aicraft Standard Values
Wst = 60500  # [N] Standard aircraft mass
ffst = 0.048  # [kg/s] Standard Fuel Flow PER ENGINE
rhost = 1.225  # [kg/m3] standard air density

# Aircraft mass
mdry = 9165 * lbstokg  # dry mass [kg]
mfuel = 750 * lbstokg  # fuel mass [kg]
p1m = 82  # mass passenger 1 [kg]
p2m = 66  # mass passenger 2 [kg]
p3m = 81  # mass passenger 3 [kg]
p4m = 69  # mass passenger 4 [kg]
p5m = 85  # mass passenger 5 [kg]
p6m = 96  # mass passenger 6 [kg]
crew = 95 + 102 + 89  # mass crew (2 pilots + coordinator) [kg]
mtot = mdry + mfuel + p1m + p2m + p3m + p4m + p5m + p6m + crew  # total mass [kg]


# Aerodynamic properties
e = 0.8  # Oswald factor [ ]
CD0 = 0.04  # Zero lift drag coefficient [ ]
CLa = 5.084  # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma  = -0.5626  # longitudinal stabilty [ ]
Cmde = -1.1642  # elevator effectiveness [ ]
Cmtc = -0.0064  # thrust moment arm []

# Aircraft geometry

S      = 30.00	          # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 0.71 * 5.968    # tail length [m]
c      = 2.0569	          # mean aerodynamic cord [m]
xcg    = 0.25 * c
lh_c   = lh / c	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabiliser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabiliser aspect ratio [ ]
Vh_V   = 1	          # [ ]
ih     = -2 * np.pi / 180   # stabiliser angle of incidence [rad]

# Constant values concerning atmosphere and gravity

rho0   = 1.2250          # air density at sea level [kg/m^3] 
lam = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)
gamma = 1.4  # specific heat ratio for dry air standard conditions
p0 = rho0 * R * Temp0

W = mtot * g            # [N]       (aircraft weight)

# Constant values concerning aircraft inertia
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 1.25 * 1.114

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * np.pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Stability derivatives

CXu    = -0.09500
CXa    = +0.47966		# Positive! (see FD lecture notes) 
CXadot = +0.08330
CXq    = -0.28170
CXde   = -0.03728

CZu    = -0.37616
CZa    = -5.74340
CZadot = -0.00350
CZq    = -5.66290
CZde   = -0.69612

Cmu    = +0.06990
Cmadot = +0.17800
Cmq    = -8.79415

CYb    = -0.7500
CYbdot =  0     
CYp    = -0.0304
CYr    = +0.8495
CYda   = -0.0400
CYdr   = +0.2300

Clb    = -0.10260
Clp    = -0.71085
Clr    = +0.23760
Clda   = -0.23088
Cldr   = +0.03440

Cnb    =  +0.1348
Cnbdot =   0     
Cnp    =  -0.0602
Cnr    =  -0.2061
Cnda   =  -0.0120
Cndr   =  -0.0939
