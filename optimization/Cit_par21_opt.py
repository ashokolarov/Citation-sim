# Citation 550 - Linear simulation

from constants import *
import numpy as np

eigen_labels = ['ph', 'sp', 'dr', 'dryd', 'apr', 'spi']

Cma, CZa, Cmq, Cmde, Cmu, CXu, CZu, KY2, CXa, CXq, CZq, CXde, CZde, CZadot, Cmadot  = [-0.39579702, -3.64301832, -8.71073962, -1.4685174,  0.03013532,
       -0.09535426, -0.78215698,  1.49040473,  0.20831055, -0.28246101,
       -5.62894729,  0.07456   , -1.79125959, -0.014     ,  0.26224654]

Clp, KX2, CYb, Cnb, Cnr, KZ2, Clb, Clr, CYdr, Clda, Cnda, Cndr, CYr, CYda, Cldr, Cnp, CXadot, CYp   = [-1.05381431,  0.03466663,  0.32332385,  0.14254302, -0.3156647 ,
        0.06323342, -0.18120948,  0.33806808, -0.46      , -0.46882385,
       -0.00808291, -0.0879196 ,  0.8495    , -0.16      ,  0.03388838,
       -0.0602    ,  0.0833    , -0.0304    ]

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
Ae     = 0.21

# Constant values concerning atmosphere and gravity

rho0   = 1.2250          # air density at sea level [kg/m^3] 
lam = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)
gamma = 1.4  # specific heat ratio for dry air standard conditions
p0 = rho0 * R * Temp0
Mair = 28.97 #grams/mole


W = mtot * g            # [N]       (aircraft weight)

# Constant values concerning aircraft inertia
KXZ    = 0.002

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * np.pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Stability derivatives



CYbdot =  0
Cnbdot =   0

