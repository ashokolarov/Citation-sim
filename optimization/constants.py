import numpy as np
import matplotlib.pyplot as plt
plt.style.use('classic')

r2d = 180 / np.pi               # rad to deg[-]
d2r = np.pi / 180               # deg to rad[-]
g = 9.80665                     # gravitational acceleration[m/s^2]
lbstokg = 0.453592              # lbs to kg[-]
ftm = 0.3048                    # feet to meter[-]
fuelfac = lbstokg / 60 / 60     # lbs/hr to kg/s
ktstoms = 0.514444              # knots to meters per second
intom = 0.0254                  # inches to meters

ref_color = 'g'
num_color = 'royalblue'
flight_color = 'r'

ref_style = (0, (3, 1, 1, 1, 1, 1))
num_style = '-'
flight_style = '--'