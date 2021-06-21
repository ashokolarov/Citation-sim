import pandas as pd
import numpy as np
from constants import *
import Cit_par21 as par
import matplotlib.pyplot as plt

def generate_thrust_data(filename, mode):
    data = pd.read_excel(filename, engine='openpyxl')
    
    stat1 = data.iloc[26:32]
    stat1 = stat1.reset_index(drop=True)
    stat1 = stat1.drop(stat1.columns[[0, 1, 2, 10, 11, 12]], axis=1)
    stat1.columns=['hp', 'IAS', 'a', 'FFl', 'FFr', 'F. used', 'TAT']
    
    trim = data.iloc[57:64]
    trim = trim.reset_index(drop=True)
    trim = trim.drop(trim.columns[[0, 1, 2]], axis=1)
    trim.columns=['hp', 'IAS', 'a', 'de', 'detr', 'Fe', 'FFl', 'FFr', 'F. used', 'TAT']
    
    df = stat1.append(trim)
    df = df.reset_index(drop=True)
    df = df.astype(float)

    hp = df['hp'] * ftm
    p = par.p0*(1+par.lam*hp/par.Temp0)**(-g/par.lam/par.R)
    M = np.sqrt(2/(par.gamma-1)*((1+par.p0/p*((1+(par.gamma-1)/2/par.gamma*par.rho0/par.p0*
       (df['IAS']*ktstoms)**2)**(par.gamma/(par.gamma-1))-1))**((par.gamma-1)/par.gamma)-1))
    
    if mode == 0:
        mf1 = df['FFl'] * fuelfac
        mf2 = df['FFr'] * fuelfac
        output = "../preflight/thrust.dat"
    elif mode == 1:
        mf1 = par.ffst/lbstokg
        mf2 = par.ffst/lbstokg
        output = "../preflight/thrust_st.dat"
    
    df['TAT'] = df['TAT'] + 273.15
    SAT = df['TAT'] / (1+(par.gamma-1)/2*M**2)
    T_isa = np.array([par.Temp0 + par.lam * h for h in hp])
    Td = SAT - T_isa

    
    with open(output, 'w+') as file:
        for i in range(len(hp)):
            if mode == 0:
                line = str(hp[i])+' '+str(M[i])+' '+str(Td[i])+' '+str(mf1[i])+' '+str(mf2[i])+'\n'
            elif mode == 1:
                line = str(hp[i])+' '+str(M[i])+' '+str(Td[i])+' '+str(mf1)+' '+str(mf2)+'\n'
            file.write(line)

if __name__ == "__main__":
    filename = "../preflight/PFD_11-03-2021fl1.xlsx"
    #filename = "../postflight/20200310_V2.xlsx" #POSTFLIGHT
    #mode 0 is for actual data, 1 for standardised data
    mode = 1
    generate_thrust_data(filename, mode)