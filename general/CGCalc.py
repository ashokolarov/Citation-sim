import numpy as np
import pandas as pd
import Cit_par21 as par
import constants as c

# FUNCTION TO CALCULATE THE CENTER OF GRAVITY AT ALL DATA POINTS
def CG(excel='../postflight/20200310_V2.xlsx'):
    # read the excel sheet
    data = pd.read_excel(excel, engine='openpyxl')
    # basic empty mass
    BEM = par.mdry/c.lbstokg
    # cg for bem
    bemcg = 291.648
    # read the masses of the passengers
    masses = data.iloc[6:15,7].reset_index(drop=True).to_numpy()/c.lbstokg
    # the moment arms for each passenger seat
    payloadarms = np.array([131,131,170,214,214,251,251,288,288]) # [in]
    # calculate the total moment caused by the payload
    payloadmoment = 0
    for i in range(len(payloadarms)):
        payloadmoment += masses[i]*payloadarms[i]
    # zero fuel mass cg
    zfmcg = (payloadmoment+BEM*bemcg)/(np.sum(masses)+BEM)
    
    # read block fuel
    fuelblock = data.iloc[16,3]
    # function to calculate the moment caused by the fuel for fuel mass
    def FuelMoment(fuelmass):
        return (2.8526*fuelmass+9.8957)*100
    
    # read the data for fuel used at each data point
    stat1 = data.iloc[26:32,8]
    trim = data.iloc[57:64,11]
    cgshift = data.iloc[74,11]  
    fused = np.append(pd.concat([stat1,trim]).reset_index(drop=True).to_numpy(),cgshift)
    
    # calculate the ramp mass and cg
    rampcg = (payloadmoment + BEM*bemcg + FuelMoment(fuelblock)) /\
                (np.sum(masses) + BEM + fuelblock)
    rampmass = (np.sum(masses) + BEM + fuelblock)
    cgs = [rampcg]
    ms = [rampmass]
    # calculate the mass and cg for all data points and store them in a list
    [cgs.append((payloadmoment + BEM*bemcg + FuelMoment(fuelblock-fused[i])) /\
                (np.sum(masses) + BEM + fuelblock)) for i in range(len(fused)-1)]
    [ms.append((np.sum(masses) + BEM + fuelblock-fused[i])) for i in range(len(fused))]
    # calculate the cg after the passenger at 3R moves front
    shiftcg = (payloadmoment + (170-payloadarms[-1])*masses[-1] + BEM*bemcg + \
               FuelMoment(fuelblock-fused[-1])) / (np.sum(masses) + BEM + fuelblock)
    cgs.append(shiftcg)
    # return the all the cgs, the shift in cg caused by the passenger movement, and the masses
    return cgs, (cgs[-2] - cgs[-1])*c.intom, ms

if __name__ == "__main__":
    cgs, cgshift, acmasses = CG()
    df = pd.DataFrame({"Center of gravity [in]":cgs, "Mass [lbs]":acmasses})
    df.to_csv("cgs&masses.csv")