import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import Cit_par21 as par
from constants import *
from scipy.optimize import curve_fit
import CGCalc
from scipy.interpolate import interp1d
from scipy.integrate import quad

'''
THIS FILE EXTRACTS FLIGHT CHARACTERISTICS FROM FLIGHT DATA
EACH CLASS AND FUNCTION IS EXPLAINED IN THIS HEADER

Class FlightData(datafile,thrustfile,thrustfilest) - calculates all stationary data
datafile - inflight excel sheet
thrustfile(st) - (standard) output of thrust.exe 
all input files are in preflight/postflight folders - pretty self explanatory

    def(stat1) - stationary series 1
        returns alpha, cl, cd, clalpha, cd0 and e
    def(stat2) - stationary series 2
        returns Vport (reduced effective vel), fred(reduced control force
        array), dred(reduceed elevator def array), cmalpha, cmdelta

Class DynamicData(csv, excel) - calculates and plots dynamic response
csv - measurement data
excel - inflight excel sheet

    def SymMode (ShortPeriod, Phugoid)
        ShortPeriod/Phugoid - True/False, True by default
        returns 3D datamatrix in order: ShortPeriod, Phugoid
            time, velocity, angèle of attack, pitch angèle, pitch rate
        
    
    def AsymMode (DutchRoll, DutchRollYD, AperRoll, Spiral)
        True/False same as SymMode
        returns 3D datamatrix in order: DutchRoll, DutchRollYD, AperRoll, Spiral
            time, yaw angèle, roll angèle, yaw rate, roll rate

Function EigProps(data, n) - models dynamic response mathematically, calculates the coefficients
and the eigenvalue
data - flight mode
n - measurement index
returns 2D array in order: [amplitude, phase difference, natural frequency, 
        damping ratio, equilibrium point]
        [eigenvalue]

'''

class FlightData:

    def __init__(self, excel='../postflight/20200310_V2.xlsx', thrustfile='../postflight/thrustvals.dat', thrustfilest='../postflight/thrustvals_st.dat'):
        # read main excel sheet (post-flight data)
        self.excel = excel
        self.data = pd.read_excel(excel, engine='openpyxl')
        # read thrust data (normal)
        self.thrustdata = pd.read_table(thrustfile, sep="\t", names=['Tl', 'Tr'])
        self.thrust = self.thrustdata['Tl'] + self.thrustdata['Tr']
        # read thrust data (standardized)
        self.thrustdatast = pd.read_table(thrustfilest, sep="\t", names=['Tl', 'Tr'])
        self.thrustst = self.thrustdatast['Tl'] + self.thrustdatast['Tr']
        # extract pax and crew masses and fuel mass from PFD
        self.masses = self.data.iloc[6:15,7].reset_index(drop=True)
        self.fuelmass = self.data.iloc[16,3]*lbstokg

    def stat1(self):
        # extract stationary 1 data
        data1 = self.data.iloc[26:32].reset_index(drop=True)
        data1 = data1.drop(data1.columns[[0, 1, 2, 10, 11, 12]], axis=1)
        data1.columns=['hp', 'IAS', 'a', 'FFl', 'FFr', 'F. used', 'TAT']
        data1 = data1.astype(float)
        # extract the thrust data for the stationary 1 data points
        thrust1 = self.thrust.iloc[:6]

        # calculate weight based on masses and amount of fuel used
        W = (par.mdry + self.masses.sum() + self.fuelmass - data1['F. used'] * lbstokg) * g
        # calculate lift coefficient
        cl = 2 * W / (par.rho0 * par.S * (data1['IAS'] * ktstoms) ** 2)
        # calculate the change in lift coefficient per change in angèle of attack
        cla = np.mean(np.diff(cl)/np.diff(data1['a']))
        # calculate the drag coefficient
        cd = 2 * thrust1 * np.cos(np.radians(data1['a'])) / \
            (par.rho0 * par.S * (data1['IAS'] * ktstoms) ** 2)
        # sort the angèle of attack, cl, and cd values by increasing angèle of attack
        sort = data1['a'].argsort()
        cl = cl[sort]
        cd = cd[sort]        
        # the last data point in flight data is an outlier, omit it from the calculations
        if self.excel=='../postflight/20200310_V2.xlsx':
            cl = cl[:-1]
            cd = cd[:-1]
        
        coef = np.polyfit((cl**2), cd, 1)
        # calculate the zero-lift drag and the oswald efficiency factor
        cd0 = coef[1]
        e = 1/(coef[0]*par.A*np.pi)
        
        # calculate pressure from altitude
        p = par.p0*(1+par.lam*(data1['hp']*ftm)/par.Temp0)**(-g/par.lam/par.R)
        # calculate mach number
        M = np.sqrt(2/(par.gamma-1)*((1+par.p0/p*((1+(par.gamma-1)/2/par.gamma*par.rho0/par.p0* \
           (data1['IAS']*ktstoms)**2)**(par.gamma/(par.gamma-1))-1))**((par.gamma-1)/par.gamma)-1))
        # convert temperature to K and into static
        data1['TAT'] = data1['TAT'] + 273.15
        SAT = data1['TAT'] / (1+(par.gamma-1)/2*M**2)
        # viscosity at each data point
        mi = 2.791 * 10**(-7) * SAT**(0.7355)
        # air density at each data point
        rho = p/(par.R*SAT)
        # ryan reynolds number at each data point
        R = rho*data1['IAS']*par.c/mi
        # print the mach and ryan reynolds ranges
        if __name__ == "__main__":
            print('Mach range:', np.amax(M), '-', np.amin(M))
            print('Reynolds range:', np.amax(R), '-', np.amin(R))
        
        return cla,cd0,e

    def stat2(self):
        # extract the stationary 2 data
        data2 = self.data.iloc[57:64].reset_index(drop=True)
        data2 = data2.drop(data2.columns[[0, 1, 2]], axis=1)
        data2.columns=['hp', 'IAS', 'a', 'de', 'detr', 'Fe', 'FFl', 'FFr', 'F. used', 'TAT']
        data2 = data2.astype(float)
        # convert temperature to kelvins
        data2['TAT'] = data2['TAT'] + 273.15
        # extract the thrust data for the stationary 2 data points
        thrust2 = self.thrust[6:].reset_index(drop=True)
        thrust2st = self.thrustst[6:].reset_index(drop=True)     

        # calculate weight based on masses and amount of fuel used
        W = (par.mdry + self.masses.sum() + self.fuelmass - data2['F. used'] * lbstokg) * g
        # calculate pressure from altitude
        p = par.p0*(1+par.lam*(data2['hp']*ftm)/par.Temp0)**(-g/par.lam/par.R)
        # calculate mach number
        M = np.sqrt(2/(par.gamma-1)*((1+par.p0/p*((1+(par.gamma-1)/2/par.gamma*par.rho0/par.p0* \
           (data2['IAS']*ktstoms)**2)**(par.gamma/(par.gamma-1))-1))**((par.gamma-1)/par.gamma)-1))
        # convert total temperature into static temperature
        SAT = data2['TAT'] / (1+(par.gamma-1)/2*M**2)
        # calculate speed of sound
        a = np.sqrt(par.gamma * par.R * SAT)
        # calculate true airspeed
        Vtrue = a * M
        # calculate equivalent airspeed
        Ve = Vtrue * np.sqrt((p / par.R / SAT) / par.rho0)
        # calculate bruno fernandes
        Vport = Ve * np.sqrt(par.Wst / W)
        # calculate his teammate fred
        fred = data2['Fe']*par.Wst/W
        
        ##### CG PART #####
        # extract the cg shift data
        datacg = self.data.iloc[73:75].reset_index(drop=True)
        datacg = datacg.drop(datacg.columns[[0, 1, 2]], axis=1)
        datacg.columns=['hp', 'IAS', 'a', 'de', 'detr', 'Fe', 'FFl', 'FFr', 'F. used', 'TAT']
        datacg = datacg.astype(float)
        # calculate weights at cg shift locations
        Wcg = (par.mdry + self.masses.sum() + self.fuelmass - datacg['F. used'] * lbstokg) * g
        #calculate normal force coefficient
        cn = np.mean(2 * Wcg / (par.rho0 * par.S * (datacg['IAS'] * ktstoms) ** 2))
        # calculate the change in elevator deflection
        Dde = np.radians(datacg['de'][0] - datacg['de'][1])
        # import the change in the centre of gravity
        Dcg = CGCalc.CG(self.excel)[1]
        # calculate the change in moment coefficient per change in elevator deflection
        cmd = -cn/Dde * Dcg/par.c
        ##### END CG PART #####
        
        # calculate the change in elevator deflection per change in angèle of attack
        ddeda = np.mean(np.diff(data2['de'])/np.diff(data2['a']))
        # calculate the change in moment coefficient per change on angèle of attack
        cma = -cmd * ddeda
        # calculate standard and real thrust coefficients
        tcs = thrust2st / (par.rho0 * par.Ae * (data2['IAS'] * ktstoms) ** 2)
        tc = thrust2 / (par.rho0 * par.Ae * (data2['IAS'] * ktstoms) ** 2)
        # calculate fred's cousin dred
        dred = -data2['de']/cmd * par.Cmtc*(tcs-tc)
        # sort the Vport, dred, and fred arrays based on Vport (ascending)
        sort = Vport.argsort()
        Vport = Vport[sort]
        dred = dred[sort]
        fred = fred[sort]
        
        return Vport.to_numpy(),fred.to_numpy(),dred.to_numpy(),cma,cmd

class DynamicData:
    
    def __init__(self, csv='../postflight/flightdata.csv', excel='../postflight/20200310_V2.xlsx'):
        # read the flight sheet to get timestamps for the eigenmotions
        self.alltimes = pd.read_excel(excel, engine='openpyxl')
        self.eigm = self.alltimes.iloc[81:].reset_index(drop=True)
        self.eigm = self.eigm.drop(self.eigm.columns[[1, 2, 5, 8, 10, 11, 12]], axis=1)
        self.eigm.columns = ['eigenmotion', 'timestamp'] * 3
        self.eigm = pd.concat([self.eigm.iloc[:,[2*n, 2*n+1]] for n in range(3)], ignore_index=True)
        # convert timestamps into seconds
        self.eigm['time [sec]'] = [float(str(self.eigm.iloc[:, 1][i])[:2]) * 3600 + \
                              float(str(self.eigm.iloc[:, 1][i])[3:5]) * 60 + \
                              float(str(self.eigm.iloc[:, 1][i])[6:]) for i in range(6)]

        # read the recorded flight data file
        self.data = pd.read_csv(csv)    
        f= open('../DSpace Parameters.txt','r')
        title = []
        for line in f.readlines():
            title.append(line.rstrip('\n'))
        self.data.columns = title
        if excel == '../postflight/20200310_V2.xlsx':
            self.color = flight_color
            self.style = flight_style
        elif excel == '../preflight/PFD_11-03-2021fl1.xlsx':
            self.color = ref_color
            self.style = ref_style    
    # symmetric modes
    def SymMode(self, ShortPeriod=True, Phugoid=True, plotting=False):
        # create empty list to store useful data depending on which modes are active
        symmode = []
        titles = []
        
        # make a function to record the time, velocity, angèle of attack, pitch
        # angèle and pitch rate for the duration of the motion
        def MakeArraySM(t0,te):
            # set the duration of motion
            ran = (self.data['time'] >= t0) & (self.data['time'] <= te)
            # store all values in an array of 5 rows, reset time index
            arr = np.array(self.data['time'][ran]-t0)
            arr = np.append([arr], [self.data['Dadc1_tas'][ran]* ktstoms], axis=0)
            arr = np.append(arr, [self.data['vane_AOA'][ran]], axis=0)
            arr = np.append(arr, [self.data['Ahrs1_Pitch'][ran]], axis=0)
            arr = np.append(arr, [self.data['Ahrs1_bPitchRate'][ran]], axis=0)
            return arr

        if ShortPeriod:
            # get initial time from excel sheet
            t0 = self.eigm['time [sec]'][1]
            # set the end time
            te = t0+10
            sharr = MakeArraySM(t0,te)
            symmode.append(sharr)
            titles.append('Short Period')
        if Phugoid:
            # get initial time from excel sheet
            t0 = self.eigm['time [sec]'][0]
            # set the end time
            te = t0+120
            pharr = MakeArraySM(t0,te)
            symmode.append(pharr)
            titles.append('Phugoid')

        # plot: velocity, angèle of attack, pitch angèle, pitch rate            
        if plotting:    
            for i in range(len(symmode)):
                fig, axs = plt.subplots(2, 2, sharex=True)
                plt.suptitle(titles[i], size='x-large', weight='bold')
                axs[0, 0].plot(symmode[i][0], symmode[i][1], c=self.color, linestyle=self.style)
                axs[0, 0].set_title('Velocity [m/s]')
                axs[0, 1].plot(symmode[i][0], symmode[i][2], c=self.color, linestyle=self.style)
                axs[0, 1].set_title('Angle of Attack [deg]')
                axs[1, 0].plot(symmode[i][0], symmode[i][3], c=self.color, linestyle=self.style)
                axs[1, 0].set_title('Pitch Angle [deg]')
                axs[1, 1].plot(symmode[i][0], symmode[i][4], c=self.color, linestyle=self.style)
                axs[1, 1].set_title('Pitch Rate [deg/s]')
                [ax.grid() for ax in axs.flatten()]
                plt.show()
            
        return symmode
            
    # asymmetric modes
    def AsymMode(self, DutchRoll=True, DutchRollYD=True, AperRoll=True, Spiral=True, plotting=False):
        # create empty list to store useful data depending on which modes are active
        asymmode = []
        titles = []
        
        # make a function to record the time, veloctiy, yaw angèle (computed),
        # roll angèle, roll rate and yaw rate for the duration of the motion
        def MakeArrayAM(t0,te):
            # set the duration of motion
            ran = (self.data['time'] >= t0) & (self.data['time'] <= te)
            # compute the yaw angèle based on its derivative, the yaw rate

            time = self.data['time'][ran]-t0
            yaw_rate = self.data['Ahrs1_bYawRate'][ran]
            interp = interp1d(time, yaw_rate, kind='cubic')
            Yaw = np.zeros(len(self.data['time'][ran]))

            dt = 0.1
            for i in range(1,len(Yaw)):
                Yaw[i] = Yaw[i-1] +  quad(interp, dt*(i-1), dt*i)[0]

            # store all values in an array of 5 rows, reset time index
            arr = np.array(self.data['time'][ran]-t0)
            arr = np.append([arr], [self.data['Dadc1_tas'][ran] * ktstoms], axis=0)
            arr = np.append(arr, [Yaw], axis=0)
            arr = np.append(arr, [self.data['Ahrs1_Roll'][ran]], axis=0)
            arr = np.append(arr, [self.data['Ahrs1_bRollRate'][ran]], axis=0)
            arr = np.append(arr, [self.data['Ahrs1_bYawRate'][ran]], axis=0)
            return arr

        if DutchRoll:
            # get initial time from excel sheet
            t0 = self.eigm['time [sec]'][2]
            # set the end time
            te = t0+15
            drarr = MakeArrayAM(t0,te)
            asymmode.append(drarr)
            titles.append('Dutch Roll')
        if DutchRollYD:
            # get initial time from excel sheet
            t0 = self.eigm['time [sec]'][3]
            # set the end time
            te = t0+10
            drydarr = MakeArrayAM(t0,te)
            asymmode.append(drydarr)
            titles.append('Dutch Roll with Yaw Damper')
        if AperRoll:
            # get initial time from excel sheet
            t0 = self.eigm['time [sec]'][4]
            # set the end time
            te = t0+15
            aparr = MakeArrayAM(t0,te)
            asymmode.append(aparr)
            titles.append('Aperiodic Roll')
        if Spiral:
            # get initial time from excel sheet
            t0 = self.eigm['time [sec]'][5]
            # set the end time
            te = t0+50
            sparr = MakeArrayAM(t0,te)
            asymmode.append(sparr)
            titles.append('Spiral')
    
        # plot: yaw angèle, roll angèle, yaw rate, roll rate
        if plotting:
            for i in range(len(asymmode)):
                fig, axs = plt.subplots(2, 2, sharex=True)
                plt.suptitle(titles[i], size='x-large', weight='bold')
                axs[0, 0].plot(asymmode[i][0], asymmode[i][2], c=self.color, linestyle=self.style)
                axs[0, 0].set_title('Yaw Angle [deg]')
                axs[0, 1].plot(asymmode[i][0], asymmode[i][3], c=self.color, linestyle=self.style)
                axs[0, 1].set_title('Roll Angle [deg]')
                axs[1, 0].plot(asymmode[i][0], asymmode[i][4], c=self.color, linestyle=self.style)
                axs[1, 0].set_title('Yaw Rate [deg/s]')
                axs[1, 1].plot(asymmode[i][0], asymmode[i][5], c=self.color, linestyle=self.style)
                axs[1, 1].set_title('Roll Rate [deg/s]')
                [ax.grid() for ax in axs.flatten()]
                plt.show()
            
        return asymmode

# define function to mathematically model the given flight mode (data) and measurement (n)
# e.g. data=symmode[1] for phugoid and n=4 for pitch rate
# AND get the related eigenvalue
def EigProps(data,n):
    # define function to model the flight mode measurements
    def func(t, A,phi,f,l,v):
    # A->amplitude, phi->phase difference, f->natural frequency, l->damping ratio, v->equilibrium point
        func = A*np.sin((t+phi)*2*np.pi*f)*np.e**(-l*t)+v
        return func
    # get appropriate values for A,phi,f,l,v from curve fit
    popt, _ = curve_fit(func, data[0], data[n])
    # calculate the real and imaginary parts of the eigenvalue
    a = -popt[3]*par.c/data[1][0]
    b = 2*np.pi*popt[2]*par.c/data[1][0]
    # define teh eigenvalue
    eigenvalue = np.complex(a,b)

    plt.plot(data[0],data[n])
    plt.plot(data[0],func(data[0],*popt))
    plt.show()
    # return the coefficients
    return popt, eigenvalue
def Trim(data,ts,te):
    trimmed = [data[n][(np.where(data[0]<te)) and (np.where(data[0]>ts))] for n in range(len(data))]
    return trimmed


if __name__ == "__main__":
    #https://www.youtube.com/watch?v=XFE02ogS17U
    # if you want to use reference data,,,
    # fd = FlightData(excel='../preflight/PFD_11-03-2021fl1.xlsx', thrustfile='../preflight/thrustvals.dat', thrustfilest='../preflight/thrustvals_st.dat')
    # dd = DynamicData(excel='../preflight/PFD_11-03-2021fl1.xlsx')
    
    # fd = FlightData()
    # stationary1 = fd.stat1()
    # stationary2 = fd.stat2()

    dd = DynamicData()

    symdata = dd.SymMode(1,1,0)
    asymdata = dd.AsymMode(1,1,1,1,0)
    print(EigProps(Trim(symdata[1],35,75),4))
    # EigProps(asymdata[0],3)