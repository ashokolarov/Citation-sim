import numpy as np
import control as ctr
import matplotlib.pyplot as plt
from constants import *
import pandas as pd
import Cit_par21_post as param


class Model:

    def __init__(self, datafile='../postflight/flightdata.csv',
                 excel_file='../postflight/20200310_V2.xlsx', dt=0.1):

        self.data = pd.read_csv(datafile)
        with open("../DSpace Parameters.txt") as namesfile:
            names = namesfile.read().splitlines()
        self.data.columns = names

        alltimes = pd.read_excel(excel_file, engine='openpyxl')
        eigm = alltimes.iloc[81:].reset_index(drop=True)
        eigm = eigm.drop(eigm.columns[[1, 2, 5, 8, 10, 11, 12]], axis=1)
        eigm.columns = ['eigenmotion', 'timestamp'] * 3
        eigm = pd.concat([eigm.iloc[:, [2 * n, 2 * n + 1]] for n in range(3)], ignore_index=True)
        # convert timestamps into seconds
        eigm['time [sec]'] = [float(str(eigm.iloc[:, 1][i])[:2]) * 3600 +
                              float(str(eigm.iloc[:, 1][i])[3:5]) * 60 +
                              float(str(eigm.iloc[:, 1][i])[6:]) for i in range(6)]

        timestamps = {param.eigen_labels[i]: eigm['time [sec]'][i] for i in range(len(param.eigen_labels))}
        self.timestamps = timestamps

        input_bandwith = [119.9, 9.99, 14.99, 10, 14.99, 49.99]

        self.states = {key: extract_state(self.data, timestamps[key], timestamps) for key in param.eigen_labels}
        self.inputs = {key: extract_input(self.data, timestamps[key], timestamps, input_bandwith[i]) for i, key in
                       enumerate(param.eigen_labels)}

        self.dt = dt
        self.labels = param.eigen_labels
        self.motion_names = ['Phugoid', 'Short period', 'Dutch roll', 'Dutch roll with yaw damper', 'Aperiodical roll',
                             'Spiral']
        self.state_space_systems = {}
        self.fill_state_space_dict()

    def symmetric(self, V0, CX0, CZ0, muc):
        C1 = np.array([[-2 * muc * (param.c / V0), 0, 0, 0],
                       [0, (param.CZadot - 2 * muc) * (param.c / V0), 0, 0],
                       [0, 0, -param.c / V0, 0],
                       [0, param.Cmadot * (param.c / V0), 0,
                        -2 * muc * param.KY2 * (param.c / V0)]], dtype=np.float64)

        C2 = np.array([[param.CXu, param.CXa, CZ0, param.CXq],
                       [param.CZu, param.CZa, -CX0, param.CZq + 2 * muc],
                       [0, 0, 0, 1],
                       [param.Cmu, param.Cma, 0, param.Cmq]], dtype=np.float64)

        C3 = np.array([param.CXde, param.CZde, 0, param.Cmde], dtype=np.float64).reshape(4, 1)

        A = -np.linalg.inv(C1) @ C2
        B = -np.linalg.inv(C1) @ C3
        C = np.eye(4)
        D = np.zeros(4).reshape(4, 1)

        return A, B, C, D

    def fill_state_space_dict(self):
        for i in range(0, 2):
            state = self.states[param.eigen_labels[i]]
            V0 = state['V0']
            CX0 = state['CX0']
            CZ0 = state['CZ0']
            muc = state['muc']
            A, B, C, D = self.symmetric(V0, CX0, CZ0, muc)
            self.state_space_systems[param.eigen_labels[i]] = ctr.ss(A, B, C, D)

        for i in range(2, 6):
            state = self.states[param.eigen_labels[i]]
            V0 = state['V0']
            CL = state['CL']
            mub = state['mub']
            A, B, C, D = self.asymmetric(V0, CL, mub)
            self.state_space_systems[param.eigen_labels[i]] = ctr.ss(A, B, C, D)

    def symmetric_nonzero_response(self, plot_eigs=False):
        ss = self.state_space_systems['ph']

        print(f"\t\t{self.motion_names[1]}")
        fdp = ctr.damp(ss)
        print()

        if plot_eigs:
            eigs_num = fdp[2]
            X = [eig.real for eig in eigs_num]
            Y = [eig.imag for eig in eigs_num]
            plt.scatter(X, Y, color='red', label='Numerical')
            plt.grid(True, which='both')
            plt.axhline(y=0, color='k')
            plt.axvline(x=0, color='k')
            plt.xlabel('Real axis', fontsize=20)
            plt.ylabel('Imaginary axis', fontsize=20)
            #plt.title('System eigenvalues in symmetrical flight', size='x-large', weight='bold')
            plt.show()

        alpha0 = 2 * d2r
        X0 = np.array([0, alpha0, 0, 0], dtype=np.float64).reshape(4, 1)
        t = np.arange(0, 15, 0.01)
        _, y = ctr.initial_response(ss, t, X0)

        V0 = self.states['sp']['V0']

        y[0] = y[0] * V0 + V0
        y[1] *= r2d
        y[2] *= r2d
        y[3] = y[3] * r2d * V0 / param.c

        fig, axs = plt.subplots(2, 2, sharex='col')

        axs[0, 0].ticklabel_format(useOffset=False)
        axs[0, 0].ticklabel_format(style='plain')
        axs[0, 0].plot(t, y[0], color=num_color, linewidth=2, linestyle=num_style)
        axs[0, 0].set_ylabel('Velocity[m/s]', fontsize=16)
        axs[0, 0].grid(True)

        axs[0, 1].plot(t, y[1], color=num_color, linewidth=2, linestyle=num_style)
        axs[0, 1].set_ylabel(r'Angle of attack $\alpha$[deg]', fontsize=16)
        axs[0, 1].grid(True)

        axs[1, 0].plot(t, y[2], color=num_color, linewidth=2, linestyle=num_style)
        axs[1, 0].set_xlabel('Time[s]', fontsize=16)
        axs[1, 0].set_ylabel(r'Pitch $\theta$[deg]', fontsize=16)
        axs[1, 0].grid(True)

        axs[1, 1].plot(t, y[3], color=num_color, linewidth=2, linestyle=num_style)
        axs[1, 1].set_xlabel('Time[s]', fontsize=16)
        axs[1, 1].set_ylabel('Pitch rate q[deg/s]', fontsize=16)
        axs[1, 1].grid(True)

        fig.suptitle('Symmetrical flight\n' r'Response to an initial condition of α=2[deg]', size='x-large',
                     weight='bold')
        plt.show()

    def response_symmetric_mode(self, label, Tf, plot=False):
        idx = self.labels.index(label)
        state = self.states[self.labels[idx]]
        V0 = state['V0']
        initial = state['X0']

        X0 = np.array([0, 0, 0, initial[-1]]).reshape(4, 1)
        ss = self.state_space_systems[label]

        t = np.arange(0, Tf, self.dt)

        u = self.inputs[self.labels[idx]]
        u = np.concatenate((u, np.ones(len(t) - len(u)) * u[0]), axis=None)
        _, y, x = ctr.forced_response(ss, t, u, X0)

        y[0] = y[0] * V0 + V0
        y[1] += initial[1]
        y[1] *= r2d
        y[2] += initial[2]
        y[2] *= r2d
        y[3] = y[3] * r2d * V0 / param.c

        if plot:
            fig, axs = plt.subplots(2, 2, sharex='col')

            axs[0, 0].ticklabel_format(useOffset=False)
            axs[0, 0].ticklabel_format(style='plain')
            axs[0, 0].plot(t, y[0], color=num_color, linewidth=2, linestyle=num_style)
            axs[0, 0].set_ylabel('Velocity[m/s]')
            axs[0, 0].grid(True)

            axs[0, 1].plot(t, y[1], color=num_color, linewidth=2, linestyle=num_style)
            axs[0, 1].set_ylabel(r'Angle of attack $\alpha$[deg]')
            axs[0, 1].grid(True)

            axs[1, 0].plot(t, y[2], color=num_color, linewidth=2, linestyle=num_style)
            axs[1, 0].set_xlabel('Time [s]')
            axs[1, 0].set_ylabel(r'Pitch $\theta$[deg]')
            axs[1, 0].grid(True)

            axs[1, 1].plot(t, y[3], color=num_color, linewidth=2, linestyle=num_style)
            axs[1, 1].set_xlabel('Time [s]')
            axs[1, 1].set_ylabel('Pitch rate q[deg/s]')
            axs[1, 1].grid(True)

            fig.suptitle(f'{self.motion_names[idx]}', size='x-large', weight='bold')
            plt.show()

        return t, y

    def asymmetric(self, V0, CL, mub):
        C1 = np.array([[(param.CYbdot - 2 * mub) * param.b / V0, 0, 0, 0],
                       [0, -param.b / (2 * V0), 0, 0],
                       [0, 0, -4 * mub * param.KX2 * param.b / V0,
                        4 * mub * param.KXZ * param.b / V0],
                       [param.Cnbdot * param.b / V0, 0, 4 * mub * param.KXZ * param.b / V0,
                        -4 * mub * param.KZ2 * param.b / V0]], dtype=np.float64)

        C2 = np.array([[param.CYb, CL, param.CYp, param.CYr - 4 * mub],
                       [0, 0, 1, 0],
                       [param.Clb, 0, param.Clp, param.Clr],
                       [param.Cnb, 0, param.Cnp, param.Cnr]], dtype=np.float64)

        C3 = np.array([[param.CYda, param.CYdr],
                       [0, 0],
                       [param.Clda, param.Cldr],
                       [param.Cnda, param.Cndr]], dtype=np.float64)

        A = -np.linalg.inv(C1) @ C2
        B = -np.linalg.inv(C1) @ C3
        return A, B, np.eye(4), np.zeros((4, 2))

    def asymmetric_nonzero_response(self, plot_eigs=False):
        ss = self.state_space_systems['spi']

        print(f"\t\t{self.motion_names[2]}")
        fdp = ctr.damp(ss)
        print()

        if plot_eigs:
            eigs_num = fdp[2]
            X = [eig.real for eig in eigs_num]
            Y = [eig.imag for eig in eigs_num]
            plt.scatter(X, Y, color='red', label='Numerical')
            plt.grid(True, which='both')
            plt.axhline(y=0, color='k')
            plt.axvline(x=0, color='k')
            plt.xlabel('Real axis', fontsize=20)
            plt.ylabel('Imaginary axis', fontsize=20)
            #plt.title('System eigenvalues in symmetrical flight', size='x-large', weight='bold')
            plt.show()

        V0 = self.states['dr']['V0']
        yaw0 = 5 * d2r
        X0 = np.array([yaw0, 0, 0, 0], dtype=np.float64).reshape(4, 1)
        t = np.arange(0, 15, 0.01)
        _, y = ctr.initial_response(ss, t, X0)

        y[0] *= r2d
        y[1] *= r2d
        y[2] = y[2] * 2 * V0 / param.b
        y[3] = y[3] * 2 * V0 / param.b

        fig, axs = plt.subplots(2, 2, sharex=True)

        axs[0, 0].plot(t, y[0], color=num_color, linewidth=2, linestyle=num_style)
        axs[0, 0].set_ylabel(fr'Sideslip $\beta$[deg]')
        axs[0, 0].grid(True)

        axs[0, 1].plot(t, y[3], color=num_color, linewidth=2, linestyle=num_style)
        axs[0, 1].set_ylabel('Yaw rate r[deg/s]')
        axs[0, 1].grid(True)

        axs[1, 0].plot(t, y[1], color=num_color, linewidth=2, linestyle=num_style)
        axs[1, 0].set_ylabel('Roll angle φ[deg]')
        axs[1, 0].set_xlabel('Time [s]')
        axs[1, 0].grid(True)

        axs[1, 1].plot(t, y[2], color=num_color, linewidth=2, linestyle=num_style)
        axs[1, 1].set_title('Roll rate p[deg/s]')
        axs[1, 1].set_xlabel('Time [s]')
        axs[1, 1].grid(True)

        fig.suptitle('Asymmetrical flight\n' r'Response to an initial condition of β=5[deg]', size='x-large',
                     weight='bold')
        plt.show()

    def response_asymmetric_mode(self, label, Tf, plot=False):
        idx = self.labels.index(label)
        state = self.states[self.labels[idx]]
        V0 = state['V0']
        X0 = state['X0']
        ss = self.state_space_systems[label]

        t = np.arange(0, Tf, self.dt)

        u = self.inputs[self.labels[idx]]

        u1 = np.concatenate((u[0], np.ones(len(t) - len(u[0])) * u[0][0]), axis=None)
        u2 = np.concatenate((u[1], np.ones(len(t) - len(u[1])) * u[1][0]), axis=None)

        u = np.array([u1, u2])

        _, y, x = ctr.forced_response(ss, t, u, X0)

        y[0] *= -r2d
        y[1] *= -r2d
        y[2] = -y[2] * r2d * 2 * V0 / param.b
        y[3] = -y[3] * r2d * 2 * V0 / param.b

        if plot:
            fig, axs = plt.subplots(2, 2, sharex=True)

            axs[0, 0].plot(t, y[0], color=num_color, linewidth=2, linestyle=num_style)
            axs[0, 0].set_ylabel(fr'Sideslip $\beta$[deg]')
            axs[0, 0].grid(True)

            axs[0, 1].plot(t, y[3], color=num_color, linewidth=2, linestyle=num_style)
            axs[0, 1].set_ylabel('Yaw rate r[deg/s]')
            axs[0, 1].grid(True)

            axs[1, 0].plot(t, y[1], color=num_color, linewidth=2, linestyle=num_style)
            axs[1, 0].set_ylabel('Roll angle φ[deg]')
            axs[1, 0].set_xlabel('Time [s]')
            axs[1, 0].grid(True)

            axs[1, 1].plot(t, y[2], color=num_color, linewidth=2, linestyle=num_style)
            axs[1, 1].set_title('Roll rate p[deg/s]')
            axs[1, 1].set_xlabel('Time [s]')
            axs[1, 1].grid(True)

            

            fig.suptitle(f'{self.motion_names[idx]}', size='x-large', weight='bold')
            plt.show()

        return t, y


def rho_alt(hp):
    return param.rho0 * (1 + param.lam * hp / param.Temp0) ** (-param.g / param.lam / param.R - 1)


def extract_state(data, tstamp, timestamps):
    row = data[data['time'] == tstamp]
    V0 = float(row['Dadc1_tas']) * ktstoms
    hp = float(row['Dadc1_bcAlt']) * ftm
    rho = rho_alt(hp)
    W = (param.mtot - float(row['lh_engine_FU']) * lbstokg - float(row['rh_engine_FU']) * lbstokg) * param.g
    theta = float(row['Ahrs1_Pitch']) * d2r
    CL = W / (0.5 * rho * V0 ** 2 * param.S)
    CX0 = W * np.sin(theta) / (0.5 * rho * V0 ** 2 * param.S)
    CZ0 = -W * np.cos(theta) / (0.5 * rho * V0 ** 2 * param.S)
    muc = (W / param.g) / (rho * param.S * param.c)
    mub = (W / param.g) / (rho * param.S * param.b)

    if (tstamp == timestamps['ph']) or (tstamp == timestamps['sp']):
        u = 0
        alpha = float(row['vane_AOA']) * d2r
        q = float(row['Ahrs1_bPitchRate']) * d2r * param.c / V0
        X0 = [u, alpha, theta, q]

    else:
        beta = 0
        phi = float(row['Ahrs1_Roll']) * d2r
        p = float(row['Ahrs1_bYawRate']) * d2r * param.b / (2*V0)
        r = float(row['Ahrs1_bRollRate']) * d2r * param.b / (2*V0)
        X0 = [beta, -phi, -r, -p]

    return {'V0': V0, 'CX0': CX0, 'CZ0': CZ0, 'CL': CL, 'muc': muc,
            'mub': mub, 'X0': X0}


def extract_input(data, tstamp, timestamps, dt):
    row = data[(data['time'] >= tstamp) & (data['time'] <= tstamp + dt)]
    if (tstamp == timestamps['ph']) or (tstamp == timestamps['sp']):
        de = row['delta_e'] * d2r
        de = de - de.iloc[0]
        return np.array(de)
    else:
        da = row['delta_a'] * d2r
        da = da - da.iloc[0]
        dr = row['delta_r'] * d2r
        dr = dr - dr.iloc[0]
        return np.array([np.array(da), np.array(dr)])


if __name__ == "__main__":
    datafile = '../postflight/flightdata.csv'
    excel_file = '../postflight/20200310_V2.xlsx'
    T = [120, 10, 15, 15, 15, 50]
    plane = Model(datafile=datafile, excel_file=excel_file)
    #plane.symmetric_nonzero_response()
    plane.asymmetric_nonzero_response()
    # plane.symmetric_nonzero_response(False)
    # plane.asymmetric_nonzero_response(False)
