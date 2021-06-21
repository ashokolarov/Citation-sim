import numpy as np
from scipy.optimize import minimize
from model import Model
from constants import *
from DataProcessing import DynamicData
import Cit_par21_opt as param
import control as ctr


def symmetric(x0, *args):
    Cma, CZa, Cmq, Cmde, Cmu, CXu, CZu, KY2, CXa, CXq, CZq, CXde, CZde, CZadot, Cmadot = x0
    V0, CX0, CZ0, muc = args

    C1 = np.array([[-2 * muc * (param.c / V0), 0, 0, 0],
                   [0, (CZadot - 2 * muc) * (param.c / V0), 0, 0],
                   [0, 0, -param.c / V0, 0],
                   [0, Cmadot * (param.c / V0), 0,
                    -2 * muc * KY2 * (param.c / V0)]], dtype=np.float64)

    C2 = np.array([[CXu, CXa, CZ0, CXq],
                   [CZu, CZa, -CX0, CZq + 2 * muc],
                   [0, 0, 0, 1],
                   [Cmu, Cma, 0, Cmq]], dtype=np.float64)

    C3 = np.array([CXde, CZde, 0, Cmde], dtype=np.float64).reshape(4, 1)

    A = -np.linalg.inv(C1) @ C2
    B = -np.linalg.inv(C1) @ C3
    C = np.eye(4)
    D = np.zeros(4).reshape(4, 1)

    return A, B, C, D


def cost_function_sym(x0_der):
    labels = ['sp', 'ph']
    T = [10, 120]
    metric = 0
    for idx, label in enumerate(labels):
        state = model.states[label]
        initial, V0, CX0, CZ0, muc, CL, mub = state['X0'], state['V0'], state['CX0'], state['CZ0'], state['muc'], state[
            'CL'], state['mub']

        A, B, C, D = symmetric(x0_der, V0, CX0, CZ0, muc)
        ss = ctr.ss(A, B, C, D)
        t = np.arange(0, T[idx] + model.dt, model.dt)

        X0 = np.array([0, 0, 0, initial[-1]]).reshape(4, 1)

        u = model.inputs[label]
        u = np.concatenate((u, np.ones(len(t) - len(u)) * u[0]), axis=None)
        _, y, x = ctr.forced_response(ss, t, u, X0)

        y[1] += initial[1]
        y[1] *= r2d
        y[2] += initial[2]
        y[2] *= r2d
        y[3] = y[3] * r2d * V0 / param.c

        fl = sym_mode_fl[idx]
        y_fl = fl[1:]

        w1 = np.array([5, 350, 55, 55])
        w1 = w1 / np.linalg.norm(w1)

        E = 0

        for i in range(4):
            E += w1[i] * np.sqrt(np.sum(np.divide((y[i]-y_fl[i]), y_fl[i]) ** 2)) / (len(y[i])-1)

        w2 = np.array([7,5])
        w2 = w2 / np.linalg.norm(w2)
        print(E)

        metric += w2[idx] * E

    return metric


def asymmetric(x0, *args):
    Clp, KX2, CYb, Cnb, Cnr, KZ2, Clb, Clr, CYdr, Clda, Cnda, Cndr, CYr, CYda, Cldr, Cnp, CXadot, CYp = x0
    V0, CL, mub = args
    C1 = np.array([[(param.CYbdot - 2 * mub) * param.b / V0, 0, 0, 0],
                   [0, -param.b / (2 * V0), 0, 0],
                   [0, 0, -4 * mub * KX2 * param.b / V0,
                    4 * mub * param.KXZ * param.b / V0],
                   [param.Cnbdot * param.b / V0, 0, 4 * mub * param.KXZ * param.b / V0,
                    -4 * mub * KZ2 * param.b / V0]], dtype=np.float64)

    C2 = np.array([[CYb, CL, param.CYp, param.CYr - 4 * mub],
                   [0, 0, 1, 0],
                   [Clb, 0, Clp, Clr],
                   [Cnb, 0, param.Cnp, Cnr]], dtype=np.float64)

    C3 = np.array([[CYda, CYdr],
                   [0, 0],
                   [Clda, Cldr],
                   [Cnda, Cndr]], dtype=np.float64)

    A = -np.linalg.inv(C1) @ C2
    B = -np.linalg.inv(C1) @ C3
    return A, B, np.eye(4), np.zeros((4, 2))


def cost_func_asym(x0):
    labels = ['dr', 'apr', 'spi']
    T = [15, 15, 50]
    metric = 0
    for idx, label in enumerate(labels):
        state = model.states[label]
        initial = state['X0']
        V0 = state['V0']
        CL = state['CL']
        mub = state['mub']
        A, B, C, D = asymmetric(x0, V0, CL, mub)
        ss = ctr.ss(A,B,C,D)

        t = np.arange(0, T[idx]+model.dt, model.dt)

        u = model.inputs[label]

        u1 = np.concatenate((u[0], np.ones(len(t) - len(u[0])) * u[0][0]), axis=None)
        u2 = np.concatenate((u[1], np.ones(len(t) - len(u[1])) * u[1][0]), axis=None)

        u = np.array([u1, u2])

        _, y, x = ctr.forced_response(ss, t, u, initial)

        y[1] *= -r2d
        y[2] = -y[2] * r2d * 2 * V0 / param.b
        y[3] = -y[3] * r2d * 2 * V0 / param.b

        if idx == 0:
            fl = asym_mode_fl[0]
        elif idx == 1:
            fl = asym_mode_fl[2]
        elif idx == 2:
            fl = asym_mode_fl[3]

        y_fl = fl[2:]

        E = 0
        w1 = np.array([2, 3, 3])
        w1 = w1 / np.linalg.norm(w1)

        for i in range(1,4):
            with np.errstate(divide='ignore', invalid='ignore'):
                a = (y[i] - y_fl[i])
                c = np.true_divide(a, y_fl[i])
                c[~ np.isfinite(c)] = 0  # -inf inf NaN
            E += w1[i-1] * np.sqrt(np.sum(c ** 2)) / (len(y[i]) - 1)
        print(E)

        metric += E

    return metric


if __name__ == "__main__":
    model = Model()
    flight = DynamicData()
    sym_mode_fl = flight.SymMode(1, 1, 0)
    asym_mode_fl = flight.AsymMode(1,1,1,1)

    x0 = np.array([param.Clp, param.KX2, param.CYb, param.Cnb, param.Cnr, param.KZ2, param.Clb, param.Clr, param.CYdr,
          param.Clda, param.Cnda, param.Cndr, param.CYr, param.CYda, param.Cldr, param.Cnp, param.CXadot,
                   param.CYp])

    cost_func_asym(x0)






