import numpy as np
import Cit_par21_post as param
from constants import *
from model import Model
import control as ctr


def short_period(muc):
    A = 4*muc**2*param.KY2
    B = -2*muc*(param.KY2*param.CZa+param.Cmadot+param.Cmq)
    C = param.CZa*param.Cmq - 2*muc*param.Cma

    roots = np.roots([A,B,C])
    return roots


def phugoid(CZ0, muc):
    A = 2*muc*(param.CZa*param.Cmq - 2*muc*param.Cma)
    B = 2*muc*(param.CXu*param.Cma - param.Cmu*param.CXa) + \
        param.Cmq*(param.CZu*param.CXa-param.CXu*param.CZa)
    C = CZ0*(param.CZa*param.Cmu-param.CZu*param.Cma)

    roots = np.roots([A,B,C])
    return roots


def aperiodic_roll(mub):
    result = param.Clp/(4*mub*param.KX2)

    return result

def dutch_roll(mub):
    A = 8*mub**2*param.KZ2
    B = -2*mub*(param.Cnr+2*param.KZ2*param.CYb)
    C = 4*mub*param.Cnb + param.CYb*param.Cnr

    roots = np.roots([A, B, C])
    return roots


def spiral(mub, CL):
    result = (2*CL*(param.Clb*param.Cnr-param.Cnb*param.Clr))\
             /(param.Clp*(4*mub*param.Cnb+param.CYb*param.Cnr)-
               param.Cnp*(4*mub*param.Clb+param.CYb*param.Clr))

    return result


if __name__ == "__main__":
    model = Model()
    mode = False

    if mode:
        state = model.states['ph']
        V0 = state['V0']
        CX0 = state['CX0']
        CZ0 = state['CZ0']
        muc = state['muc']
        CL = state['CL']
        mub = state['mub']

        A, B, C, D = model.symmetric(V0, CX0, CZ0, muc)
        print("Numerical model eigenvalues symmetric")
        ctr.damp(ctr.ss(A,B,C,D))
        eigs = np.linalg.eigvals(A)

        #Short period
        print('\nAnalytical eigenvalues from short period')
        eigs_anal = short_period(muc) * V0 / param.c
        print(eigs_anal)
        print("\nPercentage difference between num and analytical")
        diff = eigs[0] - eigs_anal[0]
        diff_real = abs(diff.real)/abs(eigs_anal[0].real) * 100
        diff_imag = abs(diff.imag) / abs(eigs_anal[0].imag) * 100
        print(f"{diff_real:.3f}%, {diff_imag:.3f}%")

        # Short period
        print('\nAnalytical eigenvalues from Phugoid')
        eigs_anal = phugoid(CZ0, muc) * V0 / param.c
        print(eigs_anal)
        print("\nPercentage difference between num and analytical")
        diff = eigs[2] - eigs_anal[0]
        diff_real = abs(diff.real) / abs(eigs_anal[0].real) * 100
        diff_imag = abs(diff.imag) / abs(eigs_anal[0].imag) * 100
        print(f"{diff_real:.3f}%, {diff_imag:.3f}%")

    else:
        state = model.states['dr']
        V0 = state['V0']
        CL = state['CL']
        mub = state['mub']
        A,B,C,D = model.asymmetric(V0, CL, mub)
        print("Numerical model eigenvalues asymmetric")
        ctr.damp(ctr.ss(A, B, C, D))
        eigs = np.linalg.eigvals(A)

        #Dutch roll
        print('\nAnalytical eigenvalues from dutch roll')
        eigs_anal = dutch_roll(mub) * V0 / param.b
        print(eigs_anal)
        print("\nPercentage difference between num and analytical")
        diff = eigs[1] - eigs_anal[0]
        diff_real = abs(diff.real) / abs(eigs_anal[0].real) * 100
        diff_imag = abs(diff.imag) / abs(eigs_anal[0].imag) * 100
        print(f"{diff_real:.3f}%, {diff_imag:.3f}%")

        # Aperiodic roll
        print('\nAnalytical eigenvalues from aperiodic roll')
        eigs_anal = aperiodic_roll(mub) * V0 / param.b
        print(eigs_anal)
        print("\nPercentage difference between num and analytical")
        diff = eigs[0] - eigs_anal
        diff_real = abs(diff) / abs(eigs_anal) * 100
        print(f"{diff_real:.3f}%")

        # Spiral
        print('\nAnalytical eigenvalues from spiral')
        eigs_anal = spiral(mub, CL) * V0 / param.b
        print(eigs_anal)
        print("\nPercentage difference between num and analytical")
        diff = eigs[3] - eigs_anal
        diff_real = abs(diff) / abs(eigs_anal) * 100
        print(f"{diff_real:.3f}%")

