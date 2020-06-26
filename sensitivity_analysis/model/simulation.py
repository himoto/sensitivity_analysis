import numpy as np
from scipy.integrate import ode

from .name2idx import C, V
from .set_model import diffeq, param_values, initial_values


def solveode(diffeq, y0, tspan, args):
    sol = ode(diffeq)
    sol.set_integrator(
        'vode', method='bdf', with_jacobian=True, min_step=1e-8
    )
    sol.set_initial_value(y0, tspan[0])
    sol.set_f_params(args)

    T = [tspan[0]]
    Y = [y0]

    while sol.successful() and sol.t < tspan[-1]:
        sol.integrate(sol.t+1.)
        T.append(sol.t)
        Y.append(sol.y)

    return np.array(T), np.array(Y)


def get_steady_state(diffeq, y0, tspan, args, steady_state_time=7200):
    sol = ode(diffeq)
    sol.set_integrator(
        'vode', method='bdf', with_jacobian=True, min_step=1e-8
    )
    sol.set_initial_value(y0, tspan[0])
    sol.set_f_params(args)

    T = [tspan[0]]
    Y = [y0]

    while sol.successful() and sol.t < steady_state_time:
        sol.integrate(steady_state_time, step=True)
        
        T.append(sol.t)
        Y.append(sol.y)

    return T[-1], Y[-1]


class Simulation(object):

    tspan = range(5401)  # Unit time: 1 sec.
    conditions = ['EGF', 'HRG']

    t = np.array(tspan)/60.  # sec. -> min. (plot_func.py)

    PMEK_cyt = np.empty((len(tspan), len(conditions)))
    PERK_cyt = np.empty((len(tspan), len(conditions)))
    PRSK_wcl = np.empty((len(tspan), len(conditions)))
    PCREB_wcl = np.empty((len(tspan), len(conditions)))
    DUSPmRNA = np.empty((len(tspan), len(conditions)))
    cFosmRNA = np.empty((len(tspan), len(conditions)))
    cFosPro = np.empty((len(tspan), len(conditions)))
    PcFos = np.empty((len(tspan), len(conditions)))

    x = param_values()
    y0 = initial_values()

    # get steady state -- preprocess
    y0[V.EGF] = 0.0
    y0[V.HRG] = 0.0
    (T_steady_state, Y_steady_state) = get_steady_state(diffeq, y0, tspan, tuple(x))
    y0 = Y_steady_state[:]
    # add ligand
    for i, condition in enumerate(conditions):
        if condition == 'EGF':
            y0[V.EGF] = 10.0
            y0[V.HRG] = 0.0
        elif condition == 'HRG':
            y0[V.EGF] = 0.0
            y0[V.HRG] = 10.0

        (T, Y) = solveode(diffeq, y0, tspan, tuple(x))

        if T[-1] < tspan[-1]:
            print('Simulation failed.')
        else:
            PMEK_cyt[:, i] = Y[:, V.ppMEKc]
            PERK_cyt[:, i] = Y[:, V.pERKc] + Y[:, V.ppERKc]
            PRSK_wcl[:, i] = Y[:, V.pRSKc] + Y[:, V.pRSKn]*(x[C.Vn]/x[C.Vc])
            PCREB_wcl[:, i] = Y[:, V.pCREBn]*(x[C.Vn]/x[C.Vc])
            DUSPmRNA[:, i] = Y[:, V.duspmRNAc]
            cFosmRNA[:, i] = Y[:, V.cfosmRNAc]
            cFosPro[:, i] = (Y[:, V.pcFOSn] + Y[:, V.cFOSn]) * (x[C.Vn]/x[C.Vc]) \
                            + Y[:, V.cFOSc] + Y[:, V.pcFOSc]
            PcFos[:, i] = Y[:, V.pcFOSn]*(x[C.Vn]/x[C.Vc]) + Y[:, V.pcFOSc]