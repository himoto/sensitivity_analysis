import sys
import numpy as np
from scipy.integrate import ode, simps

from .model.name2idx import parameters as C
from .model.name2idx import variables as V
from .model import differential_equation as ode_
from .model.param_const import f_params
from .model.initial_condition import initial_values
from .simulation import solveode, get_steady_state


def get_duration(temporal_dynamics):
    """
    Calculation of the duration as the time it takes to decline below 10% of its maximum.

    Parameters
    ----------
    temporal_dynamics: array
        Simulated time course data
    
    Returns
    -------
    duration: int

    """
    maximum_value = np.max(temporal_dynamics)
    t_max = np.argmax(temporal_dynamics)

    temporal_dynamics = temporal_dynamics - 0.1*maximum_value
    temporal_dynamics[temporal_dynamics > 0.0] = -np.inf

    duration = np.argmax(temporal_dynamics[t_max:]) + t_max

    return duration


def analyze_sensitivity(num_reaction):

    x = f_params()
    y0 = initial_values()

    conditions = ['EGF', 'HRG']
    tspan = range(5401)  # -> 90 min.

    rate = 1.01  # 1% change

    # Signaling metric
    duration_cFosmRNA = np.empty((len(conditions), num_reaction))
    integ_PcFos = np.empty((len(conditions), num_reaction))

    for j in range(num_reaction):
        ode_.perturbation = [1]*num_reaction
        ode_.perturbation[j] = rate
        # get steady state -- preprocess
        y0[V.EGF] = 0.0
        y0[V.HRG] = 0.0
        (T_steady_state, Y_steady_state) = get_steady_state(
            ode_.diffeq, y0, tspan, tuple(x)
        )
        if T_steady_state < tspan[-1]:
            print('Simulation failed.')
        else:
            y0 = Y_steady_state[:]
        # add ligand
        for i, condition in enumerate(conditions):
            if condition == 'EGF':
                y0[V.EGF] = 10.0
                y0[V.HRG] = 0.0
            elif condition == 'HRG':
                y0[V.EGF] = 0.0
                y0[V.HRG] = 10.0
            (T, Y) = solveode(ode_.diffeq, y0, tspan, tuple(x))

            cFosmRNA = Y[:, V.cfosmRNAc]
            PcFos = Y[:, V.pcFOSn]*(x[C.Vn]/x[C.Vc]) + Y[:, V.pcFOSc]

            duration_cFosmRNA[i, j] = get_duration(cFosmRNA)
            integ_PcFos[i, j] = simps(PcFos)

            sys.stdout.write('\r%d/%d' % (1+j, num_reaction))

    # Sensitivity coefficient
    s_cFosmRNA = np.log(duration_cFosmRNA /duration_cFosmRNA[:, 0][:, None]) / np.log(rate)
    s_PcFos = np.log(integ_PcFos / integ_PcFos[:, 0][:, None]) / np.log(rate)

    return s_cFosmRNA, s_PcFos
