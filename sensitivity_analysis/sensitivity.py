import sys
import numpy as np
from scipy.integrate import ode, simps

from .model.name2idx import parameters as C
from .model.name2idx import variables as V
from .model import differential_equation as ode
from .model.param_const import f_params
from .model.initial_condition import initial_values


def solveode(diffeq, y0, tspan, args):
    sol = ode(diffeq)
    sol.set_integrator('vode', method='bdf', with_jacobian=True)
    sol.set_initial_value(y0, tspan[0])
    sol.set_f_params(args)
    T = [tspan[0]]
    Y = [y0]
    while sol.successful() and sol.t < tspan[-1]:
        sol.integrate(sol.t+1.)
        T.append(sol.t)
        Y.append(sol.y)

    return np.array(T), np.array(Y)


# Calculation of the duration as the time it takes to decline below 10% of its maximum
def get_duration(time_course_vector):
    maximum_value = np.max(time_course_vector)
    t_max = np.argmax(time_course_vector)

    time_course_vector = time_course_vector - 0.1*maximum_value
    time_course_vector[time_course_vector > 0.0] = -np.inf

    duration = np.argmax(time_course_vector[t_max:]) + t_max

    return duration


def analyze_sensitivity(num_reaction):

    x = f_params()
    y0 = initial_values()

    condition = 2  # EGF & HRG
    tspan = range(5401)  # -> 90 min.

    rate = 1.01  # 1% change

    # Signaling metric
    duration_cFosmRNA = np.empty((condition, num_reaction))
    integ_PcFos = np.empty((condition, num_reaction))

    for j in range(num_reaction):
        ode.perturbation = [1]*num_reaction
        ode.perturbation[j] = rate

        for i in range(condition):
            if i == 0:
                y0[V.EGF] = 10.0
                y0[V.HRG] = 0.0
            elif i == 1:
                y0[V.EGF] = 0.0
                y0[V.HRG] = 10.0
            (T, Y) = solveode(ode.diffeq, y0, tspan, tuple(x))

            cFosmRNA = Y[:, V.cfosmRNAc]
            PcFos = Y[:, V.pcFOSn]*(x[C.Vn]/x[C.Vc]) + Y[:, V.pcFOSc]

            duration_cFosmRNA[i, j] = get_duration(cFosmRNA)
            integ_PcFos[i, j] = simps(PcFos)

            sys.stdout.write('\r%d/%d' % (1+j, num_reaction))

    # Sensitivity coefficient
    s_cFosmRNA = np.log(duration_cFosmRNA /
                        duration_cFosmRNA[:, 0][:, None]) / np.log(rate)
    s_PcFos = np.log(integ_PcFos / integ_PcFos[:, 0][:, None]) / np.log(rate)

    return s_cFosmRNA, s_PcFos
