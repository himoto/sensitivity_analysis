import numpy as np
from scipy.integrate import simps

from .model.name2idx import f_parameter as C
from .model.name2idx import f_variable as V
from .model import differential_equation as de
from .model.param_const import f_params
from .model.initial_condition import initial_values
from .ode_solver import solveode

# Signaling metric: time-integrated response
def sensitivity_analysis_integral():

    x = f_params()
    y0 = initial_values()

    condition = 2 # EGF & HRG
    tspan = range(5401) # -> 90 min.
    len_v = 57 # Num. of Rate equations (See differential_equation.py)

    rate = 1.01 # 1% change

    # Area under the curve
    integ_PcFos = np.empty((condition,len_v))

    for j in range(len_v):
        de.w = [1]*len_v
        de.w[j] = rate

        for i in range(condition):
            if i==0:
                y0[V.EGF] = 10.0
                y0[V.HRG] = 0.0
            elif i==1:
                y0[V.EGF] = 0.0
                y0[V.HRG] = 10.0

            (T,Y) = solveode(de.diffeq,y0,tspan,tuple(x))

            PcFos = Y[:,V.pcFOSn]*(x[C.Vn]/x[C.Vc]) + Y[:,V.pcFOSc]

            integ_PcFos[i,j] = simps(PcFos)

    # Sensitivity coefficient
    s_PcFos = np.log(integ_PcFos/integ_PcFos[:,0][:,None])/np.log(rate)

    return s_PcFos