from .name2idx import variables as V

def initial_values():
    y0 = [0]*V.len_f_vars

    y0[V.RsD] = 2.47e+2
    y0[V.Kin] = 8.27e+1
    y0[V.MEK] = 6.37e+2
    y0[V.A1] = 1.82e+2
    y0[V.A2] = 2.54e+1
    y0[V.A3] = 1.31e+1
    #
    y0[V.ERKc] = 9.60e+02
    y0[V.RSKc] = 3.53e+02
    y0[V.CREBn] = 1.00e+03
    y0[V.Elk1n] = 1.51e+03

    return y0