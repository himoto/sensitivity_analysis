from matplotlib import pyplot as plt

def set_reaction_module(len_v,width):

    n_module = 15
    reaction_module = [None]*n_module

    # ERK_activation
    reaction_module[0] = [i for i in range(1,7)]

    # ERK_dephosphorylation_by_DUSP
    reaction_module[1] = [i for i in range(47,57)]

    # ERK_transport
    reaction_module[2] = [i for i in range(7,10)]

    # RSK_activation
    reaction_module[3] = [24,25]

    # RSK_transport
    reaction_module[4] = [26]

    # Elk1_activation
    reaction_module[5] = [29,30]

    # CREB_activation
    reaction_module[6] = [27,28]

    # dusp_production_etc
    reaction_module[7] = [i for i in range(10,14)]

    # DUSP_transport
    reaction_module[8] = [18,19]

    # DUSP_stabilization
    reaction_module[9] = [14,15,20,21]

    # DUSP_degradation
    reaction_module[10] = [16,17,22,23]

    # cfos_production_etc
    reaction_module[11] = [i for i in range(31,35)]

    # cFos_transport
    reaction_module[12] = [40,41]

    # cFos_stabilization
    reaction_module[13] = [35,36,37,42,43,44]

    # cFos_degradation
    reaction_module[14] = [38,39,45,46]

    # ==========================================================================

    sort_idx = [0]*len_v
    nth_module = 0
    for i in range(n_module):
        for j,k in enumerate(reaction_module[i]):
            if i!=0 and j==0:
                nth_module += len(reaction_module[i-1])
            sort_idx[nth_module+j] = k

    reaction_number = [str(i) for i in sort_idx]

    # Draw a vertical span
    nth_module = 0
    for i in range(n_module):
        if i%2 == 0:
            plt.axvspan(
                nth_module-width,
                nth_module+len(reaction_module[i])-width,
                facecolor='k',alpha=0.1
            )
        nth_module += len(reaction_module[i])
    ###

    return sort_idx, reaction_number