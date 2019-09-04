import numpy as np
from matplotlib import pyplot as plt

from .reaction_module import set_reaction_module
from .duration import sensitivity_analysis_duration
from .integral import sensitivity_analysis_integral

def barplot():
    plt.figure(figsize=(12,10))
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['axes.linewidth'] = 1

    len_v = 57 # Num. of Rate Equations
    width = 0.3

    plt.subplot(2,1,1)
    s_cFosmRNA = sensitivity_analysis_duration()

    sort_idx, reaction_number = set_reaction_module(len_v,width)

    plt.bar(np.arange(len_v),s_cFosmRNA[0,sort_idx],
            width=width,color='b',align='center',label='EGF')
    plt.bar(np.arange(len_v)+width,s_cFosmRNA[1,sort_idx],
            width=width,color='r',align='center',label='HRG')

    for i in range(len_v-1):
        xp = range(len_v-1)[i]+width/2
        yp = s_cFosmRNA[np.argmax(np.abs(s_cFosmRNA[:,sort_idx[i]])),sort_idx[i]]
        if yp > 0:
            plt.text(xp,yp+0.05,reaction_number[i],
                     ha='center', va='bottom', fontsize=10, rotation=90)
        else:
            plt.text(xp,yp-0.05,reaction_number[i],
                     ha='center', va='top', fontsize=10, rotation=90)

    plt.xticks([])
    plt.ylabel('Control coefficients on\nduration ('+r'$\it{c}$'+'-'+r'$\it{fos}$'+' mRNA)')
    plt.xlim(-width,len_v-1-width)
    plt.ylim(-1.2,0.6)
    plt.yticks([-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6])
    plt.legend(loc='lower right',frameon=False)


    plt.subplot(2,1,2)
    s_PcFos = sensitivity_analysis_integral()

    sort_idx, reaction_number = set_reaction_module(len_v,width)

    plt.bar(np.arange(len_v),s_PcFos[0,sort_idx],
            width=width,color='b',align='center',label='EGF')
    plt.bar(np.arange(len_v)+width,s_PcFos[1,sort_idx],
            width=width,color='r',align='center',label='HRG')

    for i in range(len_v-1):
        xp = range(len_v-1)[i]+width/2
        yp = s_PcFos[np.argmax(np.abs(s_PcFos[:,sort_idx[i]])),sort_idx[i]]
        if yp > 0:
            plt.text(xp,yp+0.1,reaction_number[i],
                     ha='center', va='bottom', fontsize=10, rotation=90)
        else:
            plt.text(xp,yp-0.1,reaction_number[i],
                     ha='center', va='top', fontsize=10, rotation=90)

    plt.xticks([])
    plt.ylabel('Control coefficients on\nintegrated response (pc-Fos)')
    plt.xlim(-width,len_v-1-width)
    plt.ylim(-3,2)
    plt.legend(loc='lower right',frameon=False)

    plt.savefig('sensitivities.png',dpi=600,bbox_inches='tight')