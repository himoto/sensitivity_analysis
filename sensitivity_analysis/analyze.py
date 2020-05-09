import numpy as np
from matplotlib import pyplot as plt

from .sensitivity import calc_sensitivity_coefficients
from .reaction import *


def draw_vertical_span(reaction_module, num_reaction, width):
    left_end = 0
    for i, ith_module in enumerate(reaction_module):
        if i % 2 == 0:
            plt.axvspan(
                left_end - width,
                left_end + len(ith_module) - width,
                facecolor='k', alpha=0.1
            )
        left_end += len(ith_module)


def set_rc_params():
    plt.figure(figsize=(12, 5))
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['axes.linewidth'] = 1


def analyze():
    num_reaction = 57  # Num. of Rate Equations
    width = 0.3

    (s_cFosmRNA, s_PcFos) = calc_sensitivity_coefficients(num_reaction)
    reaction_module = get_reaction_module()

    sort_idx = [0]*num_reaction
    left_end = 0
    for i, ith_module in enumerate(reaction_module):
        for j, k in enumerate(ith_module):
            if i != 0 and j == 0:
                left_end += len(reaction_module[i-1])
            sort_idx[left_end+j] = k

    reaction_number = [str(i) for i in sort_idx]

    # --------------------------------------------------------------------------

    set_rc_params()

    draw_vertical_span(reaction_module, num_reaction, width)
    plt.hlines([0], -width, num_reaction-1-width, 'k', lw=1)
    plt.bar(
        np.arange(num_reaction), s_cFosmRNA[0, sort_idx], width=width,
        color='mediumblue', align='center', label='EGF'
    )
    plt.bar(
        np.arange(num_reaction)+width, s_cFosmRNA[1, sort_idx], width=width,
        color='red', align='center', label='HRG'
    )
    for i, j in enumerate(sort_idx):
        if j != 0:
            xp = i + width/2
            yp = s_cFosmRNA[np.argmax(np.abs(s_cFosmRNA[:, j])), j]
            if yp > 0:
                plt.text(
                    xp, yp+0.05, reaction_number[i], ha='center', va='bottom',
                    fontsize=10, rotation=90
                )
            else:
                plt.text(
                    xp, yp-0.05, reaction_number[i], ha='center', va='top',
                    fontsize=10, rotation=90
                )
    plt.xticks([])
    plt.ylabel(
        'Control coefficients on\nduration (' +
        r'$\it{c}$'+'-'+r'$\it{fos}$'+' mRNA)'
    )
    plt.xlim(-width, num_reaction-1-width)
    plt.ylim(-1.2, 0.6)
    plt.yticks([-1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6])
    plt.legend(loc='lower right', frameon=False)
    plt.savefig('sensitivity_cFosmRNA.png', dpi=300, bbox_inches='tight')
    plt.close()

    # --------------------------------------------------------------------------

    set_rc_params()

    draw_vertical_span(reaction_module, num_reaction, width)
    plt.hlines([0], -width, num_reaction-1-width, 'k', lw=1)
    plt.bar(
        np.arange(num_reaction), s_PcFos[0, sort_idx],
        width=width, color='mediumblue', align='center', label='EGF'
    )
    plt.bar(
        np.arange(num_reaction)+width, s_PcFos[1, sort_idx],
        width=width, color='red', align='center', label='HRG'
    )
    for i, j in enumerate(sort_idx):
        if j != 0:
            xp = i + width/2
            yp = s_PcFos[np.argmax(np.abs(s_PcFos[:, j])), j]
            if yp > 0:
                plt.text(
                    xp, yp+0.1, reaction_number[i], ha='center', va='bottom',
                    fontsize=10, rotation=90
                )
            else:
                plt.text(
                    xp, yp-0.1, reaction_number[i], ha='center', va='top',
                    fontsize=10, rotation=90
                )
    plt.xticks([])
    plt.ylabel('Control coefficients on\nintegrated response (pc-Fos)')
    plt.xlim(-width, num_reaction-1-width)
    plt.ylim(-3, 2)
    plt.legend(loc='lower right', frameon=False)

    plt.savefig('sensitivity_PcFos.png', dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    analyze()