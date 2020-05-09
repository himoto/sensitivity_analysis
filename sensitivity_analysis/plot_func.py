import numpy as np
import matplotlib.pyplot as plt


def timecourse(sim):
    plt.figure(figsize=(20, 8))
    plt.rcParams['font.size'] = 16
    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['xtick.major.width'] = 2
    plt.rcParams['ytick.major.width'] = 2
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['lines.linewidth'] = 2.5
    plt.rcParams['lines.markersize'] = 12
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'

    plt.subplots_adjust(wspace=0.5, hspace=0.5)

    plt.subplot(2, 4, 1)  # ----------------------------------------------------
    plt.plot(
        sim.t,sim.PMEK_cyt[:,0]/np.max(sim.PMEK_cyt),
        'b', label='EGF', clip_on=False
    )
    plt.plot(
        sim.t,sim.PMEK_cyt[:,1]/np.max(sim.PMEK_cyt),
        'r', label='HRG', clip_on=False
    )
    plt.plot(
        [var/60 for _, var in enumerate([0, 300, 600, 900, 1800, 2700, 3600, 5400])],
        [0.000, 0.773, 0.439, 0.252, 0.130, 0.087, 0.080, 0.066], 'bo'
    )
    plt.plot(
        [var/60 for _, var in enumerate([0, 300, 600, 900, 1800, 2700, 3600, 5400])],
        [0.000, 0.865, 1.000, 0.837, 0.884, 0.920, 0.875, 0.789], 'ro'
    )
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated MEK\n(cytoplasm)')
    plt.legend(loc='upper right',fontsize=12,frameon=False)

    plt.subplot(2, 4, 2)  # ----------------------------------------------------
    plt.plot(
        sim.t, sim.PERK_cyt[:, sim.conditions.index('EGF')]/np.max(sim.PERK_cyt), 
        'b', clip_on=False
    )
    plt.plot(
        sim.t, sim.PERK_cyt[:, sim.conditions.index('HRG')]/np.max(sim.PERK_cyt), 
        'r', clip_on=False
    )
    plt.xlim(0, 90)
    plt.xticks([0, 30, 60, 90])
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
    plt.ylim(0, 1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated ERK\n(cytoplasm)')

    plt.subplot(2, 4, 3)  # ----------------------------------------------------
    plt.plot(
        sim.t, sim.PRSK_wcl[:, sim.conditions.index('EGF')]/np.max(sim.PRSK_wcl), 
        'b', clip_on=False
    )
    plt.plot(
        sim.t, sim.PRSK_wcl[:, sim.conditions.index('HRG')]/np.max(sim.PRSK_wcl), 
        'r', clip_on=False
    )
    plt.xlim(0, 90)
    plt.xticks([0, 30, 60, 90])
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
    plt.ylim(0, 1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated RSK\n(whole cell)')

    plt.subplot(2, 4, 4)  # ----------------------------------------------------
    plt.plot(
        sim.t, sim.PCREB_wcl[:, sim.conditions.index('EGF')]/np.max(sim.PCREB_wcl), 
        'b', clip_on=False
    )
    plt.plot(
        sim.t, sim.PCREB_wcl[:, sim.conditions.index('HRG')]/np.max(sim.PCREB_wcl), 
        'r', clip_on=False
    )
    plt.xlim(0, 90)
    plt.xticks([0, 30, 60, 90])
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
    plt.ylim(0, 1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated CREB\n(whole cell)')

    plt.subplot(2, 4, 5)  # ----------------------------------------------------
    plt.plot(
        sim.t, sim.DUSPmRNA[:, sim.conditions.index('EGF')]/np.max(sim.DUSPmRNA), 
        'b', clip_on=False
    )
    plt.plot(
        sim.t, sim.DUSPmRNA[:, sim.conditions.index('HRG')]/np.max(sim.DUSPmRNA), 
        'r', clip_on=False
    )
    plt.xlim(0, 90)
    plt.xticks([0, 30, 60, 90])
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
    plt.ylim(0, 1.2)
    plt.xlabel('Time (min)')
    plt.ylabel(r'$\it{dusp}$'+' mRNA\nexpression')

    plt.subplot(2, 4, 6)  # ----------------------------------------------------
    plt.plot(
        sim.t, sim.cFosmRNA[:, sim.conditions.index('EGF')]/np.max(sim.cFosmRNA), 
        'b', clip_on=False
    )
    plt.plot(
        sim.t, sim.cFosmRNA[:, sim.conditions.index('HRG')]/np.max(sim.cFosmRNA), 
        'r', clip_on=False
    )
    plt.xlim(0, 90)
    plt.xticks([0, 30, 60, 90])
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
    plt.ylim(0, 1.2)
    plt.xlabel('Time (min)')
    plt.ylabel(r'$\it{c}$'+'-'+r'$\it{fos}$'+' mRNA\nexpression')

    plt.subplot(2, 4, 7)  # ----------------------------------------------------
    plt.plot(
        sim.t, sim.cFosPro[:, sim.conditions.index('EGF')]/np.max(sim.cFosPro), 
        'b', clip_on=False
    )
    plt.plot(
        sim.t, sim.cFosPro[:, sim.conditions.index('HRG')]/np.max(sim.cFosPro), 
        'r', clip_on=False
    )
    plt.xlim(0, 90)
    plt.xticks([0, 30, 60, 90])
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
    plt.ylim(0, 1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('c-Fos Protein\nexpression')

    plt.subplot(2, 4, 8)  # ----------------------------------------------------
    plt.plot(
        sim.t, sim.PcFos[:, sim.conditions.index('EGF')]/np.max(sim.PcFos), 
        'b', clip_on=False
    )
    plt.plot(
        sim.t, sim.PcFos[:, sim.conditions.index('HRG')]/np.max(sim.PcFos), 
        'r', clip_on=False
    )
    plt.xlim(0, 90)
    plt.xticks([0, 30, 60, 90])
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
    plt.ylim(0, 1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated c-Fos\nProtein expression')

    plt.show()