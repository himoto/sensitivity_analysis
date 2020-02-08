import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects

def timecourse(sim):
    plt.figure(figsize=(20,8))
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

    plt.subplot(2,4,1)
    plt.plot(sim.t,sim.PMEK_cyt[:,0]/np.max(sim.PMEK_cyt[:,[0,1]]),'b',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],label='EGF'
    )
    plt.plot(sim.t,sim.PMEK_cyt[:,1]/np.max(sim.PMEK_cyt[:,[0,1]]),'r',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],label='HRG',
        clip_on=False
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

    plt.subplot(2,4,2)
    plt.plot(sim.t,sim.PERK_cyt[:,0]/np.max(sim.PERK_cyt[:,1]),'b',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.plot(sim.t,sim.PERK_cyt[:,1]/np.max(sim.PERK_cyt[:,1]),'r',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated ERK\n(cytoplasm)')

    plt.subplot(2,4,3)
    plt.plot(sim.t,sim.PRSK_wcl[:,0]/np.max(sim.PRSK_wcl[:,1]),'b',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.plot(sim.t,sim.PRSK_wcl[:,1]/np.max(sim.PRSK_wcl[:,1]),'r',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated RSK\n(whole cell)')

    plt.subplot(2,4,4)
    plt.plot(sim.t,sim.PCREB_wcl[:,0]/np.max(sim.PCREB_wcl[:,1]),'b',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.plot(sim.t,sim.PCREB_wcl[:,1]/np.max(sim.PCREB_wcl[:,1]),'r',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated CREB\n(whole cell)')

    plt.subplot(2,4,5)
    plt.plot(sim.t,sim.DUSPmRNA[:,0]/np.max(sim.DUSPmRNA[:,1]),'b',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.plot(sim.t,sim.DUSPmRNA[:,1]/np.max(sim.DUSPmRNA[:,1]),'r',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel(r'$\it{dusp}$'+' mRNA\nexpression')

    plt.subplot(2,4,6)
    plt.plot(sim.t,sim.cFosmRNA[:,0]/np.max(sim.cFosmRNA[:,1]),'b',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.plot(sim.t,sim.cFosmRNA[:,1]/np.max(sim.cFosmRNA[:,1]),'r',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel(r'$\it{c}$'+'-'+r'$\it{fos}$'+' mRNA\nexpression')

    plt.subplot(2,4,7)
    plt.plot(sim.t,sim.cFosPro[:,0]/np.max(sim.cFosPro[:,1]),'b',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.plot(sim.t,sim.cFosPro[:,1]/np.max(sim.cFosPro[:,1]),'r',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('c-Fos Protein\nexpression')

    plt.subplot(2,4,8)
    plt.plot(sim.t,sim.PcFos[:,0]/np.max(sim.PcFos[:,1]),'b',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.plot(sim.t,sim.PcFos[:,1]/np.max(sim.PcFos[:,1]),'r',
        path_effects=[path_effects.SimpleLineShadow(),
        path_effects.Normal()],clip_on=False
    )
    plt.xlim(0,90)
    plt.xticks([0,30,60,90])
    plt.yticks([0,0.2,0.4,0.6,0.8,1,1.2])
    plt.ylim(0,1.2)
    plt.xlabel('Time (min)')
    plt.ylabel('Phosphorylated c-Fos\nProtein expression')

    plt.show()