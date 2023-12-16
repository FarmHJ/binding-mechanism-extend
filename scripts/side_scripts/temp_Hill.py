import matplotlib.pyplot as plt
import os
import pandas as pd

import modelling


drug_list = modelling.SD_details.drug_names

results_dir = os.path.join(modelling.PARAM_DIR, 'Lei-SD-inference')
fig_dir = os.path.join(modelling.FIG_DIR, 'Lei_SD_fit')

Lei_sim = modelling.ModelSimController('Lei')
Lei_sim.update_initial_state()

Li_sim = modelling.ModelSimController('Li')
Li_sim.update_initial_state()

Farm_sim = modelling.ModelSimController('Lei')
Farm_sim.update_initial_state()
Farm_param = pd.read_csv(os.path.join(modelling.PARAM_DIR,
                                      'Farm-SD.csv'), index_col=0)
Farm_param = Farm_param.drop(columns=['error'])

row, col = 3, 4
fig = modelling.figures.FigureStructure(figsize=(2 * col, 2 * row),
                                        gridspec=(row, col),
                                        height_ratios=[1] * row,
                                        hspace=0.3, wspace=0.15)

for d, drug in enumerate(drug_list):
    print(drug)
    conc_list = modelling.SD_details.drug_concentrations[drug]['fine']
    Li_sim.set_SD_parameters(drug)
    Lei_sim.set_SD_parameters(drug, ikr_model='Lei')
    Farm_sim.set_parameters(Farm_param.loc[[drug], :])

    Li_peaks = []
    Lei_peaks = []
    Farm_peaks = []
    for c in conc_list:
        Li_sim.set_conc(c)
        log = Li_sim.simulate()
        peak = Li_sim.extract_peak(log)
        Li_peaks.append(peak[-1])

        Lei_sim.set_conc(c)
        log = Lei_sim.simulate()
        peak = Lei_sim.extract_peak(log)
        Lei_peaks.append(peak[-1])

        Farm_sim.set_conc(c)
        log = Farm_sim.simulate()
        peak = Farm_sim.extract_peak(log)
        Farm_peaks.append(peak[-1])

    Li_peaks_norm = (Li_peaks - min(Li_peaks)) / \
        (max(Li_peaks) - min(Li_peaks))
    Lei_peaks_norm = (Lei_peaks - min(Lei_peaks)) / \
        (max(Lei_peaks) - min(Lei_peaks))
    Farm_peaks_norm = (Farm_peaks - min(Farm_peaks)) / \
        (max(Farm_peaks) - min(Farm_peaks))

    r, c = int(d / col), d % col
    fig.axs[r][c].plot(conc_list, Li_peaks_norm, 'orange', label='Li')
    fig.axs[r][c].plot(conc_list, Lei_peaks_norm, 'g', label='Lei')
    fig.axs[r][c].plot(conc_list, Farm_peaks_norm, 'k--', label='Farm')
    fig.axs[r][c].set_xscale('log')
    fig.axs[r][c].set_title(drug)

for c in range(col):
    fig.axs[2][c].set_xlabel('Drug concentration (log)')

fig.sharey(['Normalised peak current'] * 3)
fig.axs[0][0].legend()
fig.fig.savefig(os.path.join(fig_dir, f'Hills.pdf'), bbox_inches='tight')
