import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
import myokit
import numpy as np
import os
import pandas as pd

import modelling


drug_list = modelling.SD_details.drug_names
# drug_list = ['dofetilide', 'bepridil']

results_dir = os.path.join(modelling.PARAM_DIR, 'Lei-SD-inference')
fig_dir = os.path.join(modelling.FIG_DIR, 'Lei_SD_fit')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Prepare control data
general_win = np.arange(1005 + 100, 1005 + 10000, 10)

Lei_control_state = myokit.load_state(
    os.path.join(results_dir, 'control_state.csv'))
Lei_control_log = myokit.DataLog.load_csv(
    os.path.join(results_dir, 'control_log.csv'))

# Set up Lei-SD model
Lei_sim = modelling.ModelSimController('Lei')
Lei_sim.set_ikr_rescale_method('AP_duration')
Lei_sim.initial_state = Lei_control_state

Li_sim = modelling.ModelSimController('Li')
Li_sim.set_conc(0)
Li_control_log = Li_sim.simulate(prepace=999,
                                 log_var=[Li_sim.time_key, Li_sim.ikr_key],
                                 log_times=general_win)
Li_sim.initial_state = Li_sim.sim.state()

dataset = modelling.DataLibrary()
plot = modelling.figures.FigurePlot()
ref_time = 945 + 100
skip_pt = 2

####################
# Plot fit of traces
####################
Li_err_dict = {}
for drug in drug_list:

    print('#######################')
    print(drug)

    # Load Li et al.'s experimental data
    dataset.set_drug(drug)
    conc_list = dataset.concs
    signal = dataset.get_mean_signal(cache=True)
    dataset.data = signal
    conc_num = len(conc_list)

    # Set up figure structure
    fig = modelling.figures.FigureStructure(figsize=(3, 1.5 * conc_num),
                                            gridspec=(conc_num, 1),
                                            height_ratios=[1] * conc_num,
                                            hspace=0.15)

    Li_sim.set_SD_parameters(drug)
    Lei_sim.set_SD_parameters(drug, ikr_model='Lei')

    Li_error = []
    Li_error = []
    for c, conc in enumerate(conc_list):
        # Load experimental traces of a specific drug concentration
        dataset.set_conc(conc)
        sweep_signal = dataset.conc_data.loc[(dataset.conc_data['sweep'] == 1),
                                             ['time', 'frac']]
        window = sweep_signal.loc[(sweep_signal['time'] >= ref_time),
                                  'time'].values

        log_times = []
        frac_block = []
        block_std = []
        for s in range(dataset.n_sweeps):
            # Create array of time points to log simulated IKr
            log_times.extend(window + 60 + 25e3 * s)

            # Combine data of all pulses together
            sweep_data = dataset.conc_data.loc[
                dataset.conc_data['sweep'] == s + 1,
                ['time', 'frac', 'frac_std']]
            # frac_block.extend([i * -1 + 1 for i in sweep_data.values])
            frac_block.extend(sweep_data.loc[sweep_data['time'].isin(window),
                                             'frac'].values)
            block_std.extend(sweep_data.loc[sweep_data['time'].isin(window),
                                            'frac_std'].values)
        Li_control_log = Li_control_log.trim(window[0] + 60,
                                             window[-1] + 60 + 10)
        Lei_control_log = Lei_control_log.trim(window[0] + 60,
                                               window[-1] + 60 + 10)

        # Plot mean and standard deviation of current signal
        plot_times = np.arange(0, len(frac_block)) * 10
#         frac_block = [i * -1 + 1 for i in frac_block]
        fig.axs[c][0].plot(plot_times[::skip_pt], frac_block[::skip_pt],
                           alpha=0.5, label='exp')
        fig.axs[c][0].fill_between(
            plot_times[::skip_pt],
            np.array(frac_block[::skip_pt]) - np.array(block_std[::skip_pt]),
            np.array(frac_block[::skip_pt]) + np.array(block_std[::skip_pt]),
            alpha=0.25)

        Li_sim.update_initial_state(paces=0)
        Li_sim.set_conc(conc)
        Li_log = Li_sim.simulate(prepace=0, save_signal=dataset.n_sweeps,
                                 log_times=log_times,
                                 log_var=[Li_sim.time_key, Li_sim.ikr_key],
                                 reset=False)

        # Simulate optimised IKr signals of the Lei-SD model
        Lei_sim.update_initial_state(paces=0)
        Lei_sim.set_conc(conc)
        Lei_log = Lei_sim.simulate(prepace=0, save_signal=dataset.n_sweeps,
                                   log_times=log_times,
                                   log_var=[Lei_sim.time_key, Lei_sim.ikr_key],
                                   reset=False)

        # Normalise IKr signals, to match with experimental data
        Li_plot_log = []
        Lei_plot_log = []
        for s in range(dataset.n_sweeps):
            Li_plot_log.extend(list(Li_log[Li_sim.ikr_key, s] /
                                    Li_control_log[Li_sim.ikr_key]))
            Lei_plot_log.extend(list(Lei_log[Lei_sim.ikr_key, s] /
                                     Lei_control_log[Lei_sim.ikr_key]))
        Li_error.append(np.sqrt(np.sum((
            np.array(Li_plot_log) - np.array(frac_block))**2) /
            len(frac_block)))

#         Li_plot_log = [i * -1 + 1 for i in Li_plot_log]
#         fig.axs[c][0].plot(plot_times[::skip_pt], Li_plot_log[::skip_pt],
#                            label='Li')

        # Plot simulated traces
#         Lei_plot_log = [i * -1 + 1 for i in Lei_plot_log]
        fig.axs[c][0].plot(plot_times[::skip_pt], Lei_plot_log[::skip_pt],
                           label='Lei')

        # Add drug concentration into figure as text
        conc_str = "{0:.0e}".format(conc)
        base, power = conc_str.split("e")
        power = int(power)
        fig.axs[c][0].text(0.05, 0.95,
                           base + r"$\times 10^{{{:d}}} $".format(power) + 'nM',
                           fontsize=8, ha='left', va='top',
                           transform=fig.axs[c][0].transAxes)

        fig.axs[c][0].set_ylim(0, 1.1)

    Li_error = np.mean(Li_error)
    Li_err_dict.update({drug: Li_error})

    # Adjust then save figure
    fig.axs[0][0].legend(ncols=2, fontsize='small', columnspacing=1,
                         handlelength=1.5)
    fig.axs[0][0].set_title(drug)
    fig.sharex(['Time (ms)'], [(0, len(frac_block) * 10)])
    fig.savefig(os.path.join(modelling.FIG_DIR, 'Lei_SD_fit',
                             f'{drug}_check.pdf'))

#########################
# Parameter comparison
#########################
# Set up figure structure
fig = modelling.figures.FigureStructure(figsize=(4 * 2, 1.5 * 3),
                                        gridspec=(3, 2),
                                        height_ratios=[1] * 3,
                                        hspace=0.15, wspace=0.25,
                                        plot_in_subgrid=True)
subgridspecs = [(1, 2)] * 6
fig.subgrid(subgridspecs, hspace=0.08, width_ratios=[3, 0.8])

# Get list of drugs and parameter names
drug_list = modelling.SD_details.drug_names
drug_list.sort()
param_list = modelling.SD_details.SD_param_names + ['error']

# Read parameter values for Li-SD model and Lei-SD model
Li_param = pd.read_csv(os.path.join(modelling.PARAM_DIR,
                                    'Li-SD.csv'), index_col=0)
Li_param = Li_param.loc[Li_param.index.isin(drug_list), :].sort_index()
Li_param['error'] = Li_err_dict

Lei_param = pd.read_csv(os.path.join(modelling.PARAM_DIR,
                                     'Lei-SD.csv'), index_col=0)
Lei_param = Lei_param.loc[Lei_param.index.isin(drug_list), :].sort_index()

for p, param in enumerate(param_list):
    print('#######################')
    print(param)

    # Plot drug parameter values for Li-SD model and Lei-SD model
    panel = fig.axs[p]
    panel[0][0].scatter(np.arange(len(drug_list)), Li_param[param],
                        color='none', lw=1.5, ec='orange', label='Li')
    panel[0][0].scatter(np.arange(len(drug_list)), Lei_param[param],
                        marker='^', color='none', lw=1.5, ec='k',
                        label='Lei')
    bp = panel[0][1].boxplot([Li_param[param], Lei_param[param]],
                             patch_artist=True,
                             medianprops=dict(color="white"))

    # Adjust figures
    colors = ['orange', 'grey']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    if param in ['Kmax', 'Ku', 'halfmax']:
        panel[0][0].set_yscale('log')
        panel[0][1].set_yscale('log')

    panel[0][0].set_ylabel(param)
    panel[0][0].grid(axis='x', linestyle='--')
    panel[0][0].set_axisbelow(True)
    for subPanel in range(2):
        panel[0][subPanel].spines[['right', 'top']].set_visible(False)
    panel[0][1].sharey(panel[0][0])
    panel[0][1].tick_params(labelleft=False)
    if p in [4, 5]:
        panel[0][0].set_xticks(np.arange(len(drug_list)), labels=drug_list)
        plt.setp(panel[0][0].get_xticklabels(), rotation=45,
                 ha='right', rotation_mode='anchor')
        panel[0][1].set_xticks(np.arange(1, 3), labels=['Li', 'Lei'])
        plt.setp(panel[0][1].get_xticklabels(), rotation=45,
                 ha='right', rotation_mode='anchor')
    else:
        panel[0][0].set_xticks(np.arange(len(drug_list)))
        panel[0][0].tick_params(labelbottom=False)
        panel[0][1].set_xticks(np.arange(1, 3))
        panel[0][1].tick_params(labelbottom=False)
fig.axs[1][0][0].legend(ncols=2, loc='lower right', columnspacing=1.2,
                        handletextpad=0.5, bbox_to_anchor=(1.0, 1.0))
fig.savefig(os.path.join(modelling.FIG_DIR, 'Lei_SD_fit',
                         'parameters_check.pdf'))

#########################
# Hil curve comparison
#########################
row, col = 3, 4
fig = modelling.figures.FigureStructure(figsize=(2 * col, 2 * row),
                                        gridspec=(row, col),
                                        hspace=0.3, wspace=0.15)

for d, drug in enumerate(drug_list):
    conc_list = modelling.SD_details.drug_concentrations[drug]['fine']

    Li_peaks = []
    Lei_peaks = []

    Li_sim.set_SD_parameters(drug)
    Lei_sim.set_SD_parameters(drug, ikr_model='Lei')
    for i in conc_list:
        Li_sim.set_conc(i)
        Li_log = Li_sim.simulate(log_var=[Li_sim.time_key, Li_sim.ikr_key])
        peak = Li_sim.extract_peak(Li_log)
        Li_peaks.append(peak[-1])

        Lei_sim.set_conc(i)
        Lei_log = Lei_sim.simulate(log_var=[Lei_sim.time_key, Lei_sim.ikr_key])
        peak = Lei_sim.extract_peak(Lei_log)
        Lei_peaks.append(peak[-1])

    Li_peaks_norm = (Li_peaks - min(Li_peaks)) / \
        (max(Li_peaks) - min(Li_peaks))
    Lei_peaks_norm = (Lei_peaks - min(Lei_peaks)) / \
        (max(Lei_peaks) - min(Lei_peaks))

    r, c = int(d / col), d % col
    fig.axs[r][c].plot(conc_list, Li_peaks_norm, label='Li')
    fig.axs[r][c].plot(conc_list, Lei_peaks_norm, linestyle='--', label='Lei')
    fig.axs[r][c].set_xscale('log')
    fig.axs[r][c].set_title(drug)

for c in range(col):
    fig.axs[2][c].set_xlabel('Drug concentration (log)')

fig.sharey(['Normalised peak current'] * 3)
fig.axs[0][0].legend()
fig.savefig(os.path.join(fig_dir, 'Hills.pdf'))
