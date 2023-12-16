import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
import numpy as np
import os
import pandas as pd

import modelling

# Define AP model, drug and tuning method
parser = argparse.ArgumentParser(
    description="Fit parameters for Lei-SD model")
parser.add_argument('-d', '--drug', nargs='*', help='Choose drugs')
args = parser.parse_args()

if not args.drug:
    drug_list = modelling.SD_details.drug_names
else:
    drug_list = args.drug

results_dir = os.path.join(modelling.PARAM_DIR, 'Lei-SD-inference')
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)
fig_dir = os.path.join(modelling.FIG_DIR, 'Lei_SD_fit')

# Set up inference problem
sweep_num = 10
dataset = modelling.DataLibrary()

# Prepare control data
general_win = np.arange(1005, 1005 + 10000, 10)
Lei_sim = modelling.ModelSimController('Lei')
Lei_sim.set_conc(0)
Lei_control_log = Lei_sim.simulate(prepace=999,
                                   log_var=[Lei_sim.time_key, Lei_sim.ikr_key],
                                   log_times=general_win)
Lei_sim.initial_state = Lei_sim.sim.state()

Li_sim = modelling.ModelSimController('Li')
Li_sim.set_conc(0)
Li_control_log = Li_sim.simulate(prepace=999,
                                 log_var=[Li_sim.time_key, Li_sim.ikr_key],
                                 log_times=general_win)
Li_sim.initial_state = Li_sim.sim.state()

conc_num = 4
plot = modelling.figures.FigurePlot()
ref_time = 945 + 100
skip_pt = 2

###############################
# Fractional block comparison
###############################
Li_err_df = pd.DataFrame(np.zeros(len(drug_list)), index=drug_list, columns=['error'])
Lei_err_df = pd.DataFrame(np.zeros(len(drug_list)), index=drug_list, columns=['error'])

for drug in drug_list:

    print('#######################')
    print(drug)
    dataset.set_drug(drug)
    conc_list = dataset.concs

    conc_num = len(conc_list)
    fig = modelling.figures.FigureStructure(figsize=(3, 1.5 * conc_num),
                                            gridspec=(conc_num, 1),
                                            height_ratios=[1] * conc_num,
                                            hspace=0.15)

    signal = dataset.get_mean_signal(cache=True)
    dataset.data = signal

    Li_sim.set_SD_parameters(drug)
    Lei_sim.set_SD_parameters(drug, ikr_model='Lei')

    # scores_dict = pd.read_csv(os.path.join(results_dir, f'{drug}-Milnes.csv'),
    #                           index_col=0)
    # p = scores_dict.iloc[[0], :5]

    Li_error = []
    Lei_error = []
    for c, conc in enumerate(conc_list):
        # Prepare data
        dataset.set_conc(conc)

        # # Defining window to choose experimental data
        sweep_signal = dataset.conc_data.loc[(dataset.conc_data['sweep'] == 1),
                                             ['time', 'frac']]
        # ref_signal = sweep_signal.loc[(sweep_signal['time'] >= ref_time - 10)
        #                               & (sweep_signal['time'] <= ref_time + 10)]
        # if np.mean(ref_signal['frac'].values) < 0.5:
        #     window = sweep_signal.loc[(sweep_signal['frac'] < 0.5),
        #                               'time'].values
        # else:
        #     window = sweep_signal.loc[(sweep_signal['time'] >= ref_time),
        #                               'time'].values
        # window = sweep_signal['time'].values
        window = sweep_signal.loc[(sweep_signal['time'] >= ref_time),
                                  'time'].values

        log_times = []
        frac_block = []
        block_std = []
        for s in range(dataset.n_sweeps):
            log_times.extend(window + 60 + 25e3 * s)

            sweep_data = dataset.conc_data.loc[
                dataset.conc_data['sweep'] == s + 1,
                ['time', 'frac', 'frac_std']]
            # frac_block.extend([i * -1 + 1 for i in sweep_data.values])
            frac_block.extend(sweep_data.loc[sweep_data['time'].isin(window),
                                             'frac'].values)
#             block_std.extend(sweep_data.loc[sweep_data['time'].isin(window),
#                                             'frac_std'].values)
        Li_control_log = Li_control_log.trim(window[0] + 60,
                                             window[-1] + 60 + 10)
        Lei_control_log = Lei_control_log.trim(window[0] + 60,
                                               window[-1] + 60 + 10)

#         plot_times = np.arange(0, len(frac_block)) * 10
#         frac_block = [i * -1 + 1 for i in frac_block]
#         fig.axs[c][0].plot(plot_times[::skip_pt], frac_block[::skip_pt],
#                            alpha=0.5)
#         fig.axs[c][0].fill_between(
#             plot_times[::skip_pt],
#             np.array(frac_block[::skip_pt]) - np.array(block_std[::skip_pt]),
#             np.array(frac_block[::skip_pt]) + np.array(block_std[::skip_pt]),
#             alpha=0.25)

        # Li_model.set_conc(conc)
        Li_sim.update_initial_state(paces=0)
        Li_sim.set_conc(conc)
        Li_log = Li_sim.simulate(prepace=0, save_signal=sweep_num,
                                 log_times=log_times,
                                 log_var=[Li_sim.time_key, Li_sim.ikr_key],
                                 reset=False)

        Lei_sim.update_initial_state(paces=0)
        Lei_sim.set_conc(conc)
        Lei_log = Lei_sim.simulate(prepace=0, save_signal=sweep_num,
                                   log_times=log_times,
                                   log_var=[Lei_sim.time_key, Lei_sim.ikr_key],
                                   reset=False)

#         Lei_sim.set_parameters(p)
#         Lei_sim.update_initial_state(paces=0)
#         Farm_log = Lei_sim.simulate(prepace=0, save_signal=sweep_num,
#                                     log_times=log_times,
#                                     log_var=[Lei_sim.time_key, Lei_sim.ikr_key],
#                                     reset=False)
        Li_plot_log = []
        Lei_plot_log = []
        # Farm_plot_log = []
        for s in range(dataset.n_sweeps):
            Li_plot_log.extend(list(Li_log[Li_sim.ikr_key, s] /
                                    Li_control_log[Li_sim.ikr_key]))
            Lei_plot_log.extend(list(Lei_log[Lei_sim.ikr_key, s] /
                                     Lei_control_log[Lei_sim.ikr_key]))
            # Farm_plot_log.extend(list(Farm_log[Lei_sim.ikr_key, s] /
            #                           Lei_control_log[Lei_sim.ikr_key]))

        Li_error.append(np.sqrt(np.sum((np.array(Li_plot_log) - np.array(frac_block))**2) / len(frac_block)))
        Lei_error.append(np.sqrt(np.sum((np.array(Lei_plot_log) - np.array(frac_block))**2) / len(frac_block)))

#         Li_plot_log = [i * -1 + 1 for i in Li_plot_log]
#         fig.axs[c][0].plot(plot_times[::skip_pt], Li_plot_log[::skip_pt],
#                            label='Li')

#         Lei_plot_log = [i * -1 + 1 for i in Lei_plot_log]
#         fig.axs[c][0].plot(plot_times[::skip_pt], Lei_plot_log[::skip_pt],
#                            label='Lei')

#         Farm_plot_log = [i * -1 + 1 for i in Farm_plot_log]
#         fig.axs[c][0].plot(plot_times[::skip_pt], Farm_plot_log[::skip_pt],
#                            'k--', label='Farm')

#         conc_str = "{0:.0e}".format(conc)
#         base, power = conc_str.split("e")
#         power = int(power)
#         fig.axs[c][0].text(0.05, 0.95,
#                            base + r"$\times 10^{{{:d}}} $".format(power) + 'nM',
#                            fontsize=8, ha='left', va='top',
#                            transform=fig.axs[c][0].transAxes)

#         fig.axs[c][0].set_ylim(0, 1.1)

    Li_error = np.mean(Li_error)
    Lei_error = np.mean(Lei_error)
    Li_err_df['error'][drug] = Li_error
    Lei_err_df['error'][drug] = Lei_error

#     fig.axs[0][0].legend(ncols=3, fontsize='small', columnspacing=1,
#                          handlelength=1.5)
#     fig.axs[0][0].set_title(drug)
#     fig.sharex(['Time (ms)'], [(0, len(frac_block) * 10)])
#     fig.savefig(os.path.join(modelling.FIG_DIR, 'Lei_SD_fit',
#                              f'{drug}.pdf'))

#########################
# Parameter comparison
#########################

param_list = modelling.SD_details.SD_param_names
param_list.append('error')
fig = modelling.figures.FigureStructure(figsize=(4 * 2, 1.5 * 3),
                                        gridspec=(3, 2),
                                        height_ratios=[1] * 3,
                                        hspace=0.15, wspace=0.3,
                                        plot_in_subgrid=True)

subgridspecs = [(1, 2)] * 6
fig.subgrid(subgridspecs, hspace=0.08, width_ratios=[3, 1])

drug_list = modelling.SD_details.drug_names
Li_param = pd.read_csv(os.path.join(modelling.PARAM_DIR,
                                    'Li-SD.csv'), index_col=0)
Lei_param = pd.read_csv(os.path.join(modelling.PARAM_DIR,
                                     'Lei-SD.csv'), index_col=0)
Farm_param = pd.read_csv(os.path.join(modelling.PARAM_DIR,
                                      'Farm-SD.csv'), index_col=0)

Li_param = Li_param.loc[Li_param.index.isin(drug_list), :].sort_index()
Lei_param = Lei_param.loc[Lei_param.index.isin(drug_list), :].sort_index()
Farm_param = Farm_param.loc[Farm_param.index.isin(drug_list), :].sort_index()

for p, param in enumerate(param_list):

    print('#######################')
    print(param)

    # r, c = int(p / 2), p % 2
    panel = fig.axs[p]
    if param != 'error':
        panel[0][0].scatter(np.arange(len(drug_list)), Li_param[param],
                            color='none', lw=1.5, ec='orange', label='Li')
        panel[0][0].scatter(np.arange(len(drug_list)), Lei_param[param],
                            marker='^', color='none', lw=1.5, ec='g',
                            label='Lei')
        bp = panel[0][1].boxplot([Li_param[param], Lei_param[param],
                                  Farm_param[param]], patch_artist=True,
                                 medianprops=dict(color="white"))
    else:
        panel[0][0].scatter(np.arange(len(drug_list)), Li_err_df[param],
                            color='none', lw=1.5, ec='orange', label='Li')
        panel[0][0].scatter(np.arange(len(drug_list)), Lei_err_df[param],
                            marker='^', color='none', lw=1.5, ec='g',
                            label='Lei')
        bp = panel[0][1].boxplot([Li_err_df[param], Lei_err_df[param],
                                  Farm_param[param]], patch_artist=True,
                                 medianprops=dict(color="white"))
    panel[0][0].scatter(np.arange(len(drug_list)), Farm_param[param],
                        marker='s', color='none', lw=1.5, ec='k', label='Farm')

    colors = ['orange', 'g', 'grey']
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
        panel[0][1].set_xticks(np.arange(1, 4), labels=['Li', 'Lei', 'Farm'])
        plt.setp(panel[0][1].get_xticklabels(), rotation=45,
                 ha='right', rotation_mode='anchor')
    else:
        panel[0][0].set_xticks(np.arange(len(drug_list)))
        panel[0][0].tick_params(labelbottom=False)
        panel[0][1].set_xticks(np.arange(1, 4))
        panel[0][1].tick_params(labelbottom=False)
fig.axs[1][0][0].legend(ncols=3, loc='lower right', bbox_to_anchor=(1.0, 1.0))
fig.savefig(os.path.join(modelling.FIG_DIR, 'Lei_SD_fit',
                         'parameters.pdf'))


#########################
# Hil curve comparison
#########################

for drug in drug_list:
    conc_list = modelling.SD_details.drug_concentrations[drug]['fine']
    peaks = []

    Li_sim.set_SD_parameters(drug)
    for i in range(len(conc_list)):
        Li_sim.set_conc(conc_list[i])
        Li_sim.simulate()
        log = Li_sim.custom_simulation(
            param_values, drug_conc[i], 1000,
            log_var=['engine.time', 'ikr.IKr'],
            abs_tol=1e-7, rel_tol=1e-8)
        peak, _ = current_model.extract_peak(log, 'ikr.IKr')
        peaks.append(peak[-1])

    peaks_norm = (peaks - min(peaks)) / (max(peaks) - min(peaks))

    Hill_curve_SD = []
    for i in range(len(drug_conc)):
        reduction_scale = Hill_model.simulate(Hill_coef, drug_conc[i])
        Hill_curve_SD.append(reduction_scale)

    plt.figure()
    plt.plot(peaks_norm, label='Lei')
    plt.plot(Hill_curve_SD, label='SD')
    plt.legend()
    plt.savefig(fig_dir + 'LeivsSD_Hill.pdf')
