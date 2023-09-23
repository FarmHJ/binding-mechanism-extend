# Introduces the idea of trapping and justifies the use of the Milnes protocol
import argparse
import matplotlib
import myokit.lib.plots as mp
import numpy as np
import os
import pandas as pd
from scipy.optimize import minimize

import modelling

# Define AP model
parser = argparse.ArgumentParser(
    description="Tuning of IKr for AP-Li model")
parser.add_argument("APmodel", help="Name of AP model")
parser.add_argument('--method', default='all',
                    choices=['all', 'hERG_peak', 'hERG_flux', 'AP_duration'],
                    help='Method to tune the IKr ')
parser.add_argument('--noplot', action='store_false', help="Plot figures")
parser.add_argument('-x', '--cache', action='store_true', help='Use cache')
args = parser.parse_args()

plot_bool = args.noplot
tune_method = args.method
APmodel_name = args.APmodel

fig_dir = os.path.join(modelling.FIG_DIR, 'kinetics_comparison', APmodel_name)
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)
result_file = os.path.join(modelling.RESULT_DIR,
                           APmodel_name + '_conductance_scale.csv')

# Set up figure for reversal potential, AP and current contribution
if plot_bool:
    plot = modelling.figures.FigurePlot()
    fig = modelling.figures.FigureStructure(figsize=(9, 5), gridspec=(2, 4),
                                            height_ratios=[1] * 2, hspace=0.3,
                                            wspace=0.1, plot_in_subgrid=True)

    subgridspecs = [(2, 1)] * 4 + [(1, 1)] * 4
    subgs = []
    for i in range(8):
        subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], hspace=0.1))
    axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
        subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for k in
        range(len(subgs))]

    cmap = matplotlib.cm.get_cmap('tab20')

    # Figure parameters
    current_list = modelling.model_naming.current_list
    current_colours = modelling.model_naming.contribution_current_colours
    plotting_pulse_time = 800

model_keys = modelling.model_naming.model_current_keys[APmodel_name]
time_key = model_keys['time']
Vm_key = model_keys['Vm']
IKr_key = model_keys['IKr']

ikr_scale_dict = {}

#######################
#
# AP model
#
#######################
model_title = modelling.model_naming.AP_file_names[APmodel_name]['label']
APsim = modelling.ModelSimController(APmodel_name, ikr_modified=False)

log = APsim.simulate()
base_apd90 = APsim.APD90(log)
base_ikr_flux = np.trapz(log[IKr_key], x=log.time())
peak_IKr_scale = 1

# Plot AP and hERG
if plot_bool:
    AP_panel = axs[0]
    plot.add_single(AP_panel[0][0], log, Vm_key)
    plot.add_single(AP_panel[1][0], log, IKr_key)

    AP_panel[0][0].set_title(model_title)
    AP_panel[0][0].set_ylabel("AP")
    AP_panel[1][0].set_ylabel(r"$I_\mathrm{Kr}$")
    AP_panel[0][0].text(370, 0, 'APD90: ' + '{:.1f}'.format(base_apd90),
                        fontsize=8, ha='left')
    AP_y_bottom1, AP_y_top1 = AP_panel[0][0].get_ylim()
    IKr_y_bottom1, IKr_y_top1 = AP_panel[1][0].get_ylim()
    AP_panel[1][0].text(450, (IKr_y_top1 - IKr_y_bottom1) / 2, 'scale: 1',
                        fontsize=8, ha='left')
    fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
               axs=AP_panel, subgridspec=subgridspecs[0])

    # Plot current contribution
    none_key_list = [i for i in model_keys.keys() if model_keys[i] is None]
    for i in none_key_list:
        del(model_keys[i])

    plot_model_keys = [x for x in model_keys.keys() if x not in ['time', 'Vm']]
    colours = [cmap(current_colours[x]) for x in plot_model_keys]
    currents = [model_keys[x] for x in plot_model_keys]

    current_panel = axs[4]
    mp.cumulative_current(log, currents, current_panel[0][0], colors=colours,
                          normalise=True)
    current_panel[0][0].set_rasterization_zorder(2)

#######################
#
# AP-Li model
#
#######################
# Tuning method 1: scale conductance to have same peak of hERG current
# Load AP-Li model
APsim = modelling.ModelSimController(APmodel_name)

if tune_method in ['all', 'hERG_peak']:
    if args.cache:
        ikr_scale_df = pd.read_csv(result_file, index_col=[0],
                                   skipinitialspace=True)
        peak_IKr_scale = ikr_scale_df.loc['hERG_peak'].values[0]
    else:
        # Get the IKr scaling factor by matching the peak
        log1 = APsim.simulate()
        SD_IKr_peak = max(log1[IKr_key])
        del(log1)
        peak_IKr_scale = max(log[IKr_key]) / SD_IKr_peak
        print('hERG_peak: ', peak_IKr_scale)

        ikr_scale_dict.update({'hERG_peak': peak_IKr_scale})        

    APsim.set_ikr_rescale(peak_IKr_scale)
    log_tuned = APsim.simulate()
    apd90 = APsim.APD90(log_tuned)

    if plot_bool:
        # Plot AP and hERG
        AP_panel = axs[1]
        plot.add_single(AP_panel[0][0], log_tuned, Vm_key)
        plot.add_single(AP_panel[1][0], log_tuned, IKr_key)

        AP_panel[0][0].text(370, 0, 'APD90: ' + '{:.1f}'.format(apd90),
                            fontsize=8, ha='left')
        AP_y_bottom2, AP_y_top2 = AP_panel[0][0].get_ylim()
        IKr_y_bottom2, IKr_y_top2 = AP_panel[1][0].get_ylim()
        AP_panel[1][0].text(450, (IKr_y_top2 - IKr_y_bottom2) / 2,
                            'scale: ' + '{:.2f}'.format(peak_IKr_scale),
                            fontsize=8, ha='left')
        AP_panel[0][0].set_title(APmodel_name +
                                 r"-SD" "\n" r"$I_\mathrm{Kr}$ peak")
        fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
                   axs=AP_panel, subgridspec=subgridspecs[1])

        # Plot current contribution
        current_panel = axs[5]
        mp.cumulative_current(log_tuned, currents, current_panel[0][0],
                              colors=colours, normalise=True)
        current_panel[0][0].set_rasterization_zorder(2)

#######################
# Tuning method 2: scale conductance to have same APD90
#######################


def APD_problem(conductance_scale):
    APsim.set_ikr_rescale(conductance_scale)
    log = APsim.simulate(log_var=[time_key, Vm_key])

    apd90 = APsim.APD90(log)

    error = np.sqrt(np.power(apd90 - base_apd90, 2))

    return error


if tune_method in ['all', 'AP_duration']:
    if args.cache:
        ikr_scale_df = pd.read_csv(result_file, index_col=[0],
                                   skipinitialspace=True)
        conductance_scale = ikr_scale_df.loc['AP_duration'].values[0]
    else:
        # Get the IKr scaling factor by matching the APD
        initial_guess = peak_IKr_scale

        res = minimize(APD_problem, initial_guess, method='nelder-mead',
                       options={'disp': True})
        conductance_scale = res.x[0]
        print('AP_duration: ', conductance_scale)

        ikr_scale_dict.update({'AP_duration': conductance_scale})

    APsim.set_ikr_rescale(conductance_scale)
    log_tuned = APsim.simulate()
    apd90 = APsim.APD90(log_tuned)

    if plot_bool:
        # Plot AP and hERG
        AP_panel = axs[2]
        plot.add_single(AP_panel[0][0], log_tuned, Vm_key)
        plot.add_single(AP_panel[1][0], log_tuned, IKr_key)

        AP_panel[0][0].text(370, 0, 'APD90: ' + '{:.1f}'.format(apd90),
                            fontsize=8, ha='left')
        AP_y_bottom3, AP_y_top3 = AP_panel[0][0].get_ylim()
        IKr_y_bottom3, IKr_y_top3 = AP_panel[1][0].get_ylim()
        AP_panel[1][0].text(450, (IKr_y_top3 - IKr_y_bottom3) / 2,
                            'scale: ' + '{:.2f}'.format(conductance_scale),
                            fontsize=8, ha='left')
        AP_panel[0][0].set_title(APmodel_name + "-SD\nAPD")
        fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
                   axs=AP_panel, subgridspec=subgridspecs[2])

        # Plot current contribution
        current_panel = axs[6]
        mp.cumulative_current(log_tuned, currents, current_panel[0][0],
                              colors=colours, normalise=True)
        current_panel[0][0].set_rasterization_zorder(2)

#######################
# Tuning method 3: scale conductance to have same hERG current flux
#######################


def flux_problem(conductance_scale):
    APsim.set_ikr_rescale(conductance_scale)
    log = APsim.simulate(log_var=[time_key, IKr_key])

    flux = np.trapz(log[IKr_key], x=log.time())

    error = np.sqrt(np.power(flux - base_ikr_flux, 2))

    return error


if tune_method in ['all', 'hERG_flux']:
    if args.cache:
        ikr_scale_df = pd.read_csv(result_file, index_col=[0],
                                   skipinitialspace=True)
        conductance_scale = ikr_scale_df.loc['AP_duration'].values[0]
    else:
        # Get the IKr scaling factor by matching the flux of IKr
        initial_guess = peak_IKr_scale
        res = minimize(flux_problem, initial_guess, method='nelder-mead',
                       options={'disp': True})
        conductance_scale = res.x[0]
        print('hERG_flux: ', conductance_scale)

        ikr_scale_dict.update({'hERG_flux': conductance_scale})

    APsim.set_ikr_rescale(conductance_scale)
    log_tuned = APsim.simulate()
    apd90 = APsim.APD90(log_tuned)

    if plot_bool:
        # Plot AP and hERG
        AP_panel = axs[3]
        plot.add_single(AP_panel[0][0], log_tuned, Vm_key)
        plot.add_single(AP_panel[1][0], log_tuned, IKr_key)
        AP_panel[0][0].text(370, 0, 'APD90: ' + '{:.1f}'.format(apd90),
                            fontsize=8, ha='left')
        AP_y_bottom4, AP_y_top4 = AP_panel[0][0].get_ylim()
        IKr_y_bottom4, IKr_y_top4 = AP_panel[1][0].get_ylim()
        AP_panel[1][0].text(450, (IKr_y_top4 - IKr_y_bottom4) / 2,
                            'scale: ' + '{:.2f}'.format(conductance_scale),
                            fontsize=8, ha='left')
        AP_panel[0][0].set_title(APmodel_name +
                                 r"-SD" "\n" r"$I_\mathrm{Kr}$ flux")
        fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
                   axs=AP_panel, subgridspec=subgridspecs[3])

        # Plot current contribution
        current_panel = axs[7]
        mp.cumulative_current(log_tuned, currents, current_panel[0][0],
                              colors=colours, normalise=True)
        current_panel[0][0].set_rasterization_zorder(2)

#######################
# Tuning method 4: math APDs at different pacing rates
#######################

#######################
# Adjust figures
#######################
if plot_bool:
    AP_y_min = min(AP_y_bottom1, AP_y_bottom2, AP_y_bottom3, AP_y_bottom4)
    AP_y_max = max(AP_y_top1, AP_y_top2, AP_y_top3, AP_y_top4)
    IKr_y_min = min(IKr_y_bottom1, IKr_y_bottom2, IKr_y_bottom3, IKr_y_bottom4)
    IKr_y_max = max(IKr_y_top1, IKr_y_top2, IKr_y_top3, IKr_y_top4)
    for i in range(4):

        # Current panels
        axs[i + 4][0][0].set_xlabel("Time (ms)")
        axs[i + 4][0][0].set_xlim(0, plotting_pulse_time)
        axs[i + 4][0][0].set_yticklabels([])
        axs[i + 4][0][0].set_ylim(-1.02, 1.02)

        axs[i][0][0].set_ylim(AP_y_min, AP_y_max)
        axs[i][1][0].set_ylim(IKr_y_min, IKr_y_max)
        if i == 0:
            axs[i][0][0].set_ylabel('AP')
            axs[i][1][0].set_ylabel(r"$I_\mathrm{Kr}$")
        else:
            axs[i][0][0].set_yticklabels([])
            axs[i][1][0].set_yticklabels([])

    # Save figures
    # fig.savefig(os.path.join(fig_dir, 'IKr_tuning.pdf'))
    fig.savefig(os.path.join(modelling.FIG_DIR, 'checking_figures', 'IKr_tuning_test.pdf'))

if not args.cache:
    df = pd.DataFrame.from_dict(ikr_scale_dict, orient='index',
                                columns=['conductance scale'])
    df.to_csv(result_file)
