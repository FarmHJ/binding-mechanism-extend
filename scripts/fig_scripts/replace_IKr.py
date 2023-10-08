# Introduces the idea of trapping and justifies the use of the Milnes protocol
import argparse
import matplotlib
import myokit.lib.plots as mp
import os

import modelling

fig_dir = os.path.join(modelling.FIG_DIR, 'basic_sim')

parser = argparse.ArgumentParser(
    description="Comparison between AP model and the AP-IKr model")
parser.add_argument('--mode', default='current_contribution',
                    choices=['all', 'current_contribution',
                             'reversal_potential', 'AP_hERG'],
                    help='Type of comparison to plot')
args = parser.parse_args()

plot_mode = args.mode

# Set up figure for reversal potential, AP and current contribution
plot = modelling.figures.FigurePlot()

if plot_mode in ['all', 'reversal_potential']:
    fig_EK = modelling.figures.FigureStructure(figsize=(3, 2), gridspec=(1, 1))

if plot_mode in ['all', 'current_contribution']:
    fig_current = modelling.figures.FigureStructure(figsize=(12, 5),
                                                    gridspec=(2, 4),
                                                    height_ratios=[1] * 2,
                                                    hspace=0.35, wspace=0.1)
    current_colours = modelling.model_naming.contribution_current_colours

if plot_mode in ['all', 'AP_hERG']:
    fig_AP = modelling.figures.FigureStructure(
        figsize=(12, 5), gridspec=(2, 4), height_ratios=[1] * 2,
        hspace=0.35, wspace=0.1, plot_in_subgrid=True)

    subgridspecs = [(2, 1)] * 8
    subgs = []
    for i in range(8):
        subgs.append(fig_AP.gs[i].subgridspec(*subgridspecs[i], wspace=0.1,
                                              hspace=0.1))
    axs_AP = [[[fig_AP.fig.add_subplot(subgs[k][i, j]) for j in range(
        subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
        k in range(len(subgs))]

# Figure parameters
cmap = matplotlib.cm.get_cmap('tab20')
plotting_pulse_time = 800

#############
# base model
#############
base_models = modelling.model_naming.APmodel_list[:4]
modified_models = modelling.model_naming.APmodel_list[1:]
modified_models = [modified_models[i] for i in [1, 2, 3, 0]]
reversal_label = ['rev.EK', 'parameters.ek', 'reversal_potentials.E_K',
                  'reversal_potentials.EK']

AP_y_lb, AP_y_ub = -50, 0
ikr_y_lb, ikr_y_ub = 0.5, 0.5

for n, model in enumerate(base_models):
    # Load ORd-Li model
    title = modelling.model_naming.AP_file_names[model]['label']
    if model == 'ORd-Li':
        APsim = modelling.ModelSimController(model)
    else:
        APsim = modelling.ModelSimController(model, ikr_modified=False)
    log = APsim.simulate(log_var='all')

    if plot_mode in ['all', 'reversal_potential']:
        print('Reversal potential of K for ', model, ' model: ',
              log[reversal_label[n]][0])

        # Plot reversal potential of potassium
        fig_EK.axs[0][0].plot(log.time(), log[reversal_label[n]], label=model)

    if plot_mode in ['all', 'AP_hERG']:
        # Plot AP and hERG
        SD_panel = axs_AP[n]
        plot.add_single(SD_panel[0][0], log, APsim.Vm_key)
        plot.add_single(SD_panel[1][0], log, APsim.ikr_key)
        SD_panel[0][0].set_title(title)
        AP_lb_panel, AP_ub_panel = SD_panel[0][0].get_ylim()
        ikr_lb_panel, ikr_ub_panel = SD_panel[1][0].get_ylim()

        fig_AP.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
                      axs=SD_panel, subgridspec=subgridspecs[n])

        AP_y_lb = min(AP_y_lb, AP_lb_panel)
        AP_y_ub = max(AP_y_ub, AP_ub_panel)
        ikr_y_lb = min(ikr_y_lb, ikr_lb_panel)
        ikr_y_ub = max(ikr_y_ub, ikr_ub_panel)

    if plot_mode in ['all', 'current_contribution']:
        # Plot current contribution
        APSD_current_keys = modelling.model_naming.model_current_keys[model]
        none_key_list = [i for i in APSD_current_keys.keys() if
                         APSD_current_keys[i] is None]
        none_key_list.extend(['time', 'Vm'])
        for i in none_key_list:
            del APSD_current_keys[i]

        colours = [cmap(current_colours[x]) for x in APSD_current_keys.keys()]
        currents = list(APSD_current_keys.values())

        SD_current_panel = fig_current.axs[0][n]
        SD_current_panel.set_title(title)
        mp.cumulative_current(log, currents, SD_current_panel, colors=colours,
                              normalise=True)
        SD_current_panel.set_rasterization_zorder(2)

#################
# modified model
#################

for n, model in enumerate(modified_models):
    # Load AP model
    title = modelling.model_naming.AP_file_names[model]['label']
    if model == 'ORd-Lei':
        title = model
    else:
        title = model + '-Li'
    APsim = modelling.ModelSimController(model)
    log = APsim.simulate(log_var='all')

    if plot_mode in ['all', 'reversal_potential']:
        print('Reversal potential of K for ', model, ' model: ',
              log[reversal_label[n]][0])

        # Plot reversal potential of potassium
        fig_EK.axs[0][0].plot(log.time(), log[reversal_label[n]], label=title)

    if plot_mode in ['all', 'AP_hERG']:
        # Plot AP and hERG
        SD_panel = axs_AP[n + 4]
        plot.add_single(SD_panel[0][0], log, APsim.Vm_key)
        plot.add_single(SD_panel[1][0], log, APsim.ikr_key)
        SD_panel[0][0].set_title(title)
        AP_lb_panel, AP_ub_panel = SD_panel[0][0].get_ylim()
        ikr_lb_panel, ikr_ub_panel = SD_panel[1][0].get_ylim()
        fig_AP.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
                      axs=SD_panel, subgridspec=subgridspecs[n])

        AP_y_lb = min(AP_y_lb, AP_lb_panel)
        AP_y_ub = max(AP_y_ub, AP_ub_panel)
        ikr_y_lb = min(ikr_y_lb, ikr_lb_panel)
        ikr_y_ub = max(ikr_y_ub, ikr_ub_panel)

    if plot_mode in ['all', 'current_contribution']:
        # Plot current contribution
        APSD_current_keys = modelling.model_naming.model_current_keys[model]
        none_key_list = [i for i in APSD_current_keys.keys() if
                         APSD_current_keys[i] is None]
        none_key_list.extend(['time', 'Vm'])
        for i in none_key_list:
            del APSD_current_keys[i]

        colours = [cmap(current_colours[x]) for x in APSD_current_keys.keys()]
        currents = list(APSD_current_keys.values())

        SD_current_panel = fig_current.axs[1][n]
        SD_current_panel.set_title(title)
        mp.cumulative_current(log, currents, SD_current_panel, colors=colours,
                              normalise=True)
        SD_current_panel.set_rasterization_zorder(2)

#################
# Adjust figures
#################
if plot_mode in ['all', 'reversal_potential']:
    fig_EK.axs[0][0].legend()
    fig_EK.savefig(os.path.join(fig_dir, 'reversal_potential.pdf'))

if plot_mode in ['all', 'AP_hERG']:
    for i in range(6):
        axs_AP[i][1][0].set_ylim(ikr_y_lb, ikr_y_ub)
        axs_AP[i][0][0].set_ylim(AP_y_lb, AP_y_ub)
        if i == 0 or i == 3:
            axs_AP[i][0][0].set_ylabel('AP')
            axs_AP[i][1][0].set_ylabel(r"$I_\mathrm{Kr}$")
        else:
            axs_AP[i][0][0].set_yticklabels([])
            axs_AP[i][1][0].set_yticklabels([])

    fig_AP.savefig(os.path.join(fig_dir, 'AP_hERG.pdf'))

if plot_mode in ['all', 'current_contribution']:
    legend_panel = fig_current.axs[1][2]
    legend_panel.xaxis.set_visible(False)
    legend_panel.yaxis.set_visible(False)
    legend_panel.set_frame_on(False)
    lines = []
    for current, i in current_colours.items():
        lines.append(matplotlib.lines.Line2D([0], [0], color=cmap(i), lw=5))
    del current_colours['IKACh']
    del current_colours['IKATP']
    labels = [modelling.model_naming.current_labels[x]
              for x in current_colours]
    legend_panel.legend(lines, labels, loc=(1.02, 0.05), ncol=2)

    fig_current.sharex(["Time (ms)"] * 3, [(0, plotting_pulse_time)] * 3)
    fig_current.sharey(["Relative contribution"] * 2, [(-1.02, 1.02)] * 2)

    fig_current.savefig(os.path.join(fig_dir, 'current_contribution.pdf'))
