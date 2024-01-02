# Introduces the idea of trapping and justifies the use of the Milnes protocol
import matplotlib
import myokit.lib.plots as mp
import os

import modelling

fig_dir = os.path.join(modelling.FIG_DIR, 'background')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Set up figure for reversal potential, AP and current contribution
plot = modelling.figures.FigurePlot()
fig = modelling.figures.FigureStructure(figsize=(7, 3.5),
                                        gridspec=(2, 2),
                                        height_ratios=[1] * 2,
                                        hspace=0.3, wspace=0.1)

# Figure parameters
model_details = modelling.model_naming
model_list = model_details.APmodel_list[:4]
current_list = model_details.current_list
current_colours = model_details.contribution_current_colours

cmap = matplotlib.colormaps['tab20']
plotting_pulse_time = 800

for num, APmodel_name in enumerate(model_list):
    model_current_keys = model_details.model_current_keys[APmodel_name]
    model_title = model_details.AP_file_names[APmodel_name]['label']

    APsim = modelling.ModelSimController(APmodel_name)
    if APmodel_name != 'ORd-Li':
        APsim.set_ikr_rescale_method('AP_duration')
    # Simulate AP
    log = APsim.simulate(log_var='all')

    none_key_list = [i for i in model_current_keys.keys() if
                     model_current_keys[i] is None]
    none_key_list.extend(['time', 'Vm'])
    for i in none_key_list:
        del model_current_keys[i]

    # Plot current contribution
    colours = [cmap(current_colours[x]) for x in model_current_keys.keys()]
    currents = list(model_current_keys.values())

    panel = fig.axs[int(num / 2)][num % 2]
    panel.set_title(model_title)
    mp.cumulative_current(log, currents, panel, colors=colours,
                          normalise=True)
    panel.set_rasterization_zorder(2)

legend_panel = fig.axs[1][1]
lines = []
for current, i in current_colours.items():
    lines.append(matplotlib.lines.Line2D([0], [0], color=cmap(i), lw=5))
del current_colours['IKACh']
del current_colours['IKATP']
labels = [model_details.current_labels[x] for x in current_colours]
legend_panel.legend(lines, labels, loc=(1.04, 0.45), ncol=1, handlelength=1)

fig.sharex(["Time (ms)"] * 2, [(0, plotting_pulse_time)] * 2)
fig.sharey(["Relative contribution"] * 2, [(-1.02, 1.02)] * 2)

# Save figure
fig.savefig(os.path.join(fig_dir, 'base_APmodels.svg'), format='svg')
