# Introduces the idea of trapping and justifies the use of the Milnes protocol
import matplotlib
import myokit
import myokit.lib.plots as mp

import modelling

# Define protocol
pulse_time = 1000
protocol_offset = 50
protocol = myokit.pacing.blocktrain(pulse_time, 0.5, offset=protocol_offset)

# Define constants
repeats = 1000
abs_tol = 1e-7
rel_tol = 1e-8

fig_dir = '../../figures/background/'

# Set up figure for reversal potential, AP and current contribution
plot = modelling.figures.FigurePlot()
fig = modelling.figures.FigureStructure(figsize=(10, 5),
                                        gridspec=(2, 2),
                                        height_ratios=[1] * 2,
                                        hspace=0.3, wspace=0.1)

# Figure parameters
model_list = ['ORd-CiPA', 'Grandi', 'TTP', 'Tomek-Cl']
# model_list = ['Tomek', 'Tomek2']
model_details = modelling.ModelDetails()
current_list = model_details.current_list
current_colours = model_details.current_colours

cmap = matplotlib.cm.get_cmap('tab20')
plotting_pulse_time = 800

for num, APmodel_name in enumerate(model_list):
    model_current_keys = model_details.current_keys[APmodel_name]
    none_key_list = [i for i in model_current_keys.keys() if
                     model_current_keys[i] is None]
    none_key_list.extend(['time', 'Vm'])
    for i in none_key_list:
        del(model_current_keys[i])

    # Load model
    model_filenames = model_details.file_names[APmodel_name]
    APmodel = '../../' + model_filenames['AP_path']
    model_title = model_filenames['label']

    APmodel, _, x = myokit.load(APmodel)
    AP_model = modelling.Simulation(APmodel, base_constant=None)
    AP_model.protocol = protocol

    # Simulate AP
    log = AP_model.model_simulation(repeats, abs_tol=abs_tol, rel_tol=rel_tol)

    # Plot current contribution
    colours = [cmap(current_colours[x]) for x in model_current_keys.keys()]
    currents = list(model_current_keys.values())

    panel = fig.axs[int(num / 2)][num % 2]
    panel.set_title(model_title)
    mp.cumulative_current(log, currents, panel, colors=colours,
                          normalise=True)
    panel.set_rasterization_zorder(2)

legend_panel = fig.axs[1][1]
# legend_panel = axs[0][1][1]
# legend_panel.xaxis.set_visible(False)
# legend_panel.yaxis.set_visible(False)
# legend_panel.set_frame_on(False)
lines = []
for current, i in current_colours.items():
    lines.append(matplotlib.lines.Line2D([0], [0], color=cmap(i), lw=5))
del(current_colours['IKACh'])
del(current_colours['IKATP'])
labels = [modelling.ModelDetails().current_names[x] for x in current_colours]
legend_panel.legend(lines, labels, loc=(1.04, 0.45), ncol=1, handlelength=1)

fig.sharex(["Time (ms)"] * 2, [(0, plotting_pulse_time)] * 2)
fig.sharey(["Relative contribution"] * 2, [(-1.02, 1.02)] * 2)

# Save figures
fig.savefig(fig_dir + 'base_APmodels.pdf')
