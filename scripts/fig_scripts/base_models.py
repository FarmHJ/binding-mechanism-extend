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

fig_dir = '../../figures/basic_sim/'

# Set up figure for reversal potential, AP and current contribution
plot = modelling.figures.FigurePlot()
fig = modelling.figures.FigureStructure(figsize=(10, 5),
                                        gridspec=(1, 2), wspace=0.1,
                                        plot_in_subgrid=True)

subgridspecs = [(2, 2), (1, 1)]
subgs = []

for i in range(2):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.1,
                                       hspace=0.3))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(1)]

# Figure parameters
cmap = matplotlib.cm.get_cmap('tab20')
current_list = modelling.ModelDetails().current_list
plotting_pulse_time = 800
current_colours = dict({
    'IKr': 0,
    'IKs': 1,
    'Ito': 2,
    'IKb': 3,
    'IK1': 4,
    'INaK': 5,
    'INa': 16,
    'INaL': 17,
    'ICaL': 10,
    'INaCa': 12,
    'INaB': 14,
    'ICaB': 15,
    'IClCa': 6,
    'IClB': 7,
    'ICaP': 13,
    'IKACh': 18,
    'IKATP': 19,
})

#######################
#
# AP-SD model
#
#######################
# Load AP-SD model
APmodel = '../../math_model/AP_model/ohara-cipa-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, current_head_key='ikr')
AP_model.protocol = protocol

# Simulate AP
log = AP_model.model_simulation(repeats, abs_tol=abs_tol, rel_tol=rel_tol)

# Plot current contribution
APSD_current_keys = modelling.ModelDetails().current_keys['AP-SD']
none_key_list = [i for i in APSD_current_keys.keys() if
                 APSD_current_keys[i] is None]
none_key_list.extend(['time', 'Vm'])
for i in none_key_list:
    del(APSD_current_keys[i])

colours = [cmap(current_colours[x]) for x in APSD_current_keys.keys()]
currents = list(APSD_current_keys.values())

SD_current_panel = axs[0][1][0]
SD_current_panel.set_title("O'Hara-CiPA (2017)")
mp.cumulative_current(log, currents, SD_current_panel, colors=colours,
                      normalise=True)
SD_current_panel.set_rasterization_zorder(2)

#######################
#
# Grandi model
#
#######################
# Load Grandi (2010) model
APmodel = '../../math_model/AP_model/Grd-2010.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, base_constant=None)
AP_model.protocol = protocol

log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)

# Plot current contribution
Grd_current_keys = modelling.ModelDetails().current_keys['Grandi']
none_key_list = [i for i in Grd_current_keys.keys() if
                 Grd_current_keys[i] is None]
none_key_list.extend(['time', 'Vm'])
for i in none_key_list:
    del(Grd_current_keys[i])

colours = [cmap(current_colours[x]) for x in Grd_current_keys.keys()]
currents = list(Grd_current_keys.values())

Grd_current_panel = axs[0][0][0]
Grd_current_panel.set_title("Grandi (2010)")
mp.cumulative_current(log, currents, Grd_current_panel, colors=colours,
                      normalise=True)
Grd_current_panel.set_rasterization_zorder(2)

#######################
#
# ten Tusscher (2006) model
#
#######################
# Load TTP (2006) model
APmodel = '../../math_model/AP_model/TTP-2006.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, base_constant=None)
AP_model.protocol = protocol

log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)

# Plot current contribution
TTP_current_keys = modelling.ModelDetails().current_keys['TTP']
none_key_list = [i for i in TTP_current_keys.keys() if
                 TTP_current_keys[i] is None]
none_key_list.extend(['time', 'Vm'])
for i in none_key_list:
    del(TTP_current_keys[i])

colours = [cmap(current_colours[x]) for x in TTP_current_keys.keys()]
currents = list(TTP_current_keys.values())

TTP_current_panel = axs[0][0][1]
TTP_current_panel.set_title("ten Tusscher (2006)")
mp.cumulative_current(log, currents, TTP_current_panel, colors=colours,
                      normalise=True)
TTP_current_panel.set_rasterization_zorder(2)


#######################
#
# Tomek (2019) model
#
#######################
# Load Tomek (2019) model
APmodel = '../../math_model/AP_model/Tomek-2019.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, base_constant=None)
AP_model.protocol = protocol

log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)

# Plot current contribution
Tomek_current_keys = modelling.ModelDetails().current_keys['Tomek']
none_key_list = [i for i in Tomek_current_keys.keys() if
                 Tomek_current_keys[i] is None]
none_key_list.extend(['time', 'Vm'])
for i in none_key_list:
    del(Tomek_current_keys[i])

colours = [cmap(current_colours[x]) for x in Tomek_current_keys.keys()]
currents = list(Tomek_current_keys.values())

Tomek_current_panel = axs[0][1][1]
Tomek_current_panel.set_title("Tomek (2019)")
mp.cumulative_current(log, currents, Tomek_current_panel, colors=colours,
                      normalise=True)
Tomek_current_panel.set_rasterization_zorder(2)

#######################
#
# Adjust figures
#
#######################

legend_panel = axs[0][1][1]
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
# legend_panel.legend(lines, labels, loc=(0.02, 0.05), ncol=2)

fig.sharex(["Time (ms)"] * 2, [(0, plotting_pulse_time)] * 2,
           axs=axs[0], subgridspec=(2, 2))
fig.sharey(["Relative contribution"] * 2, [(-1.02, 1.02)] * 2,
           axs=axs[0], subgridspec=(2, 2))

# Save figures
fig.savefig(fig_dir + 'base_models.pdf')
