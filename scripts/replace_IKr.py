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

fig_dir = '../figures/basic_sim/'

# Set up figure for reversal potential, AP and current contribution
plot = modelling.figures.FigurePlot()
fig_EK = modelling.figures.FigureStructure(figsize=(3, 2), gridspec=(1, 1))
fig_current = modelling.figures.FigureStructure(figsize=(9, 5),
                                                gridspec=(2, 3),
                                                height_ratios=[1] * 2,
                                                hspace=0.35, wspace=0.1)
fig_AP = modelling.figures.FigureStructure(figsize=(9, 5), gridspec=(2, 3),
                                           height_ratios=[1] * 2, hspace=0.35,
                                           wspace=0.1, plot_in_subgrid=True)

subgridspecs = [(2, 1)] * 6
subgs = []
for i in range(6):
    subgs.append(fig_AP.gs[i].subgridspec(*subgridspecs[i], wspace=0.1,
                                          hspace=0.1))
axs_AP = [[[fig_AP.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(len(subgs))]

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
APmodel = '../math_model/AP_model/ohara-cipa-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, current_head_key='ikr')
AP_model.protocol = protocol

# Simulate AP
log = AP_model.model_simulation(repeats, abs_tol=abs_tol, rel_tol=rel_tol)
print('Reversal potential of K for AP-SD model: ', log['rev.EK'][0])

# Plot reversal potential of potassium
fig_EK.axs[0][0].plot(log.time(), log['rev.EK'], label='AP-SD')

# Plot AP and hERG
SD_panel = axs_AP[2]
plot.add_single(SD_panel[0][0], log, 'membrane.V')
plot.add_single(SD_panel[1][0], log, 'ikr.IKr')
SD_panel[0][0].set_title("O'Hara-CiPA (2017) epi")
AP_y_bottom1, AP_y_top1 = SD_panel[0][0].get_ylim()
IKr_y_bottom1, IKr_y_top1 = SD_panel[1][0].get_ylim()
fig_AP.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
              axs=SD_panel, subgridspec=subgridspecs[2])

# Plot current contribution
APSD_current_keys = modelling.ModelDetails().current_keys['AP-SD']
none_key_list = [i for i in APSD_current_keys.keys() if
                 APSD_current_keys[i] is None]
none_key_list.extend(['time', 'Vm'])
for i in none_key_list:
    del(APSD_current_keys[i])

colours = [cmap(current_colours[x]) for x in APSD_current_keys.keys()]
currents = list(APSD_current_keys.values())

SD_current_panel = fig_current.axs[0][2]
SD_current_panel.set_title("O'Hara-CiPA (2017) epi")
mp.cumulative_current(log, currents, SD_current_panel, colors=colours,
                      normalise=True)
SD_current_panel.set_rasterization_zorder(2)

#######################
#
# Grandi model
#
#######################
# Load Grandi (2010) model
APmodel = '../math_model/AP_model/Grd-2010.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, base_constant=None)
AP_model.protocol = protocol

log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)
print('Reversal potential of K for Grandi: ', log['parameters.ek'][0])

# Plot reversal potential of potassium
fig_EK.axs[0][0].plot(log.time(), log['parameters.ek'], label='Grandi')

# Plot AP and hERG
Grd_panel = axs_AP[0]
plot.add_single(Grd_panel[0][0], log, 'membrane_potential.V_m')
plot.add_single(Grd_panel[1][0], log, 'I_Kr.I_kr')
Grd_panel[0][0].set_title("Grandi (2010)")
AP_y_bottom2, AP_y_top2 = Grd_panel[0][0].get_ylim()
IKr_y_bottom2, IKr_y_top2 = Grd_panel[1][0].get_ylim()
fig_AP.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
              axs=Grd_panel, subgridspec=subgridspecs[0])

# Plot current contribution
Grd_current_keys = modelling.ModelDetails().current_keys['Grandi']
none_key_list = [i for i in Grd_current_keys.keys() if
                 Grd_current_keys[i] is None]
none_key_list.extend(['time', 'Vm'])
for i in none_key_list:
    del(Grd_current_keys[i])

colours = [cmap(current_colours[x]) for x in Grd_current_keys.keys()]
currents = list(Grd_current_keys.values())

Grd_current_panel = fig_current.axs[0][0]
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
APmodel = '../math_model/AP_model/TTP-2006.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, base_constant=None)
AP_model.protocol = protocol

log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)
print('Reversal potential of K for TTP: ', log['reversal_potentials.E_K'][0])

# Plot reversal potential of potassium
fig_EK.axs[0][0].plot(log.time(), log['reversal_potentials.E_K'], label='TTP')

# Plot AP and hERG
TTP_panel = axs_AP[1]
plot.add_single(TTP_panel[0][0], log, 'membrane.V')
plot.add_single(TTP_panel[1][0], log,
                'rapid_time_dependent_potassium_current.i_Kr')
TTP_panel[0][0].set_title("ten Tusscher (2006)")
AP_y_bottom3, AP_y_top3 = TTP_panel[0][0].get_ylim()
IKr_y_bottom3, IKr_y_top3 = TTP_panel[1][0].get_ylim()
fig_AP.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
              axs=TTP_panel, subgridspec=subgridspecs[1])

# Plot current contribution
TTP_current_keys = modelling.ModelDetails().current_keys['TTP']
none_key_list = [i for i in TTP_current_keys.keys() if
                 TTP_current_keys[i] is None]
none_key_list.extend(['time', 'Vm'])
for i in none_key_list:
    del(TTP_current_keys[i])

colours = [cmap(current_colours[x]) for x in TTP_current_keys.keys()]
currents = list(TTP_current_keys.values())

TTP_current_panel = fig_current.axs[0][1]
TTP_current_panel.set_title("ten Tusscher (2006)")
mp.cumulative_current(log, currents, TTP_current_panel, colors=colours,
                      normalise=True)
TTP_current_panel.set_rasterization_zorder(2)

#######################
#
# Grandi-SD model
#
#######################
# Load Grandi-SD model
APmodel = '../math_model/AP_model/Grd-2010-IKr-SD.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, current_head_key='I_Kr')
AP_model.protocol = protocol

log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)
print('Reversal potential of K for Grandi-SD: ', log['parameters.ek'][0])

# Plot reversal potential of potassium
fig_EK.axs[0][0].plot(log.time(), log['parameters.ek'], label='Grandi-SD')

# Plot AP and hERG
Grd_SD_panel = axs_AP[3]
plot.add_single(Grd_SD_panel[0][0], log, 'membrane_potential.V_m')
plot.add_single(Grd_SD_panel[1][0], log, 'I_Kr.I_kr')
Grd_SD_panel[0][0].set_title("Grandi-SD")
AP_y_bottom4, AP_y_top4 = Grd_SD_panel[0][0].get_ylim()
IKr_y_bottom4, IKr_y_top4 = Grd_SD_panel[1][0].get_ylim()
fig_AP.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
              axs=Grd_SD_panel, subgridspec=subgridspecs[3])

# Plot current contribution
Grd_current_keys = modelling.ModelDetails().current_keys['Grandi']
none_key_list = [i for i in Grd_current_keys.keys() if
                 Grd_current_keys[i] is None]
none_key_list.extend(['time', 'Vm'])
for i in none_key_list:
    del(Grd_current_keys[i])

colours = [cmap(current_colours[x]) for x in Grd_current_keys.keys()]
currents = list(Grd_current_keys.values())

Grd_SD_current_panel = fig_current.axs[1][0]
Grd_SD_current_panel.set_title("Grandi-SD")
mp.cumulative_current(log, currents, Grd_SD_current_panel,
                      colors=colours, normalise=True)
Grd_SD_current_panel.set_rasterization_zorder(2)

#######################
#
# ten Tusscher-SD model
#
#######################
# Load TTP-SD model
APmodel = '../math_model/AP_model/TTP-2006-IKr-SD.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(
    APmodel, current_head_key='rapid_time_dependent_potassium_current')
AP_model.protocol = protocol

log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)
print('Reversal potential of K for TTP: ', log['reversal_potentials.E_K'][0])

# Plot reversal potential of potassium
fig_EK.axs[0][0].plot(log.time(), log['reversal_potentials.E_K'],
                      label='TTP-SD')

# Plot AP and hERG
TTP_panel = axs_AP[4]
plot.add_single(TTP_panel[0][0], log, 'membrane.V')
plot.add_single(TTP_panel[1][0], log,
                'rapid_time_dependent_potassium_current.i_Kr')
TTP_panel[0][0].set_title("ten Tusscher-SD")
AP_y_bottom5, AP_y_top5 = TTP_panel[0][0].get_ylim()
IKr_y_bottom5, IKr_y_top5 = TTP_panel[1][0].get_ylim()
fig_AP.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
              axs=TTP_panel, subgridspec=subgridspecs[4])

# Plot current contribution
TTP_current_keys = modelling.ModelDetails().current_keys['TTP']
none_key_list = [i for i in TTP_current_keys.keys() if
                 TTP_current_keys[i] is None]
none_key_list.extend(['time', 'Vm'])
for i in none_key_list:
    del(TTP_current_keys[i])

colours = [cmap(current_colours[x]) for x in TTP_current_keys.keys()]
currents = list(TTP_current_keys.values())

TTP_SD_current_panel = fig_current.axs[1][1]
TTP_SD_current_panel.set_title("ten Tusscher-SD")
mp.cumulative_current(log, currents, TTP_SD_current_panel,
                      colors=colours, normalise=True)
TTP_SD_current_panel.set_rasterization_zorder(2)

#######################
#
# ORd-SD-Lei model
#
#######################
# Load ORd-SD-Lei model
APmodel = '../math_model/AP_model/ORd-CiPA-Lei-SD.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, current_head_key='ikr')
AP_model.protocol = protocol

log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)
print('Reversal potential of K for ORd-SD-Lei: ', log['rev.EK'][0])

# Plot reversal potential of potassium
fig_EK.axs[0][0].plot(log.time(), log['rev.EK'], label='ORd-SD-Lei')

# Plot AP and hERG
Lei_panel = axs_AP[5]
plot.add_single(Lei_panel[0][0], log, 'membrane.V')
plot.add_single(Lei_panel[1][0], log, 'ikr.IKr')
Lei_panel[0][0].set_title("ORd-SD-Lei")
AP_y_bottom6, AP_y_top6 = Lei_panel[0][0].get_ylim()
IKr_y_bottom6, IKr_y_top6 = Lei_panel[1][0].get_ylim()
fig_AP.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
              axs=Lei_panel, subgridspec=subgridspecs[5])

# Plot current contribution
Lei_current_keys = modelling.ModelDetails().current_keys['AP-SD']
none_key_list = [i for i in Lei_current_keys.keys() if
                 Lei_current_keys[i] is None]
none_key_list.extend(['time', 'Vm'])
for i in none_key_list:
    del(Lei_current_keys[i])

colours = [cmap(current_colours[x]) for x in Lei_current_keys.keys()]
currents = list(Lei_current_keys.values())

Lei_SD_current_panel = fig_current.axs[1][2]
Lei_SD_current_panel.set_title("ORd-SD-Lei")
mp.cumulative_current(log, currents, Lei_SD_current_panel,
                      colors=colours, normalise=True)
Lei_SD_current_panel.set_rasterization_zorder(2)

#######################
#
# Adjust figures
#
#######################
fig_EK.axs[0][0].legend()

AP_y_min = min(AP_y_bottom1, AP_y_bottom2, AP_y_bottom3, AP_y_bottom4,
               AP_y_bottom5, AP_y_bottom6)
AP_y_max = max(AP_y_top1, AP_y_top2, AP_y_top3, AP_y_top4, AP_y_top5,
               AP_y_top6)
IKr_y_min = min(IKr_y_bottom1, IKr_y_bottom2, IKr_y_bottom3, IKr_y_bottom4,
                IKr_y_bottom5, IKr_y_bottom6)
IKr_y_max = max(IKr_y_top1, IKr_y_top2, IKr_y_top3, IKr_y_top4, IKr_y_top5,
                IKr_y_top6)
for i in range(6):
    axs_AP[i][1][0].set_ylim(IKr_y_min, IKr_y_max)
    axs_AP[i][0][0].set_ylim(AP_y_min, AP_y_max)
    if i == 0 or i == 3:
        axs_AP[i][0][0].set_ylabel('AP')
        axs_AP[i][1][0].set_ylabel(r"$I_\mathrm{Kr}$")
    else:
        axs_AP[i][0][0].set_yticklabels([])
        axs_AP[i][1][0].set_yticklabels([])

legend_panel = fig_current.axs[1][2]
# legend_panel.xaxis.set_visible(False)
# legend_panel.yaxis.set_visible(False)
# legend_panel.set_frame_on(False)
lines = []
for current, i in current_colours.items():
    lines.append(matplotlib.lines.Line2D([0], [0], color=cmap(i), lw=5))
del(current_colours['IKACh'])
del(current_colours['IKATP'])
labels = [modelling.ModelDetails().current_names[x] for x in current_colours]
legend_panel.legend(lines, labels, loc=(1.02, 0.05), ncol=2)

fig_current.sharex(["Time (ms)"] * 3, [(0, plotting_pulse_time)] * 3)
fig_current.sharey(["Relative contribution"] * 2, [(-1.02, 1.02)] * 2)

# Save figures
fig_EK.savefig(fig_dir + 'reversal_potential.pdf')
fig_AP.savefig(fig_dir + 'AP_hERG.pdf')
fig_current.savefig(fig_dir + 'current_contribution.pdf')
