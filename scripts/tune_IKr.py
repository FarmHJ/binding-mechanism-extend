# Introduces the idea of trapping and justifies the use of the Milnes protocol
import matplotlib
import myokit
import myokit.lib.plots as mp
import numpy as np
import pandas as pd
from scipy.optimize import minimize

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
data_dir = '../simulation_data/'

# Set up figure for reversal potential, AP and current contribution
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

# Figure parameters
cmap = matplotlib.cm.get_cmap('tab20')
current_list = modelling.ModelDetails().current_list
current_keys = modelling.ModelDetails().current_keys
Grandi_current_keys = current_keys['Grandi']
plotting_pulse_time = 800
current_colours = modelling.ModelDetails().current_colours

#######################
#
# Grandi model
#
#######################
# Load Grandi (2010) model
APmodel = '../math_model/AP_model/Grd-2010.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel)
AP_model.protocol = protocol

log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)
Grandi_APD, _ = AP_model.APD90(log['membrane_potential.V_m'], protocol_offset,
                               0.1)
Grandi_flux = np.trapz(log['I_Kr.I_kr'], x=log.time())

# Plot AP and hERG
Grd_AP_panel = axs[0]
plot.add_single(Grd_AP_panel[0][0], log, 'membrane_potential.V_m')
plot.add_single(Grd_AP_panel[1][0], log, 'I_Kr.I_kr')
Grd_AP_panel[0][0].set_title("Grandi (2010)")
Grd_AP_panel[0][0].set_ylabel("AP")
Grd_AP_panel[1][0].set_ylabel(r"$I_\mathrm{Kr}$")
Grd_AP_panel[0][0].text(370, 0, 'APD90: ' + '{:.1f}'.format(Grandi_APD),
                        fontsize=8, ha='left')
AP_y_bottom1, AP_y_top1 = Grd_AP_panel[0][0].get_ylim()
IKr_y_bottom1, IKr_y_top1 = Grd_AP_panel[1][0].get_ylim()
Grd_AP_panel[1][0].text(450, (IKr_y_top1 - IKr_y_bottom1) / 2, 'scale: 1',
                        fontsize=8, ha='left')
fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
           axs=Grd_AP_panel, subgridspec=subgridspecs[0])

# Plot current contribution
none_key_list = [i for i in Grandi_current_keys.keys() if
                 Grandi_current_keys[i] is None]
for i in none_key_list:
    del(Grandi_current_keys[i])

colours = [cmap(current_colours[x]) for x in Grandi_current_keys.keys()]
currents = list(Grandi_current_keys.values())

Grd_current_panel = axs[4]
Grd_current_panel[0][0].set_xlabel("Time (ms)")
Grd_current_panel[0][0].set_xlim(0, plotting_pulse_time)
Grd_current_panel[0][0].set_ylabel("Relative contribution")
Grd_current_panel[0][0].set_ylim(-1.02, 1.02)
mp.cumulative_current(log, currents, Grd_current_panel[0][0], colors=colours,
                      normalise=True)

#######################
#
# Grandi-SD model
#
#######################
# Tuning method 1: scale conductance to have same peak of hERG current
# Load Grandi-SD model
APmodel = '../math_model/AP_model/Grd-2010-IKr-SD.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel)
AP_model.protocol = protocol

log1 = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)
SD_IKr_peak = max(log1[Grandi_current_keys['IKr']])
peak_IKr_scale = max(log[Grandi_current_keys['IKr']]) / SD_IKr_peak
log_tune1 = AP_model.model_simulation(1000,
                                      conductance_name='drug.ikr_rescale',
                                      conductance_value=peak_IKr_scale,
                                      abs_tol=abs_tol, rel_tol=rel_tol)
Grd_SD_APD, _ = AP_model.APD90(log_tune1['membrane_potential.V_m'],
                               protocol_offset, 0.1)

# Plot AP and hERG
Grd_SD_AP_panel = axs[1]
plot.add_single(Grd_SD_AP_panel[0][0], log_tune1, 'membrane_potential.V_m')
plot.add_single(Grd_SD_AP_panel[1][0], log_tune1, 'I_Kr.I_kr')
Grd_SD_AP_panel[0][0].text(370, 0, 'APD90: ' + '{:.1f}'.format(Grd_SD_APD),
                           fontsize=8, ha='left')
AP_y_bottom2, AP_y_top2 = Grd_SD_AP_panel[0][0].get_ylim()
IKr_y_bottom2, IKr_y_top2 = Grd_SD_AP_panel[1][0].get_ylim()
Grd_SD_AP_panel[1][0].text(450, (IKr_y_top2 - IKr_y_bottom2) / 2,
                           'scale: ' + '{:.2f}'.format(peak_IKr_scale),
                           fontsize=8, ha='left')
Grd_SD_AP_panel[0][0].set_title("Grandi-SD (1)")
fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
           axs=Grd_SD_AP_panel, subgridspec=subgridspecs[1])

# Plot current contribution
Grd_SD_current_panel = axs[5]
Grd_SD_current_panel[0][0].set_xlabel("Time (ms)")
Grd_SD_current_panel[0][0].set_xlim(0, plotting_pulse_time)
Grd_SD_current_panel[0][0].set_yticklabels([])
Grd_SD_current_panel[0][0].set_ylim(-1.02, 1.02)
mp.cumulative_current(log_tune1, currents, Grd_SD_current_panel[0][0],
                      colors=colours, normalise=True)

#######################
# Tuning method 2: scale conductance to have same APD90
#######################


def APD_problem(conductance_scale):
    log = AP_model.model_simulation(1000, conductance_name='drug.ikr_rescale',
                                    conductance_value=conductance_scale,
                                    abs_tol=1e-7, rel_tol=1e-8)
    APD, _ = AP_model.APD90(log['membrane_potential.V_m'],
                            protocol_offset, 0.1)

    error = np.sqrt(np.power(APD - Grandi_APD, 2))

    return error


initial_guess = peak_IKr_scale
res = minimize(APD_problem, initial_guess, method='nelder-mead',
               options={'disp': True})
APD_IKr_scale = res.x[0]
log_tune2 = AP_model.model_simulation(1000,
                                      conductance_name='drug.ikr_rescale',
                                      conductance_value=APD_IKr_scale,
                                      abs_tol=abs_tol, rel_tol=rel_tol)
Grd_SD_APD, _ = AP_model.APD90(log_tune2['membrane_potential.V_m'],
                               protocol_offset, 0.1)

# Plot AP and hERG
Grd_SD_AP_panel = axs[2]
plot.add_single(Grd_SD_AP_panel[0][0], log_tune2, 'membrane_potential.V_m')
plot.add_single(Grd_SD_AP_panel[1][0], log_tune2, 'I_Kr.I_kr')
Grd_SD_AP_panel[0][0].text(370, 0, 'APD90: ' + '{:.1f}'.format(Grd_SD_APD),
                           fontsize=8, ha='left')
AP_y_bottom3, AP_y_top3 = Grd_SD_AP_panel[0][0].get_ylim()
IKr_y_bottom3, IKr_y_top3 = Grd_SD_AP_panel[1][0].get_ylim()
Grd_SD_AP_panel[1][0].text(450, (IKr_y_top3 - IKr_y_bottom3) / 2,
                           'scale: ' + '{:.2f}'.format(APD_IKr_scale),
                           fontsize=8, ha='left')
Grd_SD_AP_panel[0][0].set_title("Grandi-SD (2)")
fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
           axs=Grd_SD_AP_panel, subgridspec=subgridspecs[2])

# Plot current contribution
Grd_SD_current_panel = axs[6]
Grd_SD_current_panel[0][0].set_xlabel("Time (ms)")
Grd_SD_current_panel[0][0].set_xlim(0, plotting_pulse_time)
Grd_SD_current_panel[0][0].set_yticklabels([])
Grd_SD_current_panel[0][0].set_ylim(-1.02, 1.02)
mp.cumulative_current(log_tune2, currents, Grd_SD_current_panel[0][0],
                      colors=colours, normalise=True)

#######################
# Tuning method 3: scale conductance to have same hERG current flux
#######################


def flux_problem(conductance_scale):
    log = AP_model.model_simulation(1000, conductance_name='drug.ikr_rescale',
                                    conductance_value=conductance_scale,
                                    abs_tol=1e-7, rel_tol=1e-8)
    flux = np.trapz(log['I_Kr.I_kr'], x=log.time())

    error = np.sqrt(np.power(flux - Grandi_flux, 2))

    return error


initial_guess = peak_IKr_scale
res = minimize(flux_problem, initial_guess, method='nelder-mead',
               options={'disp': True})
flux_IKr_scale = res.x[0]
log_tune3 = AP_model.model_simulation(1000,
                                      conductance_name='drug.ikr_rescale',
                                      conductance_value=flux_IKr_scale,
                                      abs_tol=abs_tol, rel_tol=rel_tol)
Grd_SD_APD, _ = AP_model.APD90(log_tune3['membrane_potential.V_m'],
                               protocol_offset, 0.1)

# Plot AP and hERG
Grd_SD_AP_panel = axs[3]
plot.add_single(Grd_SD_AP_panel[0][0], log_tune3, 'membrane_potential.V_m')
plot.add_single(Grd_SD_AP_panel[1][0], log_tune3, 'I_Kr.I_kr')
Grd_SD_AP_panel[0][0].text(370, 0, 'APD90: ' + '{:.1f}'.format(Grd_SD_APD),
                           fontsize=8, ha='left')
AP_y_bottom4, AP_y_top4 = Grd_SD_AP_panel[0][0].get_ylim()
IKr_y_bottom4, IKr_y_top4 = Grd_SD_AP_panel[1][0].get_ylim()
Grd_SD_AP_panel[1][0].text(450, (IKr_y_top4 - IKr_y_bottom4) / 2,
                           'scale: ' + '{:.2f}'.format(flux_IKr_scale),
                           fontsize=8, ha='left')
Grd_SD_AP_panel[0][0].set_title("Grandi-SD (3)")
fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
           axs=Grd_SD_AP_panel, subgridspec=subgridspecs[3])

# Plot current contribution
Grd_SD_current_panel = axs[7]
Grd_SD_current_panel[0][0].set_xlabel("Time (ms)")
Grd_SD_current_panel[0][0].set_xlim(0, plotting_pulse_time)
Grd_SD_current_panel[0][0].set_yticklabels([])
Grd_SD_current_panel[0][0].set_ylim(-1.02, 1.02)
mp.cumulative_current(log_tune3, currents, Grd_SD_current_panel[0][0],
                      colors=colours, normalise=True)

# Adjust figures
AP_y_min = min(AP_y_bottom1, AP_y_bottom2, AP_y_bottom3, AP_y_bottom4)
AP_y_max = max(AP_y_top1, AP_y_top2, AP_y_top3, AP_y_top4)
IKr_y_min = min(IKr_y_bottom1, IKr_y_bottom2, IKr_y_bottom3, IKr_y_bottom4)
IKr_y_max = max(IKr_y_top1, IKr_y_top2, IKr_y_top3, IKr_y_top4)
for i in range(4):
    axs[i][0][0].set_ylim(AP_y_min, AP_y_max)
    axs[i][1][0].set_ylim(IKr_y_min, IKr_y_max)
    if i == 0:
        axs[i][0][0].set_ylabel('AP')
        axs[i][1][0].set_ylabel(r"$I_\mathrm{Kr}$")
    else:
        axs[i][0][0].set_yticklabels([])
        axs[i][1][0].set_yticklabels([])

# Save figures
fig.savefig(fig_dir + 'IKr_tuning.pdf')

conductance_scale_df = pd.DataFrame(
    {'conductance scale': [peak_IKr_scale, APD_IKr_scale, flux_IKr_scale]},
    index=['hERG peak', 'AP duration', 'hERG flux'])
conductance_scale_df.to_csv(data_dir + 'Grandi_conductance_scale.csv')
