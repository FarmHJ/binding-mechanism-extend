# Introduces the idea of trapping and justifies the use of the Milnes protocol
import myokit
import os

import modelling

# Define protocol
pulse_time = 1000
protocol_offset = 10
protocol = myokit.pacing.blocktrain(pulse_time, 1, offset=protocol_offset)

# Load AP-SD model
APmodel = '../math_model/AP_model/ohara-cipa-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel)
AP_model.protocol = protocol

repeats = 1000
abs_tol = 1e-7
rel_tol = 1e-8

# Simulate AP for AP clamp protocol
log = AP_model.model_simulation(repeats, abs_tol=abs_tol, rel_tol=rel_tol)
print('Reversal potential of K for AP-SD model: ', log['rev.EK'][0])
import matplotlib.pyplot as plt
plt.figure()
plt.plot(log.time(), log['potassium.K_i'])
fig_dir = '../figures/basic_sim/'
plt.savefig(fig_dir + 'potassium_current.pdf')
print(log['potassium.K_i'])

# Plot output
# fig_hERG = modelling.figures.FigureStructure(figsize=(5, 3), gridspec=(2, 3))
plot = modelling.figures.FigurePlot()

# plot.add_single(fig_hERG.axs[0][0], log, 'stimulus.i_stim')
# plot.add_single(fig_hERG.axs[1][0], log, 'membrane.V')
# plot.state_occupancy_plot(fig.axs[2][0], log, AP_model)

# fig.sharex(['Time (ms)'], [(0, pulse_time)])
# fig.sharey(['Current\n(A/F)', 'Voltage\n(mV)'])

# fig_dir = '../figures/basic_sim/'
# fig.savefig(fig_dir + "grd_sim.pdf")

# Load Grandi (2010) model
APmodel = '../math_model/AP_model/Grd-2010.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel)
AP_model.protocol = protocol

log_Grd = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)
print('Reversal potential of K for Grandi: ', log_Grd['parameters.ek'][0])

# # Plot output
# fig = modelling.figures.FigureStructure(figsize=(5, 3), gridspec=(2, 1))
# plot = modelling.figures.FigurePlot()

# plot.add_single(fig.axs[0][0], log, 'interface.I_app')
# plot.add_single(fig.axs[1][0], log, 'membrane_potential.V_m')
# # plot.state_occupancy_plot(fig.axs[2][0], log, AP_model)

# fig.sharex(['Time (ms)'], [(0, pulse_time)])
# fig.sharey(['Current\n(A/F)', 'Voltage\n(mV)'])

# fig_dir = '../figures/basic_sim/'
# fig.savefig(fig_dir + "grd_sim.pdf")

# Load Grandi-SD model
APmodel = '../math_model/AP_model/Grd-2010-IKr-SD.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel)
AP_model.protocol = protocol

log_Grd_SD = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)
print('Reversal potential of K for Grandi-SD: ', log_Grd_SD['parameters.ek'][0])

# Plot output
fig = modelling.figures.FigureStructure(figsize=(9, 6), gridspec=(8, 3),
                                        height_ratios=[1] * 8, hspace=0.2,
                                        wspace=0.08)

current_list_Grd = ['membrane_potential.V_m', 'I_Kr.I_kr', 'I_Na.I_Na',
                'I_Ca.I_Catot', 'I_Ks.I_ks', 'I_to.I_to', 'I_Ki.I_ki',
                'I_NaK.I_nak']
current_name = ['Vm', 'I_Kr', 'I_Na', 'I_CaL', 'I_Ks',
                'I_to', 'I_K1', 'I_NaK']
current_list_SD = ['membrane.V', 'ikr.IKr', 'ina.INa', 'ical.ICaL',
                   'iks.IKs', 'ito.Ito', 'ik1.IK1', 'inak.INaK']

for i in range(len(current_list_Grd)):
    plot.add_single(fig.axs[i][0], log, current_list_SD[i])
    b_lim1, t_lim1 = fig.axs[i][0].get_ylim()
    fig.axs[i][0].spines[['right', 'top']].set_visible(False)
    # fig.axs[i][0].set_ylabel(current_name[i])

    plot.add_single(fig.axs[i][1], log_Grd, current_list_Grd[i])
    b_lim2, t_lim2 = fig.axs[i][1].get_ylim()
    fig.axs[i][1].spines[['right', 'top']].set_visible(False)

    plot.add_single(fig.axs[i][2], log_Grd_SD, current_list_Grd[i])
    b_lim3, t_lim3 = fig.axs[i][2].get_ylim()
    fig.axs[i][2].spines[['right', 'top']].set_visible(False)

    fig.axs[i][0].set_ylim(bottom=min(b_lim1, b_lim2, b_lim3),
                           top=max(t_lim1, t_lim2, t_lim3))
    fig.axs[i][1].set_ylim(bottom=min(b_lim1, b_lim2, b_lim3),
                           top=max(t_lim1, t_lim2, t_lim3))
    fig.axs[i][2].set_ylim(bottom=min(b_lim1, b_lim2, b_lim3),
                           top=max(t_lim1, t_lim2, t_lim3))

fig.axs[0][0].set_title('AP-SD model')
fig.axs[0][1].set_title('Grandi model')
fig.axs[0][2].set_title('Grandi-SD model')
fig.sharex(['Time (ms)'] * 3, [(0, 500)] * 3)
fig.sharey(current_name)

fig_dir = '../figures/basic_sim/'
fig.savefig(fig_dir + "grd_SD_sim_currents.pdf")