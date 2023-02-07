import matplotlib.pyplot as plt
import myokit
import numpy as np

import modelling

# Define protocol
pulse_time = 1000
protocol_offset = 10
protocol = myokit.pacing.blocktrain(pulse_time, 0.5, offset=protocol_offset)

#
# Grandi (2010) model
#
# # Set up AP model
# APmodel = '../math_model/AP_model/Grd-2010.mmt'
# APmodel, _, x = myokit.load(APmodel)
# AP_model = modelling.Simulation(APmodel)
# AP_model.protocol = protocol

abs_tol = 1e-7
rel_tol = 1e-8
# log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)

# AP_peak, _ = AP_model.extract_peak(log, 'membrane_potential.V_m')
# APD90 = AP_model.APD90(log['membrane_potential.V_m'], protocol_offset, timestep=0.1)
# resting_Vm = min(log['membrane_potential.V_m'][500 / 0.1:])
# print('resting membrane potential: ', resting_Vm)
# print('AP amplitude:', AP_peak[0] - resting_Vm)
# print('APD90: ', APD90)
# print('maximal upstroke velocity: ', 
#       max(log['dot(membrane_potential.V_m)']))

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

# # Plot output
# fig = modelling.figures.FigureStructure(figsize=(8, 5), gridspec=(6, 2),
#                                         height_ratios=[1] * 6, hspace=0.2)

# current_list = ['membrane_potential.V_m', 'Ca_Concentrations.Ca_i',
#                 'Ca_Concentrations.Ca_j', 'Ca_Concentrations.Ca_sl',
#                 'I_Na.I_Na', 'I_Ca.I_Catot', 'I_Kr.I_kr', 'I_Ks.I_ks',
#                 'I_to.I_to', 'I_Ki.I_ki', 'I_NCX.I_ncx', 'I_NaK.I_nak']
# current_name = ['Vm', '[Ca]i', '[Ca]jct', '[Ca]SL',
#                 'I_Na', 'I_CaL', 'I_Kr', 'I_Ks',
#                 'I_to', 'I_K1', 'I_NCX', 'I_NaK']

# for i in range(len(current_list)):
#     row = i % 6
#     col = i // 6
#     plot.add_single(fig.axs[row][col], log, current_list[i])
#     fig.axs[row][col].spines[['right', 'top']].set_visible(False)
#     fig.axs[row][col].set_ylabel(current_name[i])
# # plot.state_occupancy_plot(fig.axs[2][0], log, AP_model)

# fig.sharex(['Time (ms)'] * 2, [(0, 500)] * 2)
# # fig.sharey(['Current\n(A/F)', 'Voltage\n(mV)'])

# fig_dir = '../figures/basic_sim/'
# fig.savefig(fig_dir + "grd_sim_currents.pdf")

#
# ten Tusscher (2006) model
#
# Define protocol
timestep = 0.1
protocol_offset = 100
protocol = myokit.pacing.blocktrain(pulse_time, 0.5, offset=protocol_offset)
# Set up AP model
APmodel = '../math_model/AP_model/TTP-2006.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel)
AP_model.protocol = protocol

log = AP_model.model_simulation(1000, timestep=timestep, abs_tol=abs_tol, rel_tol=rel_tol)

# Action potential measures
AP_peak, _ = AP_model.extract_peak(log, 'membrane.V')
APD90 = AP_model.APD90(log['membrane.V'], protocol_offset, timestep=timestep)
resting_Vm = min(log['membrane.V'][500 * 10:])
print('###############')
print('resting membrane potential: ', resting_Vm)
print('AP amplitude:', AP_peak[0] - resting_Vm)
print('APD90: ', APD90)
print('maximal upstroke velocity: ', 
      max(log['dot(membrane.V)']))

# Plot TTP's action potential
fig = modelling.figures.FigureStructure(figsize=(8, 3), gridspec=(2, 1))
plot = modelling.figures.FigurePlot()

plot.add_single(fig.axs[0][0], log, 'membrane.i_Stim')
plot.add_single(fig.axs[1][0], log, 'membrane.V')
# plot.state_occupancy_plot(fig.axs[2][0], log, AP_model)

fig.sharex(['Time (ms)'], [(0, pulse_time)])
fig.sharey(['Current\n(A/F)', 'Voltage\n(mV)'])

fig_dir = '../figures/basic_sim/'
fig.savefig(fig_dir + "ttp_sim.pdf")

# Calcium dynamics measures
Ca_SS = log['calcium_dynamics.Ca_ss']
CaSS_peak_amplitude = max(Ca_SS) - min(Ca_SS)
print('CaSS peak amplitude: ', CaSS_peak_amplitude)
Ca_i = log['calcium_dynamics.Ca_i']
Cai_peak_amplitude = max(Ca_i) - min(Ca_i)
print('Cai peak amplitude: ', Cai_peak_amplitude)
# CaSS_time2peak = (np.argmax(Ca_SS) - np.argmax(log['dot(calcium_dynamics.Ca_ss)'])) * timestep
CaSS_time2peak = np.sum(log['dot(calcium_dynamics.Ca_ss)'] > 0.01 * 1) * timestep
print('CaSS time to peak: ', CaSS_time2peak)
# Cai_time2peak = (np.argmax(Ca_i) - np.argmax(log['dot(calcium_dynamics.Ca_i)'])) * timestep
Cai_time2peak = np.sum(log['dot(calcium_dynamics.Ca_i)'] > 1e-6 * 1) * timestep
print('Cai time to peak: ', Cai_time2peak)
CaSS_duration = np.sum((Ca_SS - min(Ca_SS)) > 0.1 * 1) * timestep
print('CaSS time to peak: ', CaSS_duration)
# Cai_time2peak = np.sum(log['calcium_dynamics.Ca_i'] > 1e-6 * 1) * timestep
# print('Cai time to peak: ', Cai_time2peak)

# Check the computation of time to peak
plt.figure()
plt.plot(log['calcium_dynamics.Ca_ss'])
# plt.plot((log['dot(calcium_dynamics.Ca_i)'] > 1e-5) * max(Ca_i), '--')
plt.plot((Ca_SS - min(Ca_SS)) > 0.1 * 1 * max(Ca_SS), '--')
# plt.xlim(80 * 10, 220 * 10)
fig_dir = '../figures/basic_sim/'
plt.savefig(fig_dir + 'test.pdf')

# Plot output
fig = modelling.figures.FigureStructure(figsize=(5, 2), gridspec=(1, 3),
                                        wspace=0.5)

current_list = ['membrane.V', 'calcium_dynamics.Ca_i',
                'calcium_dynamics.Ca_ss']
current_name = ['Vm', 'Cai', 'CaSS']

for i in range(len(current_list)):
    plot.add_single(fig.axs[0][i], log, current_list[i])
    fig.axs[0][i].spines[['right', 'top']].set_visible(False)
    fig.axs[0][i].set_ylabel(current_name[i])

fig.sharex(['Time (ms)'] * 3, [(0, 600), (0, 800), (0, 600)])

fig_dir = '../figures/basic_sim/'
fig.savefig(fig_dir + "ttp_sim_currents.pdf")