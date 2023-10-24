import myokit
import numpy as np
import os

import modelling


###################
# O'Hara-CiPA model
###################
APsim = modelling.ModelSimController('ORd-Li', protocol_offset=10)
log = APsim.simulate()

AP_peak = APsim.extract_peak(log)
apd90 = APsim.APD90(log)
print('APD90: ', apd90)

#
# Grandi (2010) model
#
APsim = modelling.ModelSimController('Grandi', ikr_modified=False,
                                     protocol_offset=10)
log = APsim.simulate(log_var='all')

AP_peak = APsim.extract_peak(log)
apd90 = APsim.APD90(log)
resting_Vm = min(log[APsim.Vm_key][500 / 0.1:])
print('resting membrane potential: ', resting_Vm)
print('AP amplitude:', AP_peak - resting_Vm)
print('APD90: ', apd90)
print('maximal upstroke velocity: ',
      max(log['dot(membrane_potential.V_m)']))

# Plot output
fig = modelling.figures.FigureStructure(figsize=(5, 3), gridspec=(2, 1))
plot = modelling.figures.FigurePlot()

plot.add_single(fig.axs[0][0], log, 'interface.I_app')
plot.add_single(fig.axs[1][0], log, 'membrane_potential.V_m')
# plot.state_occupancy_plot(fig.axs[2][0], log, AP_model)

fig.sharex(['Time (ms)'], [(0, APsim._cycle_length)])
fig.sharey(['Current\n(A/F)', 'Voltage\n(mV)'])

fig.savefig(os.path.join(modelling.FIG_DIR, 'basic_sim', "grd_sim.pdf"))

# Plot output
fig = modelling.figures.FigureStructure(figsize=(8, 5), gridspec=(6, 2),
                                        height_ratios=[1] * 6, hspace=0.2)

current_list = ['membrane_potential.V_m', 'Ca_Concentrations.Ca_i',
                'Ca_Concentrations.Ca_j', 'Ca_Concentrations.Ca_sl',
                'I_Na.I_Na', 'I_Ca.I_Catot', 'I_Kr.I_kr', 'I_Ks.I_ks',
                'I_to.I_to', 'I_Ki.I_ki', 'I_NCX.I_ncx', 'I_NaK.I_nak']
current_name = ['Vm', '[Ca]i', '[Ca]jct', '[Ca]SL',
                'I_Na', 'I_CaL', 'I_Kr', 'I_Ks',
                'I_to', 'I_K1', 'I_NCX', 'I_NaK']

for i in range(len(current_list)):
    row = i % 6
    col = i // 6
    plot.add_single(fig.axs[row][col], log, current_list[i])
    fig.axs[row][col].spines[['right', 'top']].set_visible(False)
    fig.axs[row][col].set_ylabel(current_name[i])
# plot.state_occupancy_plot(fig.axs[2][0], log, AP_model)

fig.sharex(['Time (ms)'] * 2, [(0, 500)] * 2)
# fig.sharey(['Current\n(A/F)', 'Voltage\n(mV)'])

fig.savefig(os.path.join(modelling.FIG_DIR, 'basic_sim',
                         "grd_sim_currents.pdf"))

#
# ten Tusscher (2006) model
#
timestep = 0.01
APsim = modelling.ModelSimController('TTP', ikr_modified=False,
                                     protocol_offset=100)
log = APsim.simulate(timestep=0.01, log_var='all')

AP_peak = APsim.extract_peak(log)
apd90 = APsim.apd90(log)

AP_Vm = log[APsim.Vm_key]
APD90_v = min(AP_Vm) + 0.1 * (max(AP_Vm) - min(AP_Vm))

# Action potential measures
resting_Vm = min(log[APsim.Vm_key][500 * 10:])
print('###############')
print('resting membrane potential: ', resting_Vm)
print('AP amplitude:', AP_peak - resting_Vm)
print('APD90: ', apd90)
print('maximal upstroke velocity: ',
      max(log['dot(membrane.V)']))

# Plot TTP's action potential
fig = modelling.figures.FigureStructure(figsize=(8, 3), gridspec=(2, 1))
plot = modelling.figures.FigurePlot()

plot.add_single(fig.axs[0][0], log, 'membrane.i_Stim')
plot.add_single(fig.axs[1][0], log, 'membrane.V')
# plot.state_occupancy_plot(fig.axs[2][0], log, AP_model)

fig.sharex(['Time (ms)'], [(0, APsim._cycle_length)])
fig.sharey(['Current\n(A/F)', 'Voltage\n(mV)'])

fig.savefig(os.path.join(modelling.FIG_DIR, 'basic_sim', "ttp_sim.pdf"))

# Calcium dynamics measures
Ca_SS = log['calcium_dynamics.Ca_ss']
CaSS_peak_amplitude = max(Ca_SS) - min(Ca_SS)
print('CaSS peak amplitude: ', CaSS_peak_amplitude)
Ca_i = log['calcium_dynamics.Ca_i']
Cai_peak_amplitude = max(Ca_i) - min(Ca_i)
print('Cai peak amplitude: ', Cai_peak_amplitude)
# CaSS_time2peak = (np.argmax(Ca_SS) -
#                   np.argmax(log['dot(calcium_dynamics.Ca_ss)'])) * timestep
CaSS_time2peak = np.sum(log['dot(calcium_dynamics.Ca_ss)'] > 0.01 * 1) * \
    timestep
print('CaSS time to peak: ', CaSS_time2peak)
# Cai_time2peak = (np.argmax(Ca_i) -
#                  np.argmax(log['dot(calcium_dynamics.Ca_i)'])) * timestep
Cai_time2peak = np.sum(log['dot(calcium_dynamics.Ca_i)'] > 1e-6 * 1) * timestep
print('Cai time to peak: ', Cai_time2peak)
CaSS_duration = np.sum((Ca_SS - min(Ca_SS)) > 0.1 * 1) * timestep
print('CaSS duration: ', CaSS_duration)
Cai_duration = np.sum((Ca_i - min(Ca_i)) > 1e-5 * 1) * timestep
print('Cai duration: ', Cai_duration)

# # Check the computation of time to peak
# plt.figure()
# plt.plot(log.time(), log['calcium_dynamics.Ca_i'])
# # plt.plot((log['dot(calcium_dynamics.Ca_i)'] > 1e-5) * max(Ca_i), '--')
# plt.plot(log.time(), ((Ca_i - min(Ca_i)) > 1e-5 * 1) * max(Ca_i), '--')
# # plt.xlim(80 * 10, 220 * 10)
# fig_dir = '../figures/basic_sim/'
# plt.savefig(fig_dir + 'test.pdf')

# # Plot output
# fig = modelling.figures.FigureStructure(figsize=(8, 2), gridspec=(1, 3),
#                                         wspace=0.5)

# current_list = ['membrane.V', 'calcium_dynamics.Ca_i',
#                 'calcium_dynamics.Ca_ss']
# current_name = ['Vm', 'Cai', 'CaSS']

# for i in range(len(current_list)):
#     fig.axs[0][i].plot(log.time(), log[current_list[i]], '--')
#     fig.axs[0][i].spines[['right', 'top']].set_visible(False)
#     fig.axs[0][i].set_ylabel(current_name[i])

# fig.sharex(['Time (ms)'] * 3, [(0, 600), (0, 800), (0, 600)])
# fig.axs[0][0].set_ylim(-90, 35)
# fig.axs[0][1].set_ylim(0, 0.001)
# fig.axs[0][2].set_ylim(0, 2)

# fig_dir = '../figures/basic_sim/'
# fig.savefig(fig_dir + "ttp_sim_currents.svg", format='svg')

# Plot measure calculation
fig = modelling.figures.FigureStructure(figsize=(8, 6), gridspec=(3, 3),
                                        height_ratios=[1, 1, 1],
                                        wspace=0.4, hspace=0.5)

# APD90
fig.axs[0][0].plot(log.time(), log['membrane.V'])
fig.axs[0][0].hlines(y=APD90_v, xmin=APsim.protocol_offset,
                     xmax=APsim.protocol_offset + apd90,
                     color='red', linestyle='--')
fig.axs[0][0].vlines(x=APsim.protocol_offset + apd90, ymin=-90, ymax=40,
                     color='red', linestyle='--')

fig.axs[0][0].set_title('APD90')

# resting membrane potential
fig.axs[0][1].plot(log.time(), log['membrane.V'])
fig.axs[0][1].hlines(y=resting_Vm, xmin=APsim.protocol_offset,
                     xmax=APsim._cycle_length, color='red', linestyle='--')
fig.axs[0][1].set_title('resting Vm')
# maximal upstroke velocity
fig.axs[0][2].plot(log.time(), log['membrane.V'])
axs_twin = fig.axs[0][2].twinx()
axs_twin.plot(log.time(), log['dot(membrane.V)'], 'r--')
fig.axs[0][2].set_title('max upstroke velocity')

for i in range(3):
    fig.axs[0][i].set_xlabel('Time (ms)')
    fig.axs[0][i].set_ylabel('Voltage (mV)')
    fig.axs[0][i].set_ylim(-90, 35)

# Ca_i peak amplitude
fig.axs[1][0].plot(log.time(), log['calcium_dynamics.Ca_i'])
fig.axs[1][0].axhline(y=min(Ca_i), color='red', linestyle='--')
fig.axs[1][0].axhline(y=max(Ca_i), color='red', linestyle='--')
fig.axs[1][0].set_title('peak amplitude')
# Ca_i time to peak
fig.axs[1][1].plot(log.time(), log['calcium_dynamics.Ca_i'])
fig.axs[1][1].plot(log.time(),
                   (log['dot(calcium_dynamics.Ca_i)'] > 1e-6) * max(Ca_i),
                   'r--')
fig.axs[1][1].set_xlim((APsim.protocol_offset - 20),
                       (APsim.protocol_offset + Cai_time2peak + 50))
fig.axs[1][1].set_title('time to peak')
# Ca_i duration
fig.axs[1][2].plot(log.time(), log['calcium_dynamics.Ca_i'])
fig.axs[1][2].plot(log.time(), ((Ca_i - min(Ca_i)) > 1e-5) * max(Ca_i), 'r--')
fig.axs[1][2].set_xlim((APsim.protocol_offset - 20),
                       (APsim.protocol_offset + Cai_duration + 100))
fig.axs[1][2].set_title('duration')

for i in range(3):
    fig.axs[1][i].set_xlabel('Time (ms)')
    fig.axs[1][i].set_ylabel('Ca_i (mM)')
    fig.axs[1][i].set_ylim(0, 0.001)

# Ca_SS peak amplitude
fig.axs[2][0].plot(log.time(), log['calcium_dynamics.Ca_ss'])
fig.axs[2][0].axhline(y=min(Ca_SS), color='red', linestyle='--')
fig.axs[2][0].axhline(y=max(Ca_SS), color='red', linestyle='--')
fig.axs[2][0].set_xlim((APsim.protocol_offset - 20),
                       (APsim.protocol_offset + CaSS_duration + 50))
fig.axs[2][0].set_title('peak amplitude')
# Ca_SS time to peak
fig.axs[2][1].plot(log.time(), log['calcium_dynamics.Ca_ss'])
fig.axs[2][1].plot(log.time(),
                   (log['dot(calcium_dynamics.Ca_ss)'] > 0.01) * max(Ca_SS),
                   'r--')
fig.axs[2][1].set_xlim((APsim.protocol_offset - 10),
                       (APsim.protocol_offset + CaSS_time2peak + 20))
fig.axs[2][1].set_title('time to peak')
# Ca_SS duration
fig.axs[2][2].plot(log.time(), log['calcium_dynamics.Ca_ss'])
fig.axs[2][2].plot(log.time(), ((Ca_SS - min(Ca_SS)) > 0.1) * max(Ca_SS),
                   'r--')
fig.axs[2][2].set_xlim((APsim.protocol_offset - 20),
                       (APsim.protocol_offset + CaSS_duration + 50))
fig.axs[2][2].set_title('duration')

for i in range(3):
    fig.axs[2][i].set_xlabel('Time (ms)')
    fig.axs[2][i].set_ylabel('Ca_SS (mM)')
    fig.axs[2][i].set_ylim(-0.05, 2.05)

# fig.sharex(['Time (ms)'] * 3, [(0, 600), (0, 800), (0, 600)])
# fig.axs[0][0].set_ylim(-90, 35)
# fig.axs[0][1].set_ylim(0, 0.001)
# fig.axs[0][2].set_ylim(0, 2)

fig.savefig(os.path.join(modelling.FIG_DIR, 'basic_sim',
                         "ttp_sim_APmeasure.svg"), format='svg')
