# Check if the drug effect on Lei model is reasonable

import matplotlib
import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import pandas as pd
import pints
import sys

import modelling

# Define drug and protocol
drug = 'dofetilide'
protocol_name = 'Milnes'
pulse_time = 25e3
protocol = modelling.ProtocolLibrary().Milnes(pulse_time)

# Define the range of drug concentration for a given drug
drug_conc_lib = modelling.DrugConcentrations()
# drug_conc = drug_conc_lib.drug_concentrations[drug]['coarse']
repeats = 1000

#
# Check the inhibition of IKr for the Lei model and its state occupancy
#

# Define directories to save simulated data
Hill_coef_dir = '../../simulation_data/kinetics_comparison/Hill_curves/' + \
    drug + '/'
fig_dir = '../../figures/checking_figures/'
Hill_coef_dir2 = '../../simulation_data/kinetics_comparison/ORd/' + \
    drug + '/'
result_filename = 'Hill_curve.txt'

# Load IKr model
model = '../../math_model/current_model/lei2019_SD.mmt'
model, _, x = myokit.load(model)

model_keys = modelling.SDModelDetails().current_keys['Lei']
current_key = model_keys['IKr']
current_head_key = current_key[:current_key.index('.')]

current_model = modelling.Simulation(model, current_head_key=current_head_key)
current_model.protocol = protocol

# Define tolerance value
abs_tol = 1e-7
rel_tol = 1e-8
log_var = [model_keys['time'], model_keys['Vm'], current_key]
states = ['I', 'CI', 'O', 'C', 'Ibound', 'Obound', 'Cbound']
color_seq = ['#7e7e7e', '#986464', '#989864', '#986496',
             '#74a9cf', '#045a8d', '#2b8cbe']

# fig = modelling.figures.FigureStructure(figsize=(12, 12), gridspec=(3, 3),
#                                         height_ratios=[1] * 3)
plot = modelling.figures.FigurePlot()

# # Simulate IKr of the SD model for a range of drug concentrations
# # Extract the peak of IKr
# peaks = []
# peaks_alt = []
# for i in range(len(drug_conc)):
#     log = current_model.drug_simulation(
#         drug, drug_conc[i], repeats, abs_tol=abs_tol,
#         rel_tol=rel_tol)
# #     ax = fig.axs[int(i / 3)][i % 3]
# #     ax_twin = ax.twinx()
# #     ax.stackplot(
# #         log.time(), *[log['ikr.' + s] for s in states],
# #         labels=['ikr.' + s for s in states],
# #         colors=color_seq[:7], zorder=-10)
# #     ax.set_rasterization_zorder(0)
# #     ax_twin.plot(log.time(), log[current_key], 'k')
# #     ax_twin.set_ylim(0, 1.3)

# # fig.savefig(fig_dir + 'Lei_state_occupancy.pdf')
#     peak, _ = current_model.extract_peak(log, current_key)
#     peaks.append(peak[-1])

#     # repeats plus one
#     log = current_model.drug_simulation(
#         drug, drug_conc[i], repeats + 1, abs_tol=abs_tol,
#         rel_tol=rel_tol)
#     peak, _ = current_model.extract_peak(log, current_key)
#     peaks_alt.append(peak[-1])

# #
# # Compare the Hill curve from Lei-SD model and the original SD model
# #

# # Normalise drug response (peak current)
# peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))
# peaks_alt = (peaks_alt - min(peaks_alt)) / (max(peaks_alt) - min(peaks_alt))

# plt.figure()
# plt.plot(peaks)
# plt.plot(peaks_alt)

# estimates = np.loadtxt(Hill_coef_dir2 + result_filename, unpack=True)
# estimates = np.array(estimates)
# Hill_coef = estimates[:2]

# # Hill_coef_df = pd.read_csv(Hill_coef_dir + 'Hill_curves.csv')
# # Hill_coef = Hill_coef_df.loc[Hill_coef_df['protocol'] == 'Milnes']
# # Hill_coef = Hill_coef.values.tolist()[0][1:-1]

Hill_model = modelling.HillModel()
# optimiser = modelling.HillModelOpt(Hill_model)

# # Hill_curve = []
# Hill_curve2 = []
# for i in range(len(drug_conc)):
#     # reduction_scale = Hill_model.simulate(Hill_coef, drug_conc[i])
#     # Hill_curve.append(reduction_scale)

#     reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
#     Hill_curve2.append(reduction_scale)

# # plt.plot(Hill_curve)
# plt.plot(Hill_curve2)
# plt.savefig(fig_dir + 'alt_repeats_Hill.pdf')

#
# Propagate to action potential
#

# Set AP model
APmodel = '../../math_model/AP_model/ORd-CiPA-Lei-SD.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, current_head_key=current_head_key)

# Define protocol
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get(current_head_key + '.gKr').value()

# # Scale conductance value
# scale_df = pd.read_csv('../../simulation_data/Lei_conductance_scale.csv',
#                        index_col=[0], skipinitialspace=True)
# conductance_scale = scale_df.loc['AP_duration'].values[0]
# AP_model.model.set_value('tune.ikr_rescale', conductance_scale)

# # Define drug concentration range for steady state APD90 comparison between
# # models
drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']
# drug_conc = drug_conc[12:-2]
# repeats = 1000
# save_signal = 2
# offset = 50
# log_var = log_var + ['ikr.' + s for s in states]

# fig = modelling.figures.FigureStructure(figsize=(12, 8), gridspec=(2, 3),
#                                         height_ratios=[1] * 2)
# fig_CS = modelling.figures.FigureStructure(figsize=(12, 8), gridspec=(2, 3),
#                                            height_ratios=[1] * 2)

# for i in range(len(drug_conc)):
#     print('simulating concentration: ' + str(drug_conc[i]))

#     # Run simulation for the AP-SD model till steady state
#     log = AP_model.drug_simulation(
#         drug, drug_conc[i], repeats, save_signal=save_signal,
#         log_var=log_var, abs_tol=abs_tol, rel_tol=rel_tol)

#     ax = fig.axs[int(i / 3)][i % 3]
#     ax_twin = ax.twinx()
#     ax.stackplot(
#         log.time(), *[log['ikr.' + s] for s in states],
#         labels=['ikr.' + s for s in states],
#         colors=color_seq[:7], zorder=-10)
#     ax.set_rasterization_zorder(0)
#     ax_twin.plot(log.time(), log[model_keys['Vm']], 'k')

#     # Compute APD90 of simulated AP
#     log = log.fold(pulse_time)
#     APD_trapping_pulse = []
#     for pulse in range(save_signal):
#         apd90 = AP_model.APD90(log[model_keys['Vm'], pulse], offset, 0.1)
#         APD_trapping_pulse.append(apd90)

#     ax.text(0.8, 0.9, '{:.2f}'.format(max(APD_trapping_pulse)),
#             transform=ax.transAxes)
#     ax.set_xlim(0, 2000)
#     ax.set_ylim(0, 1)

#     # Run simulation for the AP-CS model till steady state
#     reduction_scale = Hill_model.simulate(Hill_coef, drug_conc[i])
#     d2 = AP_model.conductance_simulation(
#         base_conductance * reduction_scale, repeats,
#         save_signal=save_signal, log_var=log_var, abs_tol=abs_tol,
#         rel_tol=rel_tol)

#     ax = fig_CS.axs[int(i / 3)][i % 3]
#     ax_twin = ax.twinx()
#     ax.stackplot(
#         d2.time(), *[d2['ikr.' + s] for s in states],
#         labels=['ikr.' + s for s in states],
#         colors=color_seq[:7], zorder=-10)
#     ax.set_rasterization_zorder(0)
#     ax_twin.plot(d2.time(), d2[model_keys['Vm']], 'k')

#     # Compute APD90 of simulated AP
#     d2 = d2.fold(pulse_time)
#     APD_conductance_pulse = []
#     for pulse in range(save_signal):
#         apd90 = AP_model.APD90(d2[model_keys['Vm'], pulse], offset, 0.1)
#         APD_conductance_pulse.append(apd90)

#     ax.text(0.8, 0.9, '{:.2f}'.format(max(APD_conductance_pulse)),
#             transform=ax.transAxes)
#     ax.set_xlim(0, 2000)
#     ax.set_ylim(0, 1)

#     print('done concentration: ' + str(drug_conc[i]))

# fig.savefig(fig_dir + 'AP_state_occupancy_drugHill.pdf')
# fig_CS.savefig(fig_dir + 'AP_CS_state_occupancy_drugHill.pdf')

#
# Check Hill curves used in parameter space exploration
#

# data_filepath = '../../simulation_data/parameter_space_exploration/'
# # Load previously saved parameter space
# sample_filepath = data_filepath + 'Hill_curves.csv'

# param_space = []
# param_values_df = pd.read_csv(sample_filepath,
#                               header=[0, 1], index_col=[0],
#                               skipinitialspace=True)

# Hill_coef = param_values_df['Hill_curve'].values[0]
# param_values = param_values_df['param_values']
# drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']
# peaks = []
# for i in range(len(drug_conc)):
#     log = current_model.custom_simulation(
#         param_values, drug_conc[i], 1000,
#         log_var=['engine.time', 'ikr.IKr'],
#         abs_tol=1e-7, rel_tol=1e-8)
#     peak, _ = current_model.extract_peak(log, 'ikr.IKr')
#     peaks.append(peak[-1])

# peaks_norm = (peaks - min(peaks)) / (max(peaks) - min(peaks))

# Hill_curve_SD = []
# for i in range(len(drug_conc)):
#     reduction_scale = Hill_model.simulate(Hill_coef, drug_conc[i])
#     Hill_curve_SD.append(reduction_scale)

# plt.figure()
# plt.plot(peaks_norm, label='Lei')
# plt.plot(Hill_curve_SD, label='SD')
# plt.legend()
# plt.savefig(fig_dir + 'LeivsSD_Hill.pdf')

#
# Check state occupancy between Lei and SD model
#

# Define drug concentration range for steady state APD90 comparison between
# models
repeats = 1000
save_signal = 2
offset = 50
log_var = log_var + ['ikr.' + s for s in states]
drug_conc = drug_conc[13:17]

fig = modelling.figures.FigureStructure(figsize=(8, 12), gridspec=(4, 2),
                                        height_ratios=[1] * 4)
CiPA_model = '../../math_model/AP_model/ohara-cipa-2017.mmt'
CiPA_model, _, x = myokit.load(CiPA_model)
CiPAmodel = modelling.Simulation(CiPA_model, current_head_key=current_head_key)

# Define protocol
pulse_time = 1000
CiPAmodel.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)

for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))

    # Run simulation for the AP-SD model till steady state
    log = AP_model.drug_simulation(
        drug, drug_conc[i], repeats, save_signal=save_signal,
        log_var=log_var, abs_tol=abs_tol, rel_tol=rel_tol)

    ax = fig.axs[i][0]
    ax_twin = ax.twinx()
    ax.stackplot(
        log.time(), *[log['ikr.' + s] for s in states],
        labels=['ikr.' + s for s in states],
        colors=color_seq[:7], zorder=-10)
    ax.set_rasterization_zorder(0)
    ax_twin.plot(log.time(), log[model_keys['Vm']], 'k')

    # Compute APD90 of simulated AP
    log = log.fold(pulse_time)
    APD_Lei = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log[model_keys['Vm'], pulse], offset, 0.1)
        APD_Lei.append(apd90)

    ax.text(0.8, 0.9, '{:.2f}'.format(max(APD_Lei)),
            transform=ax.transAxes)
    ax.set_xlim(0, 2000)
    ax.set_ylim(0, 1)

    # Run simulation for the AP-CS model till steady state
    d2 = CiPAmodel.drug_simulation(
        drug, drug_conc[i], repeats, save_signal=save_signal,
        log_var=log_var, abs_tol=abs_tol, rel_tol=rel_tol)

    ax = fig.axs[i][1]
    ax_twin = ax.twinx()
    ax.stackplot(
        d2.time(), *[d2['ikr.' + s] for s in states],
        labels=['ikr.' + s for s in states],
        colors=color_seq[:7], zorder=-10)
    ax.set_rasterization_zorder(0)
    ax_twin.plot(d2.time(), d2[model_keys['Vm']], 'k')

    # Compute APD90 of simulated AP
    d2 = d2.fold(pulse_time)
    APD_SD = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2[model_keys['Vm'], pulse], offset, 0.1)
        APD_SD.append(apd90)

    ax.text(0.8, 0.9, '{:.2f}'.format(max(APD_SD)),
            transform=ax.transAxes)
    ax.set_xlim(0, 2000)
    ax.set_ylim(0, 1)

    print('done concentration: ' + str(drug_conc[i]))

fig.savefig(fig_dir + 'AP_state_occupancy_LeivsSD.pdf')
