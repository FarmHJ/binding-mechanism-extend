# Check that the Hill curve from the Lei model is the same as from the Li model

import matplotlib
import myokit
import numpy as np
import os
import pandas as pd
import pints

import modelling


# Define drug and protocol
drug = 'dofetilide'
protocol_name = 'Milnes'
pulse_time = 25e3
protocol = modelling.ProtocolLibrary().Milnes(pulse_time)

parameter_name = ['ikr.kmax', 'ikr.ku', 'ikr.halfmax', 'ikr.hill', 'ikr.vhalf']
parameter_dir = '../../Lei-SD-fits/Milnes-data-fits/'

# for fname in glob.glob(parameter_dir +
#                       'Drug-*-lei-model12-protocol-Milnes-fit-RMSE.txt'):
drug_params = np.loadtxt(parameter_dir +
                         'Drug-dofetilide-lei-model12-protocol-Milnes-fit-RMSE.txt',
                         unpack=True)
drug_params = np.array(drug_params)
param_df = pd.DataFrame(drug_params, index=['Kmax', 'Ku', 'EC50', 'N', 'Vhalf']).T
param_df['Vhalf'] = -1 * param_df['Vhalf']
print(param_df)

# Define the range of drug concentration for a given drug
drug_conc_lib = modelling.DrugConcentrations()
drug_conc = drug_conc_lib.drug_concentrations[drug]['coarse']
repeats = 1000

# Define directories to save simulated data
data_dir = '../simulation_data/kinetics_comparison/ORd-Lei/AP_duration_match/' + \
    drug + '/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
result_filename = 'Hill_curve_Lei.txt'

# Load IKr model
model = '../math_model/current_model/lei2019_SD.mmt'
model, _, x = myokit.load(model)

model_keys = modelling.ModelDetails().current_keys['Lei']
current_key = model_keys['IKr']
current_head_key = current_key[:current_key.index('.')]

current_model = modelling.Simulation(model, current_head_key=current_head_key)
current_model.protocol = protocol

# Define tolerance value
abs_tol = 1e-7
rel_tol = 1e-8
log_var = [model_keys['time'], model_keys['Vm'], current_key]

# Simulate IKr of the SD model for a range of drug concentrations
# Extract the peak of IKr
peaks = []
SD_peaks = []
for i in range(len(drug_conc)):
    log = current_model.custom_simulation(
        param_df, drug_conc[i], repeats, log_var=log_var, abs_tol=abs_tol,
        rel_tol=rel_tol)
    peak, _ = current_model.extract_peak(log, current_key)
    peaks.append(peak[-1])

    log.save_csv(data_dir + 'SD_current_' + str(drug_conc[i]) + '.csv')

# Normalise drug response (peak current)
peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))

# Fit drug response to Hill curve
Hill_model = modelling.HillModel()
optimiser = modelling.HillModelOpt(Hill_model)
if not os.path.isfile(data_dir + result_filename):
    estimates, _ = optimiser.optimise(drug_conc, peaks)
    with open(data_dir + result_filename, 'w') as f:
        for x in estimates:
            f.write(pints.strfloat(x) + '\n')
else:
    estimates = np.loadtxt(data_dir + result_filename, unpack=True)
    estimates = np.array(estimates)

# Compare peak current
base_conductance = model.get(current_head_key + '.gKr').value()
current_model.current_head = current_model.model.get(current_head_key)
for i in range(len(drug_conc)):

    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = current_model.conductance_simulation(
        base_conductance * reduction_scale, repeats, log_var=log_var,
        abs_tol=abs_tol, rel_tol=rel_tol)

    d2.save_csv(data_dir + 'CS_current_' + str(drug_conc[i]) + '.csv')

# Plot results
fig = modelling.figures.FigureStructure(figsize=(8, 2), gridspec=(1, 2),
                                        width_ratios=[2.5, 1.2],
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(1, 2), (1, 1)]
subgs = []
for i in range(2):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.08,
                                       hspace=0.08))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
         subgridspecs[k][1])] for i in range(subgridspecs[k][0])]
       for k in range(len(subgs))]

panel1 = axs[0]

# Read files name of IKr data
SD_data_files = [f for f in os.listdir(data_dir) if
                 f.startswith('SD_current_')]
CS_data_files = [f for f in os.listdir(data_dir) if
                 f.startswith('CS_current_')]
drug_conc_SD = [float(fname[11:-4]) for fname in SD_data_files]
drug_conc_CS = [float(fname[11:-4]) for fname in CS_data_files]
conc_label_SD = [fname[11:-4] for fname in SD_data_files]

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc_SD), key=lambda x:x[1])]
sort_ind_CS = [i[0] for i in
               sorted(enumerate(drug_conc_CS), key=lambda x:x[1])]
trapping_data_files = [SD_data_files[i] for i in sort_ind]
conductance_data_files = [CS_data_files[i] for i in sort_ind_CS]
conc_label = [conc_label_SD[i] for i in sort_ind]
conc_label = [r"$10^{:d}$".format(int(np.log10(float(i)))) if float(i) >= 1e3
              else i for i in conc_label]

# Initiate constants and variables
labels = [i + ' nM' for i in conc_label]
cmap = matplotlib.cm.get_cmap('viridis')

# Load IKr data
trapping_hERG_log = []
conductance_hERG_log = []
for i in range(len(trapping_data_files)):
    trapping_hERG_log.append(myokit.DataLog.load_csv(
        data_dir + trapping_data_files[i]))
    conductance_hERG_log.append(myokit.DataLog.load_csv(
        data_dir + conductance_data_files[i]))

# hERG_trapping_plot = [e for i, e in enumerate(trapping_hERG_log)
#                       if i not in chosen_conc_ind]
# hERG_conductance_plot = [e for i, e in enumerate(conductance_hERG_log)
#                          if i not in chosen_conc_ind]
hERG_trapping_plot = trapping_hERG_log
hERG_conductance_plot = conductance_hERG_log

# Plot Ikr
plot.add_multiple(panel1[0][0], hERG_trapping_plot, current_head_key + '.IKr',
                  labels=labels, color=cmap)
plot.add_multiple(panel1[0][1], hERG_conductance_plot,
                  current_head_key + '.IKr',
                  labels=labels, color=cmap)

# Adjust figure details
panel1[0][1].legend(handlelength=0.9, ncol=2, columnspacing=0.9)
panel1[0][0].set_title('SD model')
panel1[0][1].set_title('CS model')
fig.sharex(['Time (s)'] * 2, [(0, pulse_time)] * 2,
           axs=panel1, subgridspec=subgridspecs[0])
fig.sharey(['Current (A/F)'],
           axs=panel1, subgridspec=subgridspecs[0])
panel1[0][0].spines[['right', 'top']].set_visible(False)
panel1[0][1].spines[['right', 'top']].set_visible(False)
fig.adjust_ticks(panel1[0][0], pulse_time)
fig.adjust_ticks(panel1[0][1], pulse_time)

# Plot fitting result of Hill curve
max_grid = np.ceil(np.log(drug_conc[-1]))
conc_grid = np.arange(-3, max_grid + 1, 0.5)

panel2 = axs[1]
panel2[0][0].plot(np.log(drug_conc[1:]), peaks[1:], 'o', label='peak current')
panel2[0][0].plot(conc_grid, Hill_model.simulate(estimates[:2],
                                                 np.exp(conc_grid)),
                  'k', label='fitted Hill eq')
panel2[0][0].set_xlabel('Drug concentration (log)')
panel2[0][0].set_ylabel('Normalised peak current')

fig_dir = '../figures/kinetics_comparison/Lei_reparameterise/'
fig.savefig(fig_dir + 'lei_' + drug + '.pdf')

#
# Propagate to action potential
#

# Set AP model
APmodel = '../math_model/AP_model/ORd-CiPA-Lei-SD.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, current_head_key=current_head_key)

scale_df = pd.read_csv('../simulation_data/' + 'Lei_conductance_scale.csv',
                       index_col=[0], skipinitialspace=True)
conductance_scale = scale_df.loc['AP_duration'].values[0]
AP_model.model.set_value('tune.ikr_rescale', conductance_scale)

# Define current protocol
pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get(current_head_key + '.gKr').value()

offset = 50
save_signal = 2
# Use different repeats for plotting purpose - so that EAD-like behaviour
# happens on the same pulse
if drug == 'dofetilide':
    repeats_SD = 1000 + 1
    repeats_CS = 1000
else:
    repeats_SD = 1000
    repeats_CS = 1000   

AP_conductance = []
AP_trapping = []
APD_conductance = []
APD_trapping = []

# Simulate AP of the AP-SD model and the AP-CS model
# Compute APD90
for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))
    log = AP_model.custom_simulation(
        param_df, drug_conc[i], repeats_SD, save_signal=save_signal,
        log_var=log_var, abs_tol=abs_tol, rel_tol=rel_tol)
    log.save_csv(data_dir + 'SD_AP_' + str(drug_conc[i]) + '.csv')

    APD_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log[model_keys['Vm'], pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

    AP_trapping.append(log)
    APD_trapping.append(APD_trapping_pulse)

    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats_CS,
        save_signal=save_signal, log_var=log_var, abs_tol=abs_tol,
        rel_tol=rel_tol)
    d2.save_csv(data_dir + 'CS_AP_' + str(drug_conc[i]) + '.csv')

    APD_conductance_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2[model_keys['Vm'], pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

    AP_conductance.append(d2)
    APD_conductance.append(APD_conductance_pulse)

    print('done concentration: ' + str(drug_conc[i]))

# Save simulated APD90 of both the AP-SD model and the AP-CS model
column_name = ['pulse ' + str(i) for i in range(save_signal)]
APD_trapping_df = pd.DataFrame(APD_trapping, columns=column_name)
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df.to_csv(data_dir + 'SD_APD_pulses' +
                       str(int(save_signal)) + '.csv')
APD_conductance_df = pd.DataFrame(APD_conductance, columns=column_name)
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df.to_csv(data_dir + 'CS_APD_pulses' +
                          str(int(save_signal)) + '.csv')

# Define drug concentration range for steady state APD90 comparison between
# models
drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']
repeats = 1000
save_signal = 2

APD_conductance = []
APD_trapping = []

for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))

    # Run simulation for the AP-SD model till steady state
    log = AP_model.custom_simulation(
        param_df, drug_conc[i], repeats, save_signal=save_signal,
        log_var=log_var, abs_tol=abs_tol, rel_tol=rel_tol)

    # Compute APD90 of simulated AP
    APD_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log[model_keys['Vm'], pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

    APD_trapping.append(APD_trapping_pulse)

    # Run simulation for the AP-CS model till steady state
    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats,
        save_signal=save_signal, log_var=log_var, abs_tol=abs_tol,
        rel_tol=rel_tol)

    # Compute APD90 of simulated AP
    APD_conductance_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2[model_keys['Vm'], pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

    APD_conductance.append(APD_conductance_pulse)

    print('done concentration: ' + str(drug_conc[i]))

# Compute qNet
pulse_time = 2000
AP_model.protocol = myokit.pacing.blocktrain(pulse_time, 0.5)
save_signal = 1
qNet_current_list = [model_keys['ICaL'], model_keys['INaL'], current_key,
                     model_keys['IKs'], model_keys['IK1'], model_keys['Ito']]
qNet_current_list = [i for i in qNet_current_list if i is not None]
control_log = AP_model.conductance_simulation(
    base_conductance, repeats, timestep=0.01, abs_tol=abs_tol,
    rel_tol=rel_tol)

qNet_SD_arr = []
qNet_CS_arr = []
print('Computing qNet...')
for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))

    # Run simulation for the ORd-SD model till steady state
    log = AP_model.drug_simulation(
        drug, drug_conc[i], repeats + save_signal, timestep=0.01,
        log_var=log_var[:2] + qNet_current_list, abs_tol=abs_tol,
        rel_tol=rel_tol, set_state=control_log)

    inet = 0
    for c in qNet_current_list:
        inet += log[c]  # pA/pF

    qNet = np.trapz(inet, x=log.time()) * 1e-3  # pA/pF*s
    qNet_SD_arr.append(qNet)

    # Run simulation for the ORd-CS model till steady state
    reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats + save_signal,
        timestep=0.01, log_var=log_var[:2] + qNet_current_list,
        abs_tol=abs_tol, rel_tol=rel_tol, set_state=control_log)

    inet = 0
    for c in qNet_current_list:
        inet += d2[c]  # pA/pF

    qNet = np.trapz(inet, x=d2.time()) * 1e-3  # pA/pF*s
    qNet_CS_arr.append(qNet)

    print('done concentration: ' + str(drug_conc[i]))

# Compute APD90 with AP behaviour in alternating cycles
APD_trapping = [float('nan') if np.isnan(i).any() else max(i)
                for i in APD_trapping]
APD_conductance = [float('nan') if np.isnan(i).any() else max(i)
                   for i in APD_conductance]

# Save APD90 data
APD_trapping_df = pd.DataFrame(np.array(APD_trapping), columns=['APD'])
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df['qNet'] = qNet_SD_arr
APD_trapping_df.to_csv(data_dir + 'SD_APD_fine.csv')
APD_conductance_df = pd.DataFrame(np.array(APD_conductance), columns=['APD'])
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df['qNet'] = qNet_CS_arr
APD_conductance_df.to_csv(data_dir + 'CS_APD_fine.csv')
