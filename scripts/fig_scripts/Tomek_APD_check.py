# Calibrate the ionic conductance of the CS model from the SD model and
# compare the APD90 at steady state.
# Output:
# 1. IKr at various drug concentration simulated from the SD model.
# 2. Normalised peak IKr vs drug concentration in log scale.
# 3. Fitted Hill curve over peak IKr from the SD model.
# 4. IKr at various drug concentration simulated from conductance
#    calibrated hERG model.
# 5. Comparison of peak IKr between both models, i.e. the SD model
#    and conductance hERG model.
# 6. Comparison of IKr between both models.
# 7. 2 pulses of action potential simulated from both models (steady state).
# 8. APD90 of both pulses for both models (steady state).

import myokit
import numpy as np
import os
import pandas as pd
import sys

import modelling

# Define AP model, drug and tuning method
APmodel_name = 'Tomek'
drug = 'dofetilide'
tuning_method = 'hERG_peak'

# Define the range of drug concentration for a given drug
drug_conc_lib = modelling.DrugConcentrations()
repeats = 1000

# Define directories to save simulated data
data_dir = \
    '../../simulation_data/kinetics_comparison/' + APmodel_name + '/' + \
    tuning_method + '_match/' + drug + '/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

abs_tol = 1e-7
rel_tol = 1e-8


def Hill_simulate(params, drug_conc):

    Hills_coef = params[0]
    IC50 = params[1]

    response = IC50 / (np.power(drug_conc, Hills_coef) + IC50)

    return response


Hill_coef_dir = '../../simulation_data/kinetics_comparison/Hill_curves/' + \
    drug + '/'
Hill_coef_df = pd.read_csv(Hill_coef_dir + 'Hill_curves.csv')
Hill_coef = Hill_coef_df.loc[Hill_coef_df['protocol'] == 'Milnes']
Hill_coef = Hill_coef.values.tolist()[0][1:-1]

# Propagate to action potential
# Set AP model
model_details = modelling.ModelDetails()
model_filename = model_details.file_names[APmodel_name]
APmodel = '../../' + model_filename['AP_SD_path']

APmodel, _, x = myokit.load(APmodel)
model_keys = model_details.current_keys[APmodel_name]
current_key = model_keys['IKr']
current_head_key = current_key[:current_key.index('.')]
AP_model = modelling.Simulation(APmodel, current_head_key=current_head_key)

# Scale conductance value - method 1
# Same hERG peak as Grandi model
scale_df = pd.read_csv('../../simulation_data/' + APmodel_name +
                       '_conductance_scale.csv',
                       index_col=[0], skipinitialspace=True)
conductance_scale = scale_df.loc[tuning_method].values[0]
AP_model.model.set_value('tune.ikr_rescale', conductance_scale)

# Define current protocol
pulse_time = 1000
protocol_offset = 50
protocol = myokit.pacing.blocktrain(pulse_time, 0.5, offset=protocol_offset)
AP_model.protocol = protocol
base_conductance = APmodel.get(current_head_key + '.gKr').value()
protocol_duration = protocol.characteristic_time()

offset = 50
save_signal = 2
# For plotting purpose - so that EAD-like behaviour happens on the same pulse
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

time_key = model_keys['time']
Vm_key = model_keys['Vm']
log_var = [time_key, Vm_key, current_key]

# Define drug concentration range for steady state APD90 comparison between
# models
drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']
repeats = 1000
save_signal = 2

APD_conductance = []

for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))

    # Run simulation for the ORd-CS model till steady state
    reduction_scale = Hill_simulate(Hill_coef, drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats,
        save_signal=save_signal,
        log_var=log_var, abs_tol=abs_tol, rel_tol=rel_tol)
    d2.save_csv(data_dir + 'CS_AP_' + str(drug_conc[i]) + '_check.csv')

    # Compute APD90 of simulated AP
    Vm_signal = list(d2[Vm_key, 0])
    for pulse in range(1, save_signal):
        Vm_signal += list(d2[Vm_key, pulse])
    APD_conductance_pulse = AP_model.APD90_update(
        d2.time(), Vm_signal, offset, protocol_duration)

    APD_conductance.append(APD_conductance_pulse)

    print('done concentration: ' + str(drug_conc[i]))

# Compute APD90 with AP behaviour in alternating cycles
print(APD_conductance)
APD_conductance = [float('nan') if np.isnan(i).any() else max(i)
                   for i in APD_conductance]

# Save APD90 data
APD_conductance_df = pd.DataFrame(np.array(APD_conductance), columns=['APD'])
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df.to_csv(data_dir + 'CS_APD_fine.csv')

import matplotlib

fig_dir = \
    '../../figures/kinetics_comparison/' + APmodel_name + '/' + \
    tuning_method + '_match/' + drug + '/'

# Set up structure of the figure
fig = modelling.figures.FigureStructure(figsize=(8, 5), gridspec=(1, 1))
plot = modelling.figures.FigurePlot()

# Read files name of action potential data
CS_data_files = [f for f in os.listdir(data_dir) if
                 f.startswith('CS_AP_') and
                 f.endswith('_check.csv')]
conc_label_CS = [fname[6:-10] for fname in CS_data_files]
drug_conc_CS = [float(fname[6:-10]) for fname in CS_data_files]

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc_CS), key=lambda x:x[1])]
sort_ind_CS = [i[0] for i in
               sorted(enumerate(drug_conc_CS), key=lambda x:x[1])]
drug_conc = sorted(drug_conc_CS)
conductance_data_files = [CS_data_files[i] for i in sort_ind_CS]
conc_label = [conc_label_CS[i] for i in sort_ind]
conc_label = [r"$10^{:d}$".format(int(np.log10(float(i)))) if float(i) >= 1e3
              else i for i in conc_label]

# Initiate constants and variables
labels = [i + ' nM' for i in conc_label]
cmap = matplotlib.cm.get_cmap('viridis')

# Load action potential and APD data
conductance_AP_log = []
for i in range(len(conductance_data_files)):
    conductance_AP_log.append(myokit.DataLog.load_csv(
        data_dir + conductance_data_files[i]))

model_keys = modelling.ModelDetails().current_keys[APmodel_name]
current_key = model_keys['IKr']
Vm_key = model_keys['Vm']

CS_labelname = APmodel_name + '-CS model'

# Plot AP and IKr at various drug concentrations
plot.add_multiple_continuous(fig.axs[0][0], conductance_AP_log,
                             Vm_key, cmap=cmap,
                             labels=labels)
fig.axs[0][0].set_title(CS_labelname)
fig.axs[0][0].legend(handlelength=0.9, ncol=2, columnspacing=0.9,
                     loc=(1.04, 0.05))

# Adjust axes
# fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)])
# fig.sharey(['Voltage (mV)', 'Current (A/F)'],
#            axs=panel2, subgridspec=subgridspecs[0])
fig.savefig(fig_dir + 'APD_calculation_check.pdf')
