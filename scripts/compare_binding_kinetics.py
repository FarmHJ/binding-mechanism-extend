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
APmodel_name = sys.argv[1]
drug = sys.argv[2]
tuning_method = sys.argv[3]

# Define the range of drug concentration for a given drug
drug_conc_lib = modelling.DrugConcentrations()
drug_conc = drug_conc_lib.drug_concentrations[drug]['coarse']
repeats = 1000

# Define directories to save simulated data
data_dir = \
    '../simulation_data/kinetics_comparison/' + APmodel_name + '/' + \
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


Hill_coef_dir = '../simulation_data/kinetics_comparison/Hill_curves/' + \
    drug + '/'
Hill_coef_df = pd.read_csv(Hill_coef_dir + 'Hill_curves.csv')
Hill_coef = Hill_coef_df.loc[Hill_coef_df['protocol'] == 'Milnes']
Hill_coef = Hill_coef.values.tolist()[0][1:-1]

# Propagate to action potential
# Set AP model
if APmodel_name == 'Grandi':
    APmodel = '../math_model/AP_model/Grd-2010-IKr-SD.mmt'
elif APmodel_name == 'TTP':
    APmodel = '../math_model/AP_model/TTP-2006-IKr-SD.mmt'
elif APmodel_name == 'Lei':
    APmodel = '../math_model/AP_model/ORd-CiPA-Lei-SD.mmt'
APmodel, _, x = myokit.load(APmodel)
model_keys = modelling.ModelDetails().current_keys[APmodel_name]
current_key = model_keys['IKr']
current_head_key = current_key[:current_key.index('.')]
AP_model = modelling.Simulation(APmodel, current_head_key=current_head_key)

# Scale conductance value - method 1
# Same hERG peak as Grandi model
scale_df = pd.read_csv('../simulation_data/' + APmodel_name +
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

# Simulate AP of the ORd-SD model and the ORd-CS model
# Compute APD90
for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))
    log = AP_model.drug_simulation(
        drug, drug_conc[i], repeats_SD, save_signal=save_signal,
        log_var=log_var, abs_tol=abs_tol, rel_tol=rel_tol)
    log.save_csv(data_dir + 'SD_AP_' + str(drug_conc[i]) + '.csv')

    APD_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log[Vm_key, pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

    AP_trapping.append(log)
    APD_trapping.append(APD_trapping_pulse)

    reduction_scale = Hill_simulate(Hill_coef, drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats_CS,
        save_signal=save_signal,
        log_var=log_var, abs_tol=abs_tol, rel_tol=rel_tol)
    d2.save_csv(data_dir + 'CS_AP_' + str(drug_conc[i]) + '.csv')

    APD_conductance_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2[Vm_key, pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

    AP_conductance.append(d2)
    APD_conductance.append(APD_conductance_pulse)

    print('done concentration: ' + str(drug_conc[i]))

# Save simulated APD90 of both the ORd-SD model and the ORd-CS model
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

    # Run simulation for the ORd-SD model till steady state
    log = AP_model.drug_simulation(
        drug, drug_conc[i], repeats, save_signal=save_signal,
        log_var=log_var, abs_tol=abs_tol, rel_tol=rel_tol)

    # Compute APD90 of simulated AP
    APD_trapping_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(log[Vm_key, pulse], offset, 0.1)
        APD_trapping_pulse.append(apd90)

    APD_trapping.append(APD_trapping_pulse)

    # Run simulation for the ORd-CS model till steady state
    reduction_scale = Hill_simulate(Hill_coef, drug_conc[i])
    d2 = AP_model.conductance_simulation(
        base_conductance * reduction_scale, repeats,
        save_signal=save_signal,
        log_var=log_var, abs_tol=abs_tol, rel_tol=rel_tol)

    # Compute APD90 of simulated AP
    APD_conductance_pulse = []
    for pulse in range(save_signal):
        apd90 = AP_model.APD90(d2[Vm_key, pulse], offset, 0.1)
        APD_conductance_pulse.append(apd90)

    APD_conductance.append(APD_conductance_pulse)

    print('done concentration: ' + str(drug_conc[i]))

# Compute qNet
pulse_time = 2000
AP_model.protocol = myokit.pacing.blocktrain(pulse_time, 0.5)
save_signal = 1
qNet_current_list = [model_keys['ICaL'], current_key, model_keys['IKs'],
                     model_keys['IK1'], model_keys['Ito']]
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
    reduction_scale = Hill_simulate(Hill_coef, drug_conc[i])
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
APD_trapping = [max(i) for i in APD_trapping]
APD_conductance = [max(i) for i in APD_conductance]

# Save APD90 data
APD_trapping_df = pd.DataFrame(np.array(APD_trapping), columns=['APD'])
APD_trapping_df['drug concentration'] = drug_conc
APD_trapping_df['qNet'] = qNet_SD_arr
APD_trapping_df.to_csv(data_dir + 'SD_APD_fine.csv')
APD_conductance_df = pd.DataFrame(np.array(APD_conductance), columns=['APD'])
APD_conductance_df['drug concentration'] = drug_conc
APD_conductance_df['qNet'] = qNet_CS_arr
APD_conductance_df.to_csv(data_dir + 'CS_APD_fine.csv')
