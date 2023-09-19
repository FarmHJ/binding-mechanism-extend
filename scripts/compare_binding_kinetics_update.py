# Compare the APD90 between the AP-Li-SD models and the AP-Li-CS models
# for a given drug at a chosen IKr tuning method
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

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import myokit
import numpy as np
import os
import pandas as pd
import sys

import modelling

# Define AP model, drug and tuning method
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("APmodel", help="Name of AP model")
parser.add_argument("drug", help="Drug")
parser.add_argument("tuning_method", help="Method used to tune IKr")
parser.add_argument("-m", "--mode", default="ALL", type=str,
                    help="Type of output: ALL, AP, APD or qNet")
args = vars(parser.parse_args())

APmodel_name = args['APmodel']
drug = args['drug']
tuning_method = args['tuning_method']
MODE = args["mode"]

# Define the range of drug concentration for a given drug
drug_conc_lib = modelling.DrugConcentrations()
if MODE == 'AP':
    drug_conc = drug_conc_lib.drug_concentrations[drug]['coarse']
else:
    drug_conc = drug_conc_lib.drug_concentrations[drug]['fine']

# Define directories to save simulated data
data_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                        APmodel_name, tuning_method + '_match', drug + '_new')
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

Hill_coef_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                             'Hill_curves', drug)
Hill_coef_df = pd.read_csv(os.path.join(Hill_coef_dir, 'Hill_curves.csv'))
Hill_coef = Hill_coef_df.loc[Hill_coef_df['protocol'] == 'Milnes']
Hill_coef = Hill_coef.values.tolist()[0][1:-1]

# Propagate to action potential
# Set AP model
model_keys = modelling.model_naming.model_current_keys[APmodel_name]
current_key = model_keys['IKr']
current_head_key = current_key[:current_key.index('.')]
APsim = modelling.ModelSimController(APmodel_name)

# Scale conductance value
scale_df = pd.read_csv(os.path.join(modelling.RESULT_DIR, APmodel_name +
                       '_conductance_scale.csv'),
                       index_col=[0], skipinitialspace=True)
conductance_scale = scale_df.loc[tuning_method].values[0]
APsim.set_ikr_rescale(conductance_scale)

# Define current protocol
pulse_time = 1000
protocol_offset = 50

save_signal = 2
# # For plotting purpose - so that EAD-like behaviour happens on the same pulse
# if drug == 'dofetilide':
#     repeats_SD = 1000 + 1
#     repeats_CS = 1000
# else:
#     repeats_SD = 1000
#     repeats_CS = 1000

AP_conductance = []
AP_trapping = []
APD_conductance = []
APD_trapping = []

time_key = model_keys['time']
Vm_key = model_keys['Vm']
log_var = [time_key, Vm_key, current_key]

# Simulate AP of the ORd-SD model and the ORd-CS model
# Compute APD90
HillModel = modelling.HillModel()
for i in range(len(drug_conc)):
    print('simulating concentration: ' + str(drug_conc[i]))
    APsim.set_SD_parameters(drug)
    APsim.set_conc(drug_conc[i])
    log = APsim.simulate(save_signal=save_signal, log_var=log_var)
    if MODE == "AP":
        log.save_csv(os.path.join(data_dir, 'SD_AP_' + str(drug_conc[i]) + '.csv'))

    apd90 = APsim.APD90(log)
    APD_trapping.append(apd90)

    reduction_scale = HillModel.simulate(Hill_coef, drug_conc[i])
    APsim.set_conc(0)
    APsim.set_CS_parameter(reduction_scale)
    log = APsim.simulate(save_signal=save_signal, log_var=log_var)
    APsim.set_CS_parameter(1)
    if MODE == "AP":
        log.save_csv(os.path.join(data_dir, 'CS_AP_' + str(drug_conc[i]) + '.csv'))

    apd90 = APsim.APD90(log)
    APD_conductance.append(apd90)

    print('done concentration: ' + str(drug_conc[i]))

# Save simulated APD90 of both the ORd-SD model and the ORd-CS model
column_name = ['pulse ' + str(i) for i in range(save_signal)]
APD_trapping_df = pd.DataFrame(APD_trapping, columns=column_name)
APD_trapping_df['drug concentration'] = drug_conc
if MODE == "AP":
    filename = 'SD_APD_pulses' + str(int(save_signal)) + '_n.csv'
elif MODE == "APD":
    filename = 'SD_APD_fine.csv'
APD_trapping_df.to_csv(os.path.join(data_dir, filename))

APD_conductance_df = pd.DataFrame(APD_conductance, columns=column_name)
APD_conductance_df['drug concentration'] = drug_conc
if MODE == "AP":
    filename = 'CS_APD_pulses' + str(int(save_signal)) + '_n.csv'
elif MODE == "APD":
    filename = 'CS_APD_fine.csv'
APD_conductance_df.to_csv(os.path.join(data_dir, filename))

# Compute qNet
if MODE == "qNet":
    pulse_time = 2000
    APsim.sim.set_protocol = myokit.pacing.blocktrain(pulse_time, 0.5, offset=50)
    APsim._cycle_length = pulse_time
    APsim.reset_parameters()
    control_log = AP_model.conductance_simulation(
        base_conductance, repeats, timestep=0.01)
    # Do you run drug free condition till steady state then add drug for qNet?

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
