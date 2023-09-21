# Compare the AP, APD90 and qNet between the AP-Li-SD models and the AP-Li-CS models
# for a given drug at a chosen IKr tuning method
# Output:
# AP
# 1. 2 pulses of action potential simulated from both models (steady state).
# 2. APD90 of both pulses for both models (steady state).
# APD
# 1. APD90 for both models at various drug concentration.
# qNet
# 1. qNet for both models at various drug concentration.

import argparse
import myokit
import numpy as np
import os
import pandas as pd

import modelling

# Define AP model, drug and tuning method
parser = argparse.ArgumentParser(
    description="Comparison between AP-Li-SD model and the AP-Li-CS model")
parser.add_argument("APmodel", help="Name of AP model")
parser.add_argument("drug", help="Drug")
parser.add_argument("tuning_method", help="Method used to tune IKr")
parser.add_argument("-m", "--mode", default="APD-qNet", type=str,
                    help="Type of output: AP, APD, qNet APD-qNet")
args = parser.parse_args()

APmodel_name = args.APmodel
drug = args.drug
tuning_method = args.tuning_method
MODE = args.mode

# Define directories to save simulated data
data_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                        APmodel_name, tuning_method + '_match', drug + '_new')
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

# Get the Hill coef for Li-CS model
Hill_coef_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                             'Hill_curves', drug)
Hill_coef_df = pd.read_csv(os.path.join(Hill_coef_dir, 'Hill_curves.csv'))
Hill_coef = Hill_coef_df.loc[Hill_coef_df['protocol'] == 'Milnes']
Hill_coef = Hill_coef.values.tolist()[0][1:-1]

# Set AP model
model_keys = modelling.model_naming.model_current_keys[APmodel_name]
current_key = model_keys['IKr']
APsim = modelling.ModelSimController(APmodel_name)

# Scale conductance value
scale_df = pd.read_csv(os.path.join(modelling.RESULT_DIR, APmodel_name +
                       '_conductance_scale.csv'),
                       index_col=[0], skipinitialspace=True)
conductance_scale = scale_df.loc[tuning_method].values[0]
APsim.set_ikr_rescale(conductance_scale)

# Define the range of drug concentration for a given drug
if MODE == 'AP':
    drug_conc = modelling.SD_details.drug_concentrations[drug]['coarse']
else:
    drug_conc = modelling.SD_details.drug_concentrations[drug]['fine']

# # For plotting purpose - so that EAD-like behaviour happens on the same pulse
# if drug == 'dofetilide':
#     repeats_SD = 1000 + 1
#     repeats_CS = 1000
# else:
#     repeats_SD = 1000
#     repeats_CS = 1000

time_key = model_keys['time']
Vm_key = model_keys['Vm']
log_var = [time_key, Vm_key, current_key]

HillModel = modelling.HillModel()

# Simulate AP of the AP-SD model and the AP-CS model
if MODE in ["AP", "APD", "APD-qNet"]:
    APD_conductance = []
    APD_trapping = []
    save_signal = 2

    for i in range(len(drug_conc)):
        print('simulating concentration: ' + str(drug_conc[i]))
        APsim.set_SD_parameters(drug)
        APsim.set_conc(drug_conc[i])
        log = APsim.simulate(save_signal=save_signal, log_var=log_var)
        if MODE == "AP":
            log.save_csv(os.path.join(data_dir,
                                      'SD_AP_' + str(drug_conc[i]) + '.csv'))

        apd90 = APsim.APD90(log)
        APD_trapping.append(apd90)

        reduction_scale = HillModel.simulate(Hill_coef, drug_conc[i])
        APsim.set_conc(0)
        APsim.set_CS_parameter(reduction_scale)
        log = APsim.simulate(save_signal=save_signal, log_var=log_var)
        APsim.set_CS_parameter(1)
        if MODE == "AP":
            log.save_csv(os.path.join(data_dir,
                                      'CS_AP_' + str(drug_conc[i]) + '.csv'))

        apd90 = APsim.APD90(log)
        APD_conductance.append(apd90)

        print('done concentration: ' + str(drug_conc[i]))

# Save simulated APD90 of both the ORd-SD model and the ORd-CS model
if MODE == "AP":
    filename = 'APD_pulses' + str(int(save_signal)) + '.csv'
elif MODE in ["APD", "APD-qNet"]:
    filename = 'APD_fine.csv'
    APD_trapping = [float('nan') if np.isnan(i).any() else max(i)
                    for i in APD_trapping]
    APD_conductance = [float('nan') if np.isnan(i).any() else max(i)
                       for i in APD_conductance]
if MODE in ["AP", "APD", "APD-qNet"]:
    column_name = ['pulse ' + str(i) for i in range(save_signal)]
    APD_trapping_df = pd.DataFrame(APD_trapping, columns=column_name)
    APD_trapping_df['drug concentration'] = drug_conc
    APD_trapping_df.to_csv(os.path.join(data_dir, 'SD_' + filename))

    APD_conductance_df = pd.DataFrame(APD_conductance, columns=column_name)
    APD_conductance_df['drug concentration'] = drug_conc
    APD_conductance_df.to_csv(os.path.join(data_dir, 'CS_' + filename))

# Compute qNet
if MODE in ["qNet", "APD-qNet"]:
    pulse_time = 2000
    qNet_current_list = [model_keys['ICaL'], model_keys['INaL'], current_key,
                         model_keys['IKs'], model_keys['IK1'],
                         model_keys['Ito']]
    qNet_current_list = [i for i in qNet_current_list if i is not None]

    APsim.sim.set_protocol(myokit.pacing.blocktrain(period=pulse_time,
                                                    duration=0.5,
                                                    offset=50))
    APsim._cycle_length = pulse_time
    APsim.reset_parameters()
    APsim.update_initial_state()

    qNet_SD_arr = []
    qNet_CS_arr = []
    print('Computing qNet...')
    for i in range(len(drug_conc)):
        print('simulating concentration: ' + str(drug_conc[i]))

        # Run simulation for the AP-SD model till steady state
        APsim.set_SD_parameters(drug)
        APsim.set_conc(drug_conc[i])
        log = APsim.simulate(timestep=0.01,
                             log_var=log_var[:2] + qNet_current_list)

        qNet = APsim.qNet(log)
        qNet_SD_arr.append(qNet)

        # Run simulation for the AP-CS model till steady state
        reduction_scale = HillModel.simulate(Hill_coef, drug_conc[i])
        APsim.set_conc(0)
        APsim.set_CS_parameter(reduction_scale)
        log = APsim.simulate(timestep=0.01,
                             log_var=log_var[:2] + qNet_current_list)
        APsim.set_CS_parameter(1)

        qNet = APsim.qNet(log)
        qNet_CS_arr.append(qNet)

        print('done concentration: ' + str(drug_conc[i]))

    # Save qNet
    qNet_SD_df = pd.DataFrame(np.array(qNet_SD_arr), columns=['qNet'])
    qNet_SD_df['drug concentration'] = drug_conc
    qNet_SD_df.to_csv(os.path.join(data_dir, 'SD_qNet.csv'))
    qNet_CS_df = pd.DataFrame(np.array(qNet_CS_arr), columns=['qNet'])
    qNet_CS_df['drug concentration'] = drug_conc
    qNet_CS_df.to_csv(os.path.join(data_dir, 'CS_qNet.csv'))
