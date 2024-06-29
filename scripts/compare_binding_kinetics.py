# Compare the AP, APD90 and qNet between the AP-Li-SD models and
# the AP-Li-CS models for a given drug at a chosen IKr tuning method
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

# TODO: pace model before drug

# Define AP model, drug and tuning method
parser = argparse.ArgumentParser(
    description="Comparison between AP-IKr-SD model and the AP-IKr-CS model")
parser.add_argument("APmodel", help="Name of AP model")
parser.add_argument("drug", help="Drug")
parser.add_argument("--ikr_tuning", default='AP_duration',
                    choices=['hERG_peak', 'hERG_flux', 'AP_duration'],
                    help="Method used to tune IKr")
parser.add_argument("-m", "--mode", default="APD-qNet", type=str,
                    choices=['AP', 'APD', 'qNet', 'APD-qNet'],
                    help="Type of output: AP, APD, qNet, APD-qNet")
parser.add_argument("-c", "--conc", default="conc", type=str,
                    choices=['dimless', 'conc'],
                    help='Choose to convert concentration to dimensionless')
args = parser.parse_args()

APmodel_name = args.APmodel
drug = args.drug
ikr_tuning = args.ikr_tuning
MODE = args.mode

# Define directories to save simulated data
data_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                        APmodel_name, f'{ikr_tuning}_match', drug)
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

# Get the Hill coef for Li-CS model
# Update according to file format
if APmodel_name == 'ORd-Lei':
    ikr_model = 'Lei'
else:
    ikr_model = 'Li'

# Hill_coef = modelling.BindingParameters().load_published_Hill_eq(
#     drug, ikr_model=ikr_model, channel='IKr')
Hill_coef = modelling.BindingParameters().load_Hill_eq(
    drug, ikr_model=ikr_model)

# Set AP model
model_keys = modelling.model_naming.model_current_keys[APmodel_name]
current_key = model_keys['IKr']
APsim = modelling.ModelSimController(APmodel_name)

# Scale conductance value
if APmodel_name != 'ORd-Li':
    APsim.set_ikr_rescale_method(ikr_tuning)

# Define the range of drug concentration for a given drug
if MODE == 'AP':
    drug_conc = modelling.SD_details.drug_concentrations[drug]['coarse']
    if args.conc == 'dimless':
        dimless_input = np.array([10 ** i for i in
                                  np.linspace(-10, np.log10(0.9), 10)])
else:
    drug_conc = modelling.SD_details.drug_concentrations[drug]['fine']
    if args.conc == 'dimless':
        dimless_input = np.array([10 ** i for i in
                                  np.linspace(-10, np.log10(0.9), 20)])
    # if drug == 'dofetilide':
    #     drug_conc = 10.0**np.linspace(3, 6, 10)
    # elif drug == 'verapamil':
    #     drug_conc = 10.0**np.linspace(5, 8, 10)


def from_dimless(dimless_conc, n, halfmax):
    expon = np.log(dimless_conc / (1 - dimless_conc)) / n
    return np.power(halfmax, 1 / n) * np.exp(expon)


HillModel = modelling.HillModel()

# Simulate AP of the AP-SD model and the AP-CS model
if MODE in ["AP", "APD", "APD-qNet"]:
    APD_conductance = []
    APD_trapping = []
    save_signal = 2

    APsim.set_SD_parameters(drug, ikr_model=ikr_model)
    if args.conc == 'dimless':
        SD_params = modelling.BindingParameters(
            ikr_model=ikr_model).load_SD_parameters(drug)
        drug_conc = from_dimless(dimless_input, SD_params['n'].values,
                                 SD_params['halfmax'].values)
    for i in range(len(drug_conc)):
        print('simulating concentration: ' + str(drug_conc[i]))

        APsim.set_conc(drug_conc[i])
        log = APsim.simulate(save_signal=save_signal)
        if MODE == "AP":
            if args.conc == 'dimless':
                log.save_csv(os.path.join(data_dir,
                                          f'SD_AP_dimless_{drug_conc[i]}.csv'))
            else:
                log.save_csv(os.path.join(data_dir,
                                          f'SD_AP_conc_{drug_conc[i]}.csv'))

        apd90 = APsim.APD90(log)
        APD_trapping.append(apd90)

        reduction_scale = HillModel.simulate(Hill_coef, drug_conc[i])
        APsim.set_conc(0)
        APsim.set_CS_parameter(reduction_scale)
        log = APsim.simulate(save_signal=save_signal)
        APsim.set_CS_parameter(1)
        if MODE == "AP":
            if args.conc == 'dimless':
                log.save_csv(os.path.join(data_dir,
                                          f'CS_AP_dimless_{drug_conc[i]}.csv'))
            else:
                log.save_csv(os.path.join(data_dir,
                                          f'CS_AP_conc_{drug_conc[i]}.csv'))

        apd90 = APsim.APD90(log)
        APD_conductance.append(apd90)

        print('done concentration: ' + str(drug_conc[i]))

# Save simulated APD90 of both the ORd-SD model and the ORd-CS model
if MODE == "AP":
    if args.conc == 'dimless':
        filename = f'APD_pulses{int(save_signal)}_dimless.csv'
    else:
        filename = f'APD_pulses{int(save_signal)}.csv'
elif MODE in ["APD", "APD-qNet"]:
    if args.conc == 'dimless':
        filename = 'APD_fine_dimless.csv'
    else:
        filename = 'APD_fine_litHill.csv'
    APD_trapping = [float('nan') if np.isnan(i).any() else max(i)
                    for i in APD_trapping]
    APD_conductance = [float('nan') if np.isnan(i).any() else max(i)
                       for i in APD_conductance]
if MODE in ["AP", "APD", "APD-qNet"]:
    if MODE == "AP":
        column_name = ['pulse ' + str(i) for i in range(save_signal)]
    else:
        column_name = ['APD']
    APD_trapping_df = pd.DataFrame(APD_trapping, columns=column_name)
    APD_trapping_df['drug concentration'] = drug_conc
    APD_trapping_df.to_csv(os.path.join(data_dir, f'SD_{filename}'))

    APD_conductance_df = pd.DataFrame(APD_conductance, columns=column_name)
    APD_conductance_df['drug concentration'] = drug_conc
    APD_conductance_df.to_csv(os.path.join(data_dir, f'CS_{filename}'))

# Compute qNet
if MODE in ["qNet", "APD-qNet"]:
    pulse_time = 2000
    base_log_key = [APsim.time_key, APsim.Vm_key]
    qNet_current_list = [model_keys['ICaL'], model_keys['INaL'], current_key,
                         model_keys['IKs'], model_keys['IK1'],
                         model_keys['Ito']]
    qNet_current_list = [i for i in qNet_current_list if i is not None]

    APsim.sim.set_protocol(myokit.pacing.blocktrain(period=pulse_time,
                                                    duration=0.5,
                                                    offset=50))
    APsim._cycle_length = pulse_time
    APsim.reset_parameters()
    states = myokit.load_state(
        os.path.join(modelling.RESULT_DIR, 'steady_states',
                     f'{APmodel_name}_steadystate_qNetprot.csv'))
    APsim.set_initial_state(states)

    qNet_SD_arr = []
    qNet_CS_arr = []
    print('Computing qNet...')
    APsim.set_SD_parameters(drug, ikr_model=ikr_model)
    if args.conc == 'dimless':
        SD_params = modelling.BindingParameters(
            ikr_model=ikr_model).load_SD_parameters(drug)
        drug_conc = from_dimless(dimless_input, SD_params['n'].values,
                                 SD_params['halfmax'].values)

    for i in range(len(drug_conc)):
        print(f'simulating concentration: {drug_conc[i]}')

        # Run simulation for the AP-SD model till steady state
        APsim.set_conc(drug_conc[i])
        log = APsim.simulate(save_signal=2, timestep=0.01,
                             log_var=base_log_key + qNet_current_list)

        # Check for EAD-like behaviour
        apd90 = APsim.APD90(log)
        if np.isnan(apd90).any():
            qNet = float('nan')
        else:
            qNet = APsim.qNet(log, period=2)
        # qNet = APsim.qNet(log)
        qNet_SD_arr.append(qNet)

        # Run simulation for the AP-CS model till steady state
        reduction_scale = HillModel.simulate(Hill_coef, drug_conc[i])
        APsim.set_conc(0)
        APsim.set_CS_parameter(reduction_scale)
        log = APsim.simulate(save_signal=2, timestep=0.01,
                             log_var=base_log_key + qNet_current_list)
        APsim.set_CS_parameter(1)

        # Check for EAD-like behaviour
        apd90 = APsim.APD90(log)
        if np.isnan(apd90).any():
            qNet = float('nan')
        else:
            qNet = APsim.qNet(log, period=2)
        # qNet = APsim.qNet(log)
        qNet_CS_arr.append(qNet)

    # Save qNet
    qNet_SD_df = pd.DataFrame(np.array(qNet_SD_arr), columns=['qNet'])
    qNet_SD_df['drug concentration'] = drug_conc
    qNet_CS_df = pd.DataFrame(np.array(qNet_CS_arr), columns=['qNet'])
    qNet_CS_df['drug concentration'] = drug_conc
    if args.conc == 'dimless':
        qNet_SD_df.to_csv(os.path.join(data_dir, 'SD_qNet_EAD_dimless.csv'))
        qNet_CS_df.to_csv(os.path.join(data_dir, 'CS_qNet_EAD_dimless.csv'))
    else:
        qNet_SD_df.to_csv(os.path.join(data_dir, 'SD_qNet_EAD_litHill.csv'))
        qNet_CS_df.to_csv(os.path.join(data_dir, 'CS_qNet_EAD_litHill.csv'))
