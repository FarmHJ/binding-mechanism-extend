#
# Compute the APD90 differences between the AP-SD model and the AP-CS model
# for all synthetic drug.
#

import myokit
import numpy as np
import os
import pandas as pd
import sys

import modelling

APmodel_name = sys.argv[1]

# Define directories to save simulation data
data_dir = '../simulation_data/parameter_SA/'

# Load IKr model and set up protocol
model_filepath = '../math_model/current_model/ohara-cipa-2017-IKr.mmt'
model, _, x = myokit.load(model_filepath)
Milnes_protocol = modelling.ProtocolLibrary().Milnes(25e3)
current_model = modelling.Simulation(model, protocol=Milnes_protocol,
                                     current_head_key='ikr')

# Load AP model and set up protocol
# if APmodel_name == 'Grandi':
#     AP_model_filepath = '../math_model/AP_model/Grd-2010-IKr-SD.mmt'
# elif APmodel_name == 'TTP':
#     AP_model_filepath = '../math_model/AP_model/TTP-2006-IKr-SD.mmt'
model_details = modelling.ModelDetails()
AP_model_filepath = '../' + model_details.file_names[
    APmodel_name]['AP_SD_path']
APmodel, _, x = myokit.load(AP_model_filepath)
model_keys = model_details.current_keys[APmodel_name]
current_key = model_keys['IKr']
time_key = model_keys['time']
Vm_key = model_keys['Vm']
current_head_key = current_key[:current_key.index('.')]
AP_model = modelling.Simulation(APmodel, current_head_key=current_head_key)

pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)

# Define constants for simulations
offset = 50
save_signal = 2
repeats = 1000
APD_points = 20

# Get list of synthetic drug
param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

# Get name of parameters
param_names = modelling.SDModelDetails().param_names


def param_evaluation(param_values, drug):

    # Define parameter values of synthetic drug
    print('Running for drug: ', drug)
    orig_half_effect_conc = param_values['EC50'][0]
    param_values.loc[0, 'EC50'] = 1
    ComparisonController.drug_param_values = param_values

    # Calculate the normalising constant
    Hill_n = param_values['N'][0]
    norm_constant = np.power(orig_half_effect_conc, 1 / Hill_n)

    # Compute Hill curve of the synthetic drug with the SD model
    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        ComparisonController.compute_Hill(current_model,
                                          norm_constant=norm_constant,
                                          parallel=False)
    # The parameters of Hill curve are based on the normalised drug
    # concentration.
    # Hill coefficient remains the same but IC50 -> IC50/EC50

    # Define drug concentration range similar to the drug concentration used
    # to infer Hill curve
    drug_conc_AP = 10**np.linspace(np.log10(drug_conc_Hill[1]),
                                   np.log10(max(drug_conc_Hill)),
                                   APD_points)

    # Simulate APs and APD90s of the AP-SD model and the AP-CS model
    APD_trapping, APD_conductance, drug_conc_AP = \
        ComparisonController.APD_sim(
            AP_model, Hill_curve_coefs, drug_conc=drug_conc_AP,
            IKr_tuning_factor=scaling_factor, EAD=True)

    # Calculate RMSD and MD of simulated APD90 of the two models
    RMSError = ComparisonController.compute_RMSE(APD_trapping,
                                                 APD_conductance)
    MAError = ComparisonController.compute_ME(APD_trapping,
                                              APD_conductance)

    # Create dataframe to save results
    conc_Hill_ind = ['conc_' + str(i) for i, _ in
                     enumerate(drug_conc_Hill)]
    conc_AP_ind = ['conc_' + str(i) for i, _ in enumerate(drug_conc_AP)]
    index_dict = {'drug': ['drug'],
                  'drug_conc_Hill': conc_Hill_ind,
                  'peak_current': conc_Hill_ind,
                  'Hill_curve': ['Hill_coef', 'IC50'],
                  'param_values': param_names, 'drug_conc_AP': conc_AP_ind,
                  'APD_trapping': conc_AP_ind,
                  'APD_conductance': conc_AP_ind, 'RMSE': ['RMSE'],
                  'ME': ['ME']}
    all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
    index = pd.MultiIndex.from_tuples(all_index)

    param_values.loc[0, 'EC50'] = orig_half_effect_conc
    big_df = pd.DataFrame(
        [drug] + list(drug_conc_Hill) + list(peaks_norm) +
        list(Hill_curve_coefs) + list(param_values.values[0]) +
        list(drug_conc_AP) + APD_trapping + APD_conductance + [RMSError] +
        [MAError], index=index)

    return big_df


# Load IKr tuning factor
scaling_df_filepath = '../simulation_data/' + APmodel_name + \
    '_conductance_scale.csv'
scaling_df = pd.read_csv(scaling_df_filepath, index_col=[0])
if APmodel_name == 'Lei':
    scaling_factor = scaling_df.loc['AP_duration']['conductance scale']
else:
    scaling_factor = scaling_df.loc['hERG_peak']['conductance scale']

# Determine completed simulations so that same simulations are not repeated
filename = 'SA_alldrugs_' + APmodel_name + '_tuned.csv'
if os.path.exists(data_dir + filename):
    results_df = pd.read_csv(data_dir + filename, header=[0, 1], index_col=[0],
                             skipinitialspace=True)
    ran_drugs = results_df['drug']['drug'].values
else:
    ran_drugs = []

drug_list = [i for i in drug_list if i not in ran_drugs]

first_iter = True
for drug in drug_list:

    # Get parameter values of each synthetic drug
    Vhalf = param_lib.binding_parameters[drug]['Vhalf']
    Kmax = param_lib.binding_parameters[drug]['Kmax']
    Ku = param_lib.binding_parameters[drug]['Ku']
    Hill_n = param_lib.binding_parameters[drug]['N']
    half_effect_conc = param_lib.binding_parameters[drug]['EC50']

    all_params = [Vhalf, Kmax, Ku, Hill_n, half_effect_conc]

    # Define parameter values input to the system
    orig_param_values = pd.DataFrame(all_params, index=param_names)
    orig_param_values = orig_param_values.T
    ComparisonController = modelling.ModelComparison(orig_param_values,
                                                     [time_key, Vm_key])

    # Evaluate the RMSD and MD between APD90s of a synthetic drug from the
    # AP-SD model and the AP-CS model
    if not os.path.exists(data_dir + filename):
        results_df = param_evaluation(orig_param_values, drug)
        results_df = results_df.T
        first_iter = False
    else:
        results_df = pd.read_csv(data_dir + filename, header=[0, 1],
                                 index_col=[0], skipinitialspace=True)
        big_df = param_evaluation(orig_param_values, drug)
        results_df = pd.concat([results_df, big_df.T])

    results_df.to_csv(data_dir + filename)
