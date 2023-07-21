#
# Compute the RMSD between APD90s of the AP-SD model and the AP-CS model for
# each synthetic drug with varying parameter N.
#

import myokit
import numpy as np
import os
import pandas as pd
import pints
import sys

import modelling

APmodel_name = sys.argv[1]

# Define directory to save simulation data
data_filepath = '../simulation_data/parameter_SA/APD90diff_N/' + \
    APmodel_name + '/'
if not os.path.isdir(data_filepath):
    os.makedirs(data_filepath)

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
model_keys = modelling.ModelDetails().current_keys[APmodel_name]
current_key = model_keys['IKr']
time_key = model_keys['time']
Vm_key = model_keys['Vm']
current_head_key = current_key[:current_key.index('.')]
AP_model = modelling.Simulation(APmodel, current_head_key=current_head_key)

pulse_time = 1000
AP_model.protocol = modelling.ProtocolLibrary().current_impulse(pulse_time)
base_conductance = APmodel.get(current_head_key + '.gKr').value()

# Define parameters used in simulations
offset = 50
save_signal = 2
repeats = 1000
APD_points = 20

# Getn list of synthetic drugs and the name of parameters
param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds
param_names = modelling.SDModelDetails().param_names
parameter_interest = 'N'

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df,
                                                 [time_key, Vm_key])


def param_evaluation(param, param_values):

    # Define parameter values of virtual drug
    param_values.loc[0, parameter_interest] = param
    orig_half_effect_conc = param_values['EC50'][0]
    param_values.loc[0, 'EC50'] = 1
    ComparisonController.drug_param_values = param_values

    Hill_n = param_values['N'][0]
    norm_constant = np.power(orig_half_effect_conc, 1 / Hill_n)

    # Compute Hill curve of the synthetic drug with the SD model
    Hill_curve_coefs, drug_conc_Hill, peaks_norm = \
        ComparisonController.compute_Hill(current_model,
                                          norm_constant=norm_constant,
                                          max_counter=30,
                                          parallel=False)

    # Define drug concentration range similar to the drug concentration used
    # to infer Hill curve
    drug_conc_AP = 10**np.linspace(np.log10(drug_conc_Hill[1]),
                                   np.log10(max(drug_conc_Hill)),
                                   APD_points)
    try:
        # Simulate APs and APD90s of the ORd-SD model and the ORd-CS model
        APD_trapping, APD_conductance, drug_conc_AP = \
            ComparisonController.APD_sim(
                AP_model, Hill_curve_coefs, drug_conc=drug_conc_AP,
                IKr_tuning_factor=scaling_factor, EAD=True,)
                # abs_tol=1e-8, rel_tol=1e-9)

        # Calculate RMSD and MD of simulated APD90 of the two models
        RMSError = ComparisonController.compute_RMSE(APD_trapping,
                                                     APD_conductance)
        MAError = ComparisonController.compute_ME(APD_trapping,
                                                  APD_conductance)
    except myokit.SimulationError:
        APD_trapping = [float("Nan")] * APD_points
        APD_conductance = [float("Nan")] * APD_points
        RMSError = float("Nan")
        MAError = float("Nan")

    # Create dataframe to save results
    conc_Hill_ind = ['conc_' + str(i) for i, _ in
                     enumerate(drug_conc_Hill)]
    conc_AP_ind = ['conc_' + str(i) for i, _ in enumerate(drug_conc_AP)]
    index_dict = {'drug_conc_Hill': conc_Hill_ind,
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
        list(drug_conc_Hill) + list(peaks_norm) + list(Hill_curve_coefs) +
        list(param_values.values[0]) + list(drug_conc_AP) + APD_trapping +
        APD_conductance + [RMSError] + [MAError], index=index)

    return big_df


scaling_df_filepath = '../simulation_data/' + APmodel_name + \
    '_conductance_scale.csv'
scaling_df = pd.read_csv(scaling_df_filepath, index_col=[0])
if APmodel_name == 'Lei':
    scaling_factor = scaling_df.loc['AP_duration']['conductance scale']
else:
    scaling_factor = scaling_df.loc['hERG_peak']['conductance scale']

interest_param_list = []
for i in drug_list:
    interest_param_list.append(
        param_lib.binding_parameters[i][parameter_interest])

min_value = min(interest_param_list)
max_value = max(interest_param_list)
param_fullrange = np.linspace(min_value, max_value, 20)

# Load previously saved parameter space
data_mainpath = '../simulation_data/parameter_SA/'
sample_filepath = data_mainpath + 'SA_alldrugs.csv'

drug_space = []
# if os.path.exists(sample_filepath):
drug_space_df = pd.read_csv(sample_filepath,
                            header=[0, 1], index_col=[0],
                            skipinitialspace=True)

# drug_list = [drug_list[-2]]

for drug in drug_list:
    print(drug)
    # Get parameter values of each synthetic drug
    param_values = drug_space_df.loc[drug_space_df[('drug', 'drug')] == drug][
        'param_values']

    # Define parameter values to the system
    ComparisonController = modelling.ModelComparison(
        param_values, [time_key, Vm_key])

    # Check for completed simulations to prevent repetition
    filename = 'SA_' + drug + '_' + parameter_interest + '.csv'
    if os.path.exists(data_filepath + filename):
        saved_results_df = pd.read_csv(data_filepath + filename,
                                       header=[0, 1], index_col=[0],
                                       skipinitialspace=True)
        ran_values = saved_results_df['param_values'][
            parameter_interest].values
    else:
        ran_values = []

    param_range = [i for i in param_fullrange if i not in ran_values]
    print(param_range)

    # Evaluate the RMSD and MD between APD90s of a synthetic drug with
    # changing Hill coefficient from the ORd-SD model and the ORd-CS model
    n_workers = 8
    evaluator = pints.ParallelEvaluator(param_evaluation,
                                        n_workers=n_workers,
                                        args=[param_values])
    for i in range(int(np.ceil(len(param_range) / n_workers))):
        print('Running samples ', n_workers * i, 'to',
              n_workers * (i + 1) - 1)
        big_df = evaluator.evaluate(
            param_range[i * n_workers: (i + 1) * n_workers])

        if os.path.exists(data_filepath + filename):
            combined_df = pd.read_csv(data_filepath + filename,
                                      header=[0, 1], index_col=[0],
                                      skipinitialspace=True)
            for i in range(len(big_df)):
                combined_df = pd.concat([combined_df, big_df[i].T])
        else:
            combined_df = big_df[0].T
            for i in range(1, len(big_df)):
                combined_df = pd.concat([combined_df, big_df[i].T])

        combined_df.to_csv(data_filepath + filename)
