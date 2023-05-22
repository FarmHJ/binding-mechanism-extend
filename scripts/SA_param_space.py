#
# Explore the parameter space of drug-related parameters (Vhalf, Kmax and Ku).
# Compute the APD90 differences between the AP-SD model and the AP-CS model
# for a given virtual drug.
#

import myokit
import numpy as np
import os
import pandas as pd
import pints
import sys
import time

import modelling

# Define directory to save simulation data
data_filepath = '../simulation_data/parameter_space_exploration/'
if not os.path.exists(data_filepath):
    os.makedirs(data_filepath)

APmodel_name = sys.argv[1]
model_details = modelling.ModelDetails()

# Model directory
AP_model_filepath = '../' + model_details.file_names[
    APmodel_name]['AP_SD_path']

# Load AP model and set current protocol
APmodel, _, x = myokit.load(AP_model_filepath)
model_keys = model_details.current_keys[APmodel_name]
current_key = model_keys['IKr']
time_key = model_keys['time']
Vm_key = model_keys['Vm']
current_head_key = current_key[:current_key.index('.')]
AP_model = modelling.Simulation(APmodel, current_head_key=current_head_key)

pulse_time = 1000
protocol_offset = 50
protocol = myokit.pacing.blocktrain(pulse_time, 0.5, offset=protocol_offset)
AP_model.protocol = protocol

# Define constants for simulations
offset = 50
save_signal = 2
repeats = 1000
APD_points = 20

# Get name of parameters
param_names = modelling.SDModelDetails().param_names

starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
ComparisonController = modelling.ModelComparison(starting_param_df,
                                                 [time_key, Vm_key])


def param_evaluation(inputs):

    # Define parameter values of virtual drug
    param_id = inputs[('param_id', 'param_id')][0]
    param_values = inputs['param_values']
    drug_conc_Hill = inputs['drug_conc_Hill'].values[0]
    drug_conc_Hill = drug_conc_Hill[~np.isnan(drug_conc_Hill)]

    ComparisonController.drug_param_values = param_values

    Hill_curve_coefs = inputs['Hill_curve'].values[0]

    # Define drug concentration range similar to the drug concentration used
    # to infer Hill curve
    drug_conc_AP = 10**np.linspace(np.log10(drug_conc_Hill[1]),
                                   np.log10(max(drug_conc_Hill)),
                                   APD_points)

    try:
        # Simulate APs and APD90s of the AP-SD model and the AP-CS model
        print('running simulation')
        APD_trapping, APD_conductance, drug_conc_AP = \
            ComparisonController.APD_sim(
                AP_model, Hill_curve_coefs, drug_conc=drug_conc_AP,
                IKr_tuning_factor=scaling_factor, EAD=True)
        print('simulation done')

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
        print('simulation error')

    # # Create dataframe to save results
    conc_Hill_ind = ['conc_' + str(i) for i, _ in
                     enumerate(drug_conc_Hill)]
    conc_AP_ind = ['conc_' + str(i) for i, _ in enumerate(drug_conc_AP)]
    index_dict = {'param_id': ['param_id'],
                  'drug_conc_Hill': conc_Hill_ind,
                  'Hill_curve': ['Hill_coef', 'IC50'],
                  'param_values': param_names, 'drug_conc_AP': conc_AP_ind,
                  'APD_trapping': conc_AP_ind,
                  'APD_conductance': conc_AP_ind, 'RMSE': ['RMSE'],
                  'ME': ['ME']}
    all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
    index = pd.MultiIndex.from_tuples(all_index)

    big_df = pd.DataFrame(
        [param_id] + list(drug_conc_Hill) + list(Hill_curve_coefs) +
        list(param_values.values[0]) + list(drug_conc_AP) +
        APD_trapping + APD_conductance + [RMSError] + [MAError],
        index=index)

    return big_df


# Load IKr tuning factor
scaling_df_filepath = '../simulation_data/' + APmodel_name + \
    '_conductance_scale.csv'
scaling_df = pd.read_csv(scaling_df_filepath, index_col=[0])
if APmodel_name == 'Lei':
    scaling_factor = scaling_df.loc['AP_duration']['conductance scale']
else:
    scaling_factor = scaling_df.loc['hERG_peak']['conductance scale']

# Assuming drug concentration are all normalised, the EC50 value in the model
# becomes 1.
# Since Hill coefficient, N, does not affect APD difference behaviour, it
# can be fixed at any value.
# For simplicity, let N = 1.

# Save defined parameter space or load previously saved parameter space
sample_filepath = data_filepath + 'Hill_curves.csv'

param_space = []
# if os.path.exists(sample_filepath):
param_values_df = pd.read_csv(sample_filepath,
                              header=[0, 1], index_col=[0],
                              skipinitialspace=True)
for i in range(len(param_values_df.index)):
    param_space.append(param_values_df.iloc[[i]])

# Set up variables for data saving
total_samples = len(param_space)
samples_per_save = 1000
samples_split_n = int(np.ceil(total_samples / samples_per_save))
total_saving_file_num = np.arange(samples_split_n)

# Determine completed simulations so that it is not repeated
file_prefix = 'SA_allparam_'
data_dir = data_filepath + 'SA_space/' + APmodel_name + '/'
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
evaluation_result_files = [f for f in os.listdir(data_dir) if
                           f.startswith(file_prefix)]
if len(evaluation_result_files) == 0:
    file_id_dict = {}
    for i in range(samples_split_n):
        file_id_dict[i] = param_values_df[('param_id', 'param_id')].values[
            i * samples_per_save: (i + 1) * samples_per_save]
    saving_file_dict = {'file_num': total_saving_file_num,
                        'sample_id_each_file': file_id_dict}
else:
    # Check for missing results file
    result_files_num = [int(fname[len(file_prefix):-4]) for fname in
                        evaluation_result_files]
    file_num_to_run = []
    file_id_dict = {}
    missing_file = [i for i in total_saving_file_num
                    if i not in result_files_num]
    for i in missing_file:
        file_id_dict[i] = param_values_df[('param_id', 'param_id')].values[
            i * samples_per_save: (i + 1) * samples_per_save]
        file_num_to_run.append(i)
    for file in evaluation_result_files:
        file_num = int(file[len(file_prefix):-4])
        saved_results_df = pd.read_csv(data_dir + file,
                                       header=[0, 1], index_col=[0],
                                       skipinitialspace=True)
        ran_values = saved_results_df['param_id']['param_id'].values
        expected_ids = param_values_df['param_id']['param_id'].values[
            file_num * samples_per_save: (file_num + 1) * samples_per_save]
        param_space_id = [i for i in expected_ids if i not in ran_values]
        if len(param_space) != 0:
            file_num_to_run.append(file_num)
            file_id_dict[file_num] = param_space_id
    saving_file_dict = {'file_num': sorted(file_num_to_run),
                        'sample_id_each_file': file_id_dict}

# Use PINTS' parallel evaluator to evaluate the APD90 difference for each
# virtual drugs in the parameter space
n_workers = 8
evaluator = pints.ParallelEvaluator(param_evaluation,
                                    n_workers=n_workers)
for file_num in saving_file_dict['file_num']:
    print('Starting function evaluation for file number: ', file_num)
    current_time = time.strftime("%H:%M:%S", time.localtime())
    print('Starting time: ', current_time)
    samples_to_run = saving_file_dict['sample_id_each_file'][file_num]
    samples_num = len(samples_to_run)
    filename = file_prefix + str(file_num) + '.csv'

    for i in range(int(np.ceil(samples_num / n_workers))):
        subset_samples_to_run = samples_to_run[
            n_workers * i:n_workers * (i + 1)]
        print('Running samples ', int(subset_samples_to_run[0]), ' to ',
              int(subset_samples_to_run[-1]))
        subset_param_space = param_values_df.loc[
            param_values_df[('param_id', 'param_id')].isin(
                subset_samples_to_run)]
        param_space = []
        for i in range(len(subset_param_space.index)):
            param_space.append(subset_param_space.iloc[[i]])

        big_df = evaluator.evaluate(param_space)

        if os.path.exists(data_dir + filename):
            combined_df = pd.read_csv(data_dir + filename,
                                      header=[0, 1], index_col=[0],
                                      skipinitialspace=True)
            for i in range(len(big_df)):
                combined_df = pd.concat([combined_df, big_df[i].T])
        else:
            combined_df = big_df[0].T
            for i in range(1, len(big_df)):
                combined_df = pd.concat([combined_df, big_df[i].T])

        combined_df.to_csv(data_dir + filename)
