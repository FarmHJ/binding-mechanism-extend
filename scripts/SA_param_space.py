# Explore the parameter space of drug-related parameters (Vhalf, Kmax and Ku).
# Compute the APD90 differences between the AP-SD model and the AP-CS model
# for a given virtual drug.

import argparse
import glob
import itertools
import myokit
import numpy as np
import os
import pandas as pd
import pints
import time

import modelling

parser = argparse.ArgumentParser(
    description="Parameter space exploration for AP-IKr-SD model and "
    "the AP-IKr-CS model")
parser.add_argument("APmodel", help="Name of AP model")
args = parser.parse_args()

APmodel = args.APmodel

# Define directory to save simulation data
data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration')
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

APsim = modelling.ModelSimController(APmodel)
states = myokit.load_state(
    os.path.join(modelling.RESULT_DIR, 'steady_states',
                 f'{APmodel}_steadystate_APDprot.csv'))
APsim.set_initial_state(states)

if APmodel != 'ORd-Li':
    APsim.set_ikr_rescale_method('AP_duration')

# Get name of parameters
param_names = modelling.SD_details.SD_param_names
if APmodel == 'ORd-Lei':
    IKrmodel = 'Lei'
    IKr_sim = modelling.ModelSimController(IKrmodel)
    states = myokit.load_state(
        os.path.join(modelling.RESULT_DIR, 'steady_states',
                     f'{IKrmodel}_steadystate_Milnes.csv'))
    IKr_sim.set_initial_state(states)
    ComparisonController = modelling.ModelComparison(APsim,
                                                     IKr_simulator=IKr_sim)
else:
    IKrmodel = 'Li'
    ComparisonController = modelling.ModelComparison(APsim)


def param_evaluation(inputs, skip_ikr):

    # Prepare the inputs for simulation
    ComparisonController.prepare_inputs(inputs, skip_ikr=skip_ikr)
    if not skip_ikr:
        start = time.time()
        ComparisonController.get_drug_effect(parallel=False)
        print('Get Hill curve: ', time.time() - start)

    try:    
        # Simulate APs and APD90s of the AP-SD model and the AP-CS model
        start = time.time()
        ComparisonController.get_APD()
        print('Get APD: ', time.time() - start)

        # Calculate RMSD and MD of simulated APD90 of the two models
        ComparisonController.RMSE()
        ComparisonController.ME()

    except myokit.SimulationError:
        ComparisonController.APD_trapping = [float("Nan")] * \
            len(ComparisonController.drug_conc_AP)
        ComparisonController.APD_conductance = [float("Nan")] * \
            len(ComparisonController.drug_conc_AP)
        ComparisonController.RMSError = float("Nan")
        ComparisonController.MAError = float("Nan")

    outcome_df = ComparisonController.process_data()

    return outcome_df


# Assuming drug concentration are all normalised, the EC50 value in the model
# becomes 1.
# Since Hill coefficient, N, does not affect APD difference behaviour, it
# can be fixed at any value.
# For simplicity, let N = 1.

# Save defined parameter space or load previously saved parameter space
param_space_dir = os.path.join(data_dir, 'parameter_space')
if not os.path.isdir(param_space_dir):
    os.makedirs(param_space_dir)

param_space_fpath = os.path.join(param_space_dir,
                                 f'parameter_space_{IKrmodel}.csv')

param_space = []
if os.path.exists(param_space_fpath):
    param_dict = pd.read_csv(param_space_fpath, header=[0, 1], index_col=[0],
                             skipinitialspace=True)
    param_dict = param_dict.rename(columns={"N": "n", "EC50": "halfmax"})
    param_dict = param_dict.to_dict(orient='list')
else:
    if APmodel == 'ORd-Lei':
        Vhalf_fullrange = np.linspace(-200, 0, 20)
        Kmax_fullrange = 10**np.linspace(1, 10, 20)
        Ku_fullrange = 10**np.linspace(-5, 1, 20)

    counter = 0
    col_header = ['param_id'] + ['param_values'] * 5
    dict_keys = ['param_id'] + param_names
    param_dict = {(i, j): [] for (i, j) in zip(col_header, dict_keys)}
    for Vhalf, Kmax, Ku in itertools.product(
            Vhalf_fullrange, Kmax_fullrange, Ku_fullrange):

        param_values = [counter, Kmax, Ku, 1, 1, Vhalf]
        for i in range(len(dict_keys)):
            param_dict[(col_header[i], dict_keys[i])].append(param_values[i])

        counter += 1

    param_df = pd.DataFrame.from_dict(param_dict)
    param_df = param_df.set_axis(
        pd.MultiIndex.from_arrays([col_header, dict_keys]), axis=1)
    param_df.to_csv(param_space_fpath)
    del param_df

# Set up variables for data saving
id_key = ('param_id', 'param_id')
total_samples = len(param_dict[id_key])
samples_per_save = 1000
samples_split_n = int(np.ceil(total_samples / samples_per_save))
total_saving_file_num = np.arange(samples_split_n)

# Determine completed simulations so that it is not repeated
file_prefix = 'SA_allparam_'
results_dir = os.path.join(data_dir, 'SA_space', APmodel)
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)
existing_result_files = glob.glob(os.path.join(results_dir,
                                               f'{file_prefix}*.csv'))
if len(existing_result_files) == 0:
    file_id_dict = {}
    for i in range(samples_split_n):
        file_id_dict[i] = param_dict[id_key][
            i * samples_per_save: (i + 1) * samples_per_save]
    saving_file_dict = {'file_num': total_saving_file_num,
                        'sample_id_each_file': file_id_dict}
else:
    # Check for missing results file
    result_files_num = [int(os.path.basename(fname)[len(file_prefix):-4])
                        for fname in existing_result_files]
    file_num_to_run = []
    file_id_dict = {}
    missing_file = [i for i in total_saving_file_num
                    if i not in result_files_num]
    for i in missing_file:
        file_id_dict[i] = param_dict[id_key][
            i * samples_per_save: (i + 1) * samples_per_save]
        file_num_to_run.append(i)
    for file in existing_result_files:
        file_num = int(os.path.basename(file)[len(file_prefix):-4])
        saved_results_df = pd.read_csv(os.path.join(data_dir, file),
                                       header=[0, 1], index_col=[0],
                                       skipinitialspace=True)
        ran_values = saved_results_df.index.values
        expected_ids = param_dict[id_key][
            file_num * samples_per_save: (file_num + 1) * samples_per_save]
        param_space_id = [i for i in expected_ids if i not in ran_values]
        if len(param_space_id) != 0:
            file_num_to_run.append(file_num)
            file_id_dict[file_num] = param_space_id
    saving_file_dict = {'file_num': sorted(file_num_to_run),
                        'sample_id_each_file': file_id_dict}

# Use PINTS' parallel evaluator to evaluate the APD90 difference for each
# virtual drugs in the parameter space
n_workers = 6
skip_ikr = False if APmodel == 'ORd-Lei' else True
evaluator = pints.ParallelEvaluator(param_evaluation,
                                    n_workers=n_workers,
                                    args=[skip_ikr])
for file_num in saving_file_dict['file_num']:
    print('Starting function evaluation for file number: ', file_num)
    current_time = time.strftime("%H:%M:%S", time.localtime())
    print('Starting time: ', current_time)

    samples_to_run = saving_file_dict['sample_id_each_file'][file_num]
    samples_num = len(samples_to_run)
    filename = f'{file_prefix}{file_num}.csv'
    filepath = os.path.join(results_dir, filename)

    for i in range(int(np.ceil(samples_num / n_workers))):
        subset_samples_to_run = samples_to_run[
            n_workers * i:n_workers * (i + 1)]
        start = time.time()
        print('Running samples ', int(subset_samples_to_run[0]), ' to ',
              int(subset_samples_to_run[-1]))
        param_space = []
        for j in subset_samples_to_run:
            ind = param_dict[id_key].index(j)
            param_values = {n: param_dict[('param_values', n)][ind]
                            for n in param_names}
            input = {'id': j, 'param_values': param_values}
            if APmodel != 'ORd-Lei':
                key_list = [k for k in param_dict.keys()
                            if 'drug_conc_Hill' in k]
                drug_conc = [param_dict[k][ind] for k in key_list]
                key_list = [k for k in param_dict.keys() if 'Hill_curve' in k]
                Hill_curve = {k[1]: param_dict[k][ind] for k in key_list}
                input.update({'drug_conc_Hill': drug_conc,
                              'Hill_curve': Hill_curve})
            param_space.append(input)

        big_df = evaluator.evaluate(param_space)
        print('Time taken in minutes: ', (time.time() - start) / 60)

        if os.path.exists(filepath):
            combined_df = pd.read_csv(filepath, header=[0, 1], index_col=[0],
                                      skipinitialspace=True)
            for i in range(len(big_df)):
                combined_df = pd.concat([combined_df, big_df[i]])
        else:
            combined_df = big_df[0]
            for i in range(1, len(big_df)):
                combined_df = pd.concat([combined_df, big_df[i]])

        combined_df.to_csv(filepath)
