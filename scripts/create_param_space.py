# Explore the parameter space of drug-related parameters (Vhalf, Kmax and Ku).
# Compute the APD90 differences between the AP-SD model and the AP-CS model
# for a given virtual drug.

import argparse
import glob
import itertools
import numpy as np
import os
import pandas as pd

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

# Get name of parameters
param_names = modelling.SD_details.SD_param_names
IKrmodel = 'Lei' if APmodel == 'ORd-Lei' else 'Li'

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
    missing_ids = []
    for i in range(samples_split_n):
        missing_ids += param_dict[id_key][
            i * samples_per_save: (i + 1) * samples_per_save]
else:
    # Check for missing results file
    result_files_num = [int(os.path.basename(fname)[len(file_prefix):-4])
                        for fname in existing_result_files]
    missing_file = [i for i in total_saving_file_num
                    if i not in result_files_num]
    missing_ids = []
    for i in missing_file:
        missing_ids += param_dict[id_key][
            i * samples_per_save: (i + 1) * samples_per_save]
    for file in existing_result_files:
        file_num = int(os.path.basename(file)[len(file_prefix):-4])
        saved_results_df = pd.read_csv(os.path.join(data_dir, file),
                                       header=[0, 1], index_col=[0],
                                       skipinitialspace=True)
        ran_values = saved_results_df.index.values
        expected_ids = param_dict[id_key][
            file_num * samples_per_save: (file_num + 1) * samples_per_save]
        param_space_id = [i for i in expected_ids if i not in ran_values]
        missing_ids += param_space_id

missing_ids = sorted(missing_ids)
fpath = os.path.join(data_dir, 'SA_space', APmodel, 'missing_id.txt')
with open(fpath, 'w') as fp:
    for i in missing_ids:
        fp.write("%d\n" % i)
