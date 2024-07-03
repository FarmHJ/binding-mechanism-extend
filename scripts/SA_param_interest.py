#
# Compute the RMSD between APD90s of the AP-SD model and the AP-CS model for
# each synthetic drug with varying parameter of interest.
#

import argparse
import myokit
import numpy as np
import os
import pandas as pd
import pints

import modelling

parser = argparse.ArgumentParser(
    description="Parameter space exploration for AP-IKr-SD model and "
    "the AP-IKr-CS model")
parser.add_argument("APmodel", help="Name of AP model")
args = parser.parse_args()

APmodel = args.APmodel

# Define directory to save simulation data
data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA')
if args.mode == 'parameter_SA':
    data_dir = os.path.join(data_dir, 'parameter_n', APmodel)
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

# Get list of synthetic drugs and the name of parameters
drug_list = modelling.SD_details.drug_names

APsim = modelling.ModelSimController(APmodel)
states = myokit.load_state(
    os.path.join(modelling.RESULT_DIR, 'steady_states',
                 f'{APmodel}_steadystate_APDprot.csv'))
APsim.set_initial_state(states)
APsim.set_ikr_rescale_method('AP_duration')

# Get name of parameters
param_names = modelling.SD_details.SD_param_names
if APmodel == 'ORd-Lei':
    ikr_model = 'Lei'
    IKr_sim = modelling.ModelSimController(ikr_model)
    states = myokit.load_state(
        os.path.join(modelling.RESULT_DIR, 'steady_states',
                     f'{ikr_model}_steadystate_Milnes.csv'))
    IKr_sim.set_initial_state(states)
    ComparisonController = modelling.ModelComparison(APsim,
                                                     IKr_simulator=IKr_sim)
else:
    ikr_model = 'Li'
    ComparisonController = modelling.ModelComparison(APsim)


def param_evaluation(inputs, counter):

    # Prepare the inputs for simulation
    ComparisonController.prepare_inputs(inputs)

    # Set the normalising constant
    ComparisonController.normalise_drug_conc()

    # Compute Hill curve of the synthetic drug with the SD model
    ComparisonController.get_drug_effect(max_counter=counter,
                                         parallel=False)
    try:
        # Simulate APs and APD90s of the ORd-SD model and the ORd-CS model
        ComparisonController.get_APD()

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

    outcome_df = ComparisonController.process_data(save_orig_halfmax=True)

    return outcome_df


SD_params = modelling.BindingParameters(ikr_model=ikr_model)
SD_params.load_SD_parameters()

########################################
# Run model comparison for 12 CiPA drugs
filename = f'SA_alldrugs_{APmodel}.csv'
fpath = os.path.join(data_dir, filename)

# Load completed simulation to reduce repetition
if os.path.exists(fpath):
    results_df = pd.read_csv(fpath, header=[0, 1], index_col=[0],
                             skipinitialspace=True)
    ran_drugs = results_df.index
else:
    ran_drugs = []
drug_list = [i for i in drug_list if i not in ran_drugs]

# Combine parameter combinations of 12 CiPA drugs
param_space = []
for drug in drug_list:
    param_values = SD_params.get_SD_parameters(drug)
    input = {'id': drug, 'param_values': param_values}
    param_space.append(input)

# Evaluate the RMSD and MD for 12 CiPA drugs
n_workers = 8
evaluator = pints.ParallelEvaluator(param_evaluation,
                                    n_workers=n_workers,
                                    args=[20])
for i in range(int(np.ceil(len(param_space) / n_workers))):
    print('Running samples ', n_workers * i, 'to',
          n_workers * (i + 1) - 1)
    sub_param_space = param_space[i * n_workers: (i + 1) * n_workers]

    big_df = evaluator.evaluate(sub_param_space)

    if os.path.exists(fpath):
        combined_df = pd.read_csv(fpath, header=[0, 1], index_col=[0],
                                  skipinitialspace=True)
        for i in range(len(big_df)):
            combined_df = pd.concat([combined_df, big_df[i]])
    else:
        combined_df = big_df[0]
        for i in range(1, len(big_df)):
            combined_df = pd.concat([combined_df, big_df[i]])

    combined_df.to_csv(fpath)

###############################################
# Run model comparison with varying parameter n
# Take min and max of parameter n from the 12 CiPA drugs
param_value_list = [SD_params.get_SD_parameters(i)['n'].values[0]
                    for i in drug_list]
param_fullrange = np.linspace(min(param_value_list),
                              max(param_value_list), 20)

# Load previously saved parameter space
drug_space_df = pd.read_csv(
    os.path.join(modelling.RESULT_DIR, 'parameter_SA',
                 f'SA_alldrugs_{APmodel}.csv'),
    header=[0, 1], index_col=[0], skipinitialspace=True)

for drug in drug_list:
    print('parameter SA: ', drug)

    # Check for completed simulations to prevent repetition
    filename = f'SA_{drug}_n.csv'
    filepath = os.path.join(data_dir, filename)
    if os.path.exists(filepath):
        saved_results_df = pd.read_csv(filepath, header=[0, 1],
                                       index_col=[0],
                                       skipinitialspace=True)
        ran_values = saved_results_df['param_values']['n'].values
    else:
        ran_values = []

    param_range = [i for i in param_fullrange if i not in ran_values]
    param_space = []
    for v in param_range:
        # Get parameter values of each synthetic drug
        param_values = drug_space_df.loc[[drug]]['param_values']
        param_values['n'] = v
        input = {'id': drug, 'param_values': param_values}
        param_space.append(input)

    # Evaluate the RMSD and MD between APD90s of a synthetic drug with
    # changing Hill coefficient
    n_workers = 8
    evaluator = pints.ParallelEvaluator(param_evaluation,
                                        n_workers=n_workers,
                                        args=[30])
    for i in range(int(np.ceil(len(param_space) / n_workers))):
        print('Running samples ', n_workers * i, 'to', n_workers * (i + 1) - 1)
        sub_param_space = param_space[i * n_workers: (i + 1) * n_workers]

        # Evaluate RMSD and MD for a subset of the parameter combination list
        big_df = evaluator.evaluate(sub_param_space)

        # Save results
        if os.path.exists(filepath):
            combined_df = pd.read_csv(filepath, header=[0, 1],
                                      index_col=[0], skipinitialspace=True)
            for i in range(len(big_df)):
                combined_df = pd.concat([combined_df, big_df[i]])
        else:
            combined_df = big_df[0]
            for i in range(1, len(big_df)):
                combined_df = pd.concat([combined_df, big_df[i]])

        combined_df.to_csv(filepath)
