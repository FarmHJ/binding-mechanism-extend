# Explore the parameter space of drug-related parameters (Vhalf, Kmax and Ku).
# Compute the APD90 differences between the AP-SD model and the AP-CS model
# for a given virtual drug.

import argparse
import myokit
import os
import pandas as pd
import time

import modelling

parser = argparse.ArgumentParser(
    description="Parameter space exploration for AP-IKr-SD model and "
    "the AP-IKr-CS model")
parser.add_argument("APmodel", help="Name of AP model")
parser.add_argument('id_num', help='ID of parameter values')
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

# Set up AP models and their corresponding IKr model
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


# Define function to evaluate the RMSD between the AP-SD model and
# the AP-CS model
def param_evaluation(inputs, skip_ikr):

    # Prepare the inputs for simulation
    ComparisonController.prepare_inputs(inputs, skip_ikr=skip_ikr)
    if not skip_ikr:
        ComparisonController.get_drug_effect(parallel=False)

    try:
        # Simulate APs and APD90s of the AP-SD model and the AP-CS model
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

    # Post-process the data
    outcome_df = ComparisonController.process_data()

    return outcome_df


# Load previously saved parameter space
param_space_fpath = os.path.join(data_dir, 'parameter_space',
                                 f'parameter_space_{IKrmodel}.csv')
param_dict = pd.read_csv(param_space_fpath, header=[0, 1], index_col=[0],
                         skipinitialspace=True)
param_dict = param_dict.rename(columns={"N": "n", "EC50": "halfmax"})
param_dict = param_dict.to_dict(orient='list')

skip_ikr = False if APmodel == 'ORd-Lei' else True

j = args.id_num
ind = param_dict[('param_id', 'param_id')].index(int(j))
param_values = {n: param_dict[('param_values', n)][ind]
                for n in param_names}
input = {'id': j, 'param_values': param_values}

# Use result from previous publication to reduce replicating simulations
if APmodel != 'ORd-Lei':
    key_list = [k for k in param_dict.keys()
                if 'drug_conc_Hill' in k]
    drug_conc = [param_dict[k][ind] for k in key_list]
    key_list = [k for k in param_dict.keys() if 'Hill_curve' in k]
    Hill_curve = {k[1]: param_dict[k][ind] for k in key_list}
    input.update({'drug_conc_Hill': drug_conc,
                  'Hill_curve': Hill_curve})

start = time.time()
output = param_evaluation(input, skip_ikr=skip_ikr)
print('Single run time: ', time.time() - start)

# Save output
file_prefix = 'SA_paramid_'
results_dir = os.path.join(data_dir, 'SA_space', APmodel)
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)
filename = f'{file_prefix}{args.id_num}.csv'
output.to_csv(os.path.join(results_dir, filename))
