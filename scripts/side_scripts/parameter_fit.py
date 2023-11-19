import argparse
import glob
import matplotlib
import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import pandas as pd
import pints

import modelling


# Define AP model, drug and tuning method
parser = argparse.ArgumentParser(
    description="Fit parameters for Lei-SD model")
# parser.add_argument("protocol", help="Name of protocol",
#                     choices=['CIPA', 'Pharm', 'Milnes', 'staircase'])
# # parser.add_argument("--drug", default='all', help="Drug")
# parser.add_argument('--plot', action='store_true', help="Plot synthetic data")
parser.add_argument("--cache", action='store_true',
                    help='Use generated synthetic data')
args = parser.parse_args()

dataset = modelling.DataLibrary()
drug_list = modelling.SD_details.drug_names

results_dir = os.path.join(modelling.PARAM_DIR, 'Lei-SD-inference')
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)


# PINTS model
class InfModel(pints.ForwardModel):
    def __init__(self, conc):
        # , control_pulses=10
        self.sim = modelling.ModelSimController('Lei')
        # TODO: have to rescale conductance? Not for now
        self.sim.sim.set_state(control_state)
        self.sim.set_conc(conc)

    def n_outputs(self):
        return 1

    def n_parameters(self):
        return 5

    def simulate(self, parameters, times):

        param_df = pd.DataFrame(parameters,
                                index=modelling.SD_details.SD_param_names).T
        self.sim.set_parameters(param_df)

        # TODO: add pre-simulation of -80mV for 10s
        log_times = []
        for i in range(sweep_num):
            log_times.extend(window + 25e3 * i)
        log = self.sim.simulate(prepace=0, save_signal=sweep_num,
                                log_times=log_times,
                                log_var=[self.sim.time_key, self.sim.ikr_key],
                                reset=False)

        output = []
        for s in range(sweep_num):
            output += list(log[self.sim.ikr_key, s] / control_log[self.sim.ikr_key])
        time = np.arange(0, len(output)) * 10
        assert len(time) == len(times)

        return output


# Set up inference problem
drug = 'dofetilide'
sweep_num = 10

# Load data
dataset = modelling.DataLibrary()
# for drug in drug_list:
drug = drug_list[0]
dataset.set_drug(drug)
conc_list = dataset.concs

signal = dataset.get_mean_signal(cache=True)
dataset.data = signal.loc[signal['time'] > 1000]

# Prepare data
dataset.set_conc(conc_list[0])
frac_block = []
for s in range(dataset.n_sweeps):
    frac_block.extend(dataset.data.loc[dataset.data['sweep'] == s + 1,
                                       'frac'].values)

print(len(frac_block))
plt.figure()
plt.plot(frac_block, 'k--')
plt.savefig(os.path.join(modelling.FIG_DIR, 'checking_figures',
                         'exp_data.pdf'))

time = np.arange(0, len(frac_block)) * 10
window = dataset.data.loc[dataset.data['sweep'] == s + 1, 'time'].values

# Prepare control data
if args.cache:
    control_state = myokit.load_state(
        os.path.join(results_dir, 'control_state.csv'))
    control_log = myokit.DataLog.load_csv(os.path.join(results_dir,
                                                       'control_log.csv'))

    plt.figure()
    plt.plot(control_log.time(), control_log['ikr.IKr'])
    plt.savefig(os.path.join(modelling.FIG_DIR, 'checking_figures',
                             'control_log.pdf'))
else:
    # Simulate control condition
    sim = modelling.ModelSimController('Lei')
    sim.set_conc(0)
    control_log = sim.simulate(prepace=999, log_var=[sim.time_key, sim.ikr_key],
                               log_times=window)

    control_state = sim.sim.state()
    control_log.save_csv(os.path.join(results_dir, 'control_log.csv'))
    myokit.save_state(os.path.join(results_dir, 'control_state.csv'),
                      control_state)
    del sim

# TODO: optimise for all concentrations
# Instantiate forward model
model = InfModel(conc_list[0])
# Create single output problem
problem = pints.SingleOutputProblem(model, time, frac_block)

# Error function
error_fn = pints.RootMeanSquaredError(problem)

# Transform parameters
transform = pints.ComposedTransformation(
    pints.LogTransformation(error_fn.n_parameters() - 2),
    pints.IdentityTransformation(2),
)

# prior
# SD_param_names = ['Kmax', 'Ku', 'halfmax', 'n', 'Vhalf']
# prior_list = [pints.UniformLogPrior(1e-1, 1e8),
#               pints.UniformLogPrior(1e-5, 1e1),
#               pints.UniformLogPrior(1e0, 1e9),
#               pints.UniformLogPrior(1e-1, 1e8)]
prior_list = [pints.UniformLogPrior(1e5, 1e10),
              pints.UniformLogPrior(1e-7, 1e-5),
              pints.UniformLogPrior(1e7, 1e9)]
# print(prior_list[1].sample())
guess = [p.sample()[0][0] for p in prior_list]
guess.append(np.random.uniform(0, 1))
guess.append(np.random.uniform(-100, 0))

opt = pints.OptimisationController(error_fn, guess, transformation=transform,
                                   method=pints.CMAES)
opt.set_log_to_file(os.path.join(results_dir, drug + '.txt'))
opt.set_parallel(True)

p, s = opt.run()
print(p)

output = pd.DataFrame(p, index=modelling.SD_details.SD_param_names).T
output['error'] = s
output.to_csv(os.path.join(results_dir, f'{drug}-Milnes-conc{conc_list[0]}.csv'))

log = model.simulate(p, time)
plt.figure()
plt.plot(time, frac_block)
plt.plot(time, log)
plt.savefig(os.path.join(modelling.FIG_DIR, 'checking_figures',
                         'opt_result.pdf'))
