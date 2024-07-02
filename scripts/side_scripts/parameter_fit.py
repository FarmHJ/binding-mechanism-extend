import myokit
import numpy as np
import os
import pandas as pd
import pints

import modelling


# Set up directory to save results
results_dir = os.path.join(modelling.PARAM_DIR, 'Lei-SD-inference')
if not os.path.isdir(results_dir):
    os.makedirs(results_dir)


# PINTS model
class InfModel(pints.ForwardModel):
    def __init__(self):
        # Set up model
        self.sim = modelling.ModelSimController('Lei')
        self.sim.set_ikr_rescale_method('AP_duration')
        self.sim.initial_state = control_state

    def n_outputs(self):
        return 1

    def n_parameters(self):
        return 5

    def set_conc(self, conc):
        self.conc = conc

    def simulate(self, parameters, times):

        self.sim.update_initial_state(paces=0)
        param_values = {name: parameters[n] for n, name in
                        enumerate(modelling.SD_details.SD_param_names)}
        self.sim.set_parameters(param_values)
        self.sim.set_conc(self.conc)

        # TODO: add pre-simulation of -80mV for 10s
        log = self.sim.simulate(prepace=0, save_signal=sweep_num,
                                log_times=times,
                                log_var=[self.sim.time_key, self.sim.ikr_key],
                                reset=False)

        output = []
        for s in range(sweep_num):
            output += list(log[self.sim.ikr_key, s] /
                           control_log_win[self.sim.ikr_key])

        return output


# Prepare control data
dataset = modelling.DataLibrary()
general_win = np.arange(1005 + 100, 1005 + 10000, 10)
if os.path.isdir(os.path.join(results_dir, 'control_state.csv')) and \
   os.path.isdir(os.path.join(results_dir, 'control_log.csv')):
    control_state = myokit.load_state(
        os.path.join(results_dir, 'control_state.csv'))
    control_log = myokit.DataLog.load_csv(
        os.path.join(results_dir, 'control_log.csv'))

else:
    # Simulate control condition
    sim = modelling.ModelSimController('Lei')
    sim.set_ikr_rescale_method('AP_duration')
    sim.set_conc(0)
    control_log = sim.simulate(prepace=999,
                               log_var=[sim.time_key, sim.ikr_key],
                               log_times=general_win)

    control_state = sim.sim.state()
    control_log.save_csv(os.path.join(results_dir, 'control_log.csv'))
    myokit.save_state(os.path.join(results_dir, 'control_state.csv'),
                      control_state)
    del sim

# Set up inference problem
sweep_num = 10
ref_time = 945 + 100
drug_list = modelling.SD_details.drug_names
for drug in drug_list:
    print(drug)
    # Get drug concentrations and mean of IKr from Li et al.'s experimental
    # data
    dataset.set_drug(drug)
    conc_list = dataset.concs

    signal = dataset.get_mean_signal(cache=True)
    dataset.data = signal

    errors = []
    for c in conc_list:
        dataset.set_conc(c)

        # Define window of interest, to avoid the channel opening process
        # from interfering with the drug block development process, as
        # mentioned in Li et al.
        sweep_signal = dataset.conc_data.loc[(dataset.conc_data['sweep'] == 1),
                                             ['time', 'frac']]
        window = sweep_signal.loc[(sweep_signal['time'] >= ref_time),
                                  'time'].values

        log_times = []
        frac_block = []
        for s in range(dataset.n_sweeps):
            # Create times to log for simulated output
            log_times.extend(window + 60 + 25e3 * s)

            # Combine data of all pulses together
            sweep_data = dataset.conc_data.loc[
                dataset.conc_data['sweep'] == s + 1, ['time', 'frac']]
            frac_block.extend(sweep_data.loc[sweep_data['time'].isin(window),
                                             'frac'].values)
        control_log_win = control_log.trim(window[0] + 60,
                                           window[-1] + 60 + 10)

        # Instantiate forward model
        model = InfModel()
        model.set_conc(c)
        # Create single output problem
        problem = pints.SingleOutputProblem(model, log_times, frac_block)

        # Define error function
        errors.append(pints.RootMeanSquaredError(problem))

    error_fn = pints.SumOfErrors(errors, [1 / len(conc_list)] * len(conc_list))
    # Transform parameters
    transform = pints.ComposedTransformation(
        pints.LogTransformation(error_fn.n_parameters() - 2),
        pints.IdentityTransformation(2),
    )

    # Define prior and set boundaries
    # SD_param_names = ['Kmax', 'Ku', 'halfmax', 'n', 'Vhalf']
    prior_list = [pints.UniformLogPrior(1e5, 1e10),
                  pints.UniformLogPrior(1e-7, 1e-5),
                  pints.UniformLogPrior(1e7, 1e9)]
    boundaries = pints.RectangularBoundaries([1e-1, 1e-8, 1e0, 0, -200],
                                             [1e10, 1e1, 1e10, 2, 0])

    reps = 5
    param_scores = []
    for i in range(reps):
        # Generate random initial guess
        np.random.seed((i + 1) * 100)
        guess = [p.sample()[0][0] for p in prior_list]
        guess.append(np.random.uniform(0, 1))
        guess.append(np.random.uniform(-100, 0))

        # Run optimisation
        print('initial guess: ', guess)
        opt = pints.OptimisationController(error_fn, guess,
                                           boundaries=boundaries,
                                           transformation=transform,
                                           method=pints.CMAES)
        opt.set_parallel(True)

        p, s = opt.run()
        param_scores.append(list(p) + [s])

    # Save optimised parameters
    scores_dict = pd.DataFrame(param_scores,
                               columns=modelling.SD_details.SD_param_names +
                               ['error'])
    scores_dict = scores_dict.sort_values(by=['error'])
    scores_dict.to_csv(os.path.join(results_dir, f'{drug}.csv'))
