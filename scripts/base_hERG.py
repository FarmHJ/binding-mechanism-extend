# Generate an AP clamp protocol from ORd-Li model and subject it to the Li and
# the Lei model
import os

import modelling

data_dir = os.path.join(modelling.RESULT_DIR, 'background')
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

model_list = modelling.model_naming.APmodel_list[:-1]
for APmodel in model_list:
    # Load AP models
    APsim = modelling.ModelSimController(APmodel)
    if APmodel != 'ORd-Li':
        APsim.set_ikr_rescale_method('AP_duration')

    # Simulate AP for AP clamp protocol
    APclamp = APsim.simulate(log_var=[APsim.time_key, APsim.Vm_key])
    APclamp.save_csv(os.path.join(data_dir, f'APclamp_{APmodel}.csv'))

    # Set up AP clamp protocol
    times = APclamp[APsim.time_key]
    voltages = APclamp[APsim.Vm_key]

    # Load IKr model
    model_list = ['Li', 'Lei']
    for model_name in model_list:
        IKr_sim = modelling.ModelSimController(model_name)
        IKr_sim.set_protocol(times, voltages)
        log = IKr_sim.simulate(log_var='all')

        log.save_csv(os.path.join(
            data_dir, f'{APmodel}_APclamp_{model_name}_current.csv'))
