# Generate an AP clamp protocol from ORd-Li model and subject it to the Li and
# the Lei model
import os

import modelling

data_dir = os.path.join(modelling.RESULT_DIR, 'background')

# Load AP models
APsim = modelling.ModelSimController('ORd-Li')

# Simulate AP for AP clamp protocol
APclamp = APsim.simulate(log_var=[APsim.time_key, APsim.Vm_key])
APclamp.save_csv(os.path.join(data_dir, 'APclamp.csv'))

# Set up AP clamp protocol
times = APclamp[APsim.time_key]
voltages = APclamp[APsim.Vm_key]

# Load IKr model
model_list = ['Li', 'Lei']
for model_name in model_list:
    IKr_sim = modelling.ModelSimController(model_name)
    IKr_sim.set_protocol(times, voltages)
    log = IKr_sim.simulate()

    log.save_csv(os.path.join(data_dir, model_name + '_APclamp_current.csv'))
