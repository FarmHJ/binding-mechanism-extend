import myokit
import os

import modelling

# Simulating the SD model with AP clamp protocol for 1000 pulses till steady
# state under drug free, addition of example drug T and example drug N
# conditions

data_dir = '../simulation_data/background/'

# Load AP models
model_details = modelling.ModelDetails()

model_filenames = model_details.file_names['ORd-CiPA']
APmodel = '../' + model_filenames['AP_SD_path']
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, base_constant=None)

# Define current protocol
pulse_time = 1000
protocol_offset = 50
protocol = myokit.pacing.blocktrain(pulse_time, 0.5, offset=protocol_offset)
AP_model.protocol = protocol

# Simulation constants
repeats = 1000
abs_tol = 1e-7
rel_tol = 1e-8

# Simulate AP for AP clamp protocol
APclamp = AP_model.model_simulation(repeats, abs_tol=abs_tol, rel_tol=rel_tol)
APclamp.save_csv(data_dir + 'APclamp.csv')

# Set up AP clamp protocol
times = APclamp['engine.time']
voltages = APclamp['membrane.V']
tmax = times[-1] + 1

# Load IKr model
model_list = ['ORd-CiPA', 'Lei']
for model_name in model_list:
    model_filenames = model_details.file_names[model_name]
    model = '../' + model_filenames['IKr_path']

    model, _, x = myokit.load(model)
    current_model = modelling.Simulation(model, current_head_key='ikr',
                                         base_constant=None)

    # Simulate SD model with AP clamp protocol under drug free condition
    log_control = current_model.drug_APclamp(times, voltages, tmax, repeats,
                                             abs_tol=abs_tol, rel_tol=rel_tol)
    log_control.save_csv(data_dir + model_name + '_APclamp_current.csv')
