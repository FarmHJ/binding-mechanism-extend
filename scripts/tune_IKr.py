# Introduces the idea of trapping and justifies the use of the Milnes protocol
import numpy as np
import os
import pandas as pd
from scipy.optimize import minimize

import modelling

# Define AP model
model_list = modelling.model_naming.APmodel_list[1:]


def APD_match(conductance_scale):
    APsim.set_ikr_rescale(conductance_scale)
    log = APsim.simulate(log_var=[APsim.time_key, APsim.Vm_key])
    apd90 = APsim.APD90(log)

    error = np.sqrt(np.power(apd90 - base_apd90, 2))

    return error


ikr_scale_dict = {}
for APmodel in model_list:
    model_keys = modelling.model_naming.model_current_keys[APmodel]

    base_APsim = modelling.ModelSimController(APmodel, ikr_modified=False)
    log = base_APsim.simulate(log_var='all')
    base_apd90 = base_APsim.APD90(log)
    base_peak = max(log[base_APsim.ikr_key])

    APsim = modelling.ModelSimController(APmodel)
    # Get peak ratio as initial guess
    log1 = APsim.simulate()
    ikr_peak_ratio = base_peak / max(log1[APsim.ikr_key])
    del log1

    # Get the IKr scaling factor by matching the APD
    initial_guess = ikr_peak_ratio
    res = minimize(APD_match, initial_guess, method='nelder-mead',
                   options={'disp': True})
    conductance_scale = res.x[0]
    print('AP_duration: ', conductance_scale)

    ikr_scale_dict.update({APmodel: conductance_scale})

df = pd.DataFrame.from_dict(ikr_scale_dict, orient='index',
                            columns=['AP_duration'])
result_file = os.path.join(modelling.RESULT_DIR, 'ikr_calibration.csv')
df.to_csv(result_file)
