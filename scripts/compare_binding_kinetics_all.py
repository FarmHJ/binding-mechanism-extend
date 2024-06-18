# Compare the AP, APD90 and qNet between the AP-Li-SD models and
# the AP-Li-CS models for dofetilide and verapamil
# Output:
# 1. 2 pulses of action potential simulated from both models (steady state).
# 2. APD90 of both pulses for both models (steady state).
# 3. qNet for both models.

import myokit
import numpy as np
import os
import pandas as pd

import modelling

model_details = model_details = modelling.model_naming
APmodel_list = model_details.APmodel_list[1:]
drug_list = ['dofetilide', 'verapamil']
model_drug_pair = [(m, d) for m in APmodel_list for d in drug_list]

for APmodel, drug in model_drug_pair:
    print('##############')
    print(APmodel, drug)
    # Define directories to save simulated data
    data_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                            APmodel, 'AP_duration_match', drug)
    if not os.path.isdir(data_dir):
        os.makedirs(data_dir)

    # Set AP model
    model_keys = modelling.model_naming.model_current_keys[APmodel]
    current_key = model_keys['IKr']
    APsim = modelling.ModelSimController(APmodel)

    # Get the Hill coef for Li-CS model
    if APmodel == 'ORd-Lei':
        ikr_model = 'Lei'
    else:
        APsim.set_ikr_rescale_method('AP_duration')
        ikr_model = 'Li'

    HillModel = modelling.HillModel()
    APsim.update_initial_state()

    Hill_coef = modelling.BindingParameters().load_Hill_eq(
        drug, ikr_model=ikr_model)
    # Hill_coef = modelling.BindingParameters().load_published_Hill_eq(
    #     drug, ikr_model=ikr_model, channel='IKr')
    # Define the range of drug concentration for a given drug
    drug_conc = modelling.SD_details.drug_concentrations[drug]['fine']

    # Simulate AP of the AP-SD model and the AP-CS model
    APD_conductance = []
    APD_trapping = []
    save_AP = 3
    EAD_present = False

    APsim.set_SD_parameters(drug, ikr_model=ikr_model)

    for i in range(len(drug_conc)):
        print('simulating concentration: ', i, drug_conc[i])

        APsim.set_conc(drug_conc[i])
        log = APsim.simulate(save_signal=2)
        if i % save_AP == 0:
            conc_str = "{0:.3e}".format(drug_conc[i])
            log.save_csv(os.path.join(data_dir,
                                      f'SD_AP_litHill_{conc_str}.csv'))

        apd90 = APsim.APD90(log)
        APD_trapping.append(apd90)

        reduction_scale = HillModel.simulate(Hill_coef, drug_conc[i])
        APsim.set_conc(0)
        APsim.set_CS_parameter(reduction_scale)
        log = APsim.simulate(save_signal=2)
        APsim.set_CS_parameter(1)
        if i % save_AP == 0:
            log.save_csv(os.path.join(data_dir,
                                      f'CS_AP_litHill_{conc_str}.csv'))

        apd90 = APsim.APD90(log)
        APD_conductance.append(apd90)

    # Save simulated APD90 of both the ORd-SD model and the ORd-CS model
    APD_SD_max = [float('nan') if np.isnan(i).any() else max(i)
                  for i in APD_trapping]
    APD_CS_max = [float('nan') if np.isnan(i).any() else max(i)
                  for i in APD_conductance]

    column_name = ['pulse 0', 'pulse 1']
    trapping_df = pd.DataFrame(APD_trapping, columns=column_name)
    trapping_df['APD'] = APD_SD_max
    trapping_df['drug concentration'] = drug_conc

    conductance_df = pd.DataFrame(APD_conductance, columns=column_name)
    conductance_df['APD'] = APD_CS_max
    conductance_df['drug concentration'] = drug_conc

    # Compute qNet
    pulse_time = 2000
    base_log_key = [APsim.time_key, APsim.Vm_key]
    qNet_current_list = [model_keys['ICaL'], model_keys['INaL'],
                         current_key, model_keys['IKs'],
                         model_keys['IK1'], model_keys['Ito']]
    qNet_current_list = [i for i in qNet_current_list if i is not None]

    APsim.sim.set_protocol(myokit.pacing.blocktrain(period=pulse_time,
                                                    duration=0.5,
                                                    offset=50))
    APsim._cycle_length = pulse_time
    APsim.reset_parameters()
    states = myokit.load_state(
        os.path.join(modelling.RESULT_DIR, 'steady_states',
                     f'{APmodel}_steadystate_qNetprot.csv'))
    APsim.set_initial_state(states)

    qNet_SD_arr = []
    qNet_CS_arr = []
    print('Computing qNet...')
    APsim.set_SD_parameters(drug, ikr_model=ikr_model)

    for i in range(len(drug_conc)):
        print(f'simulating concentration: {drug_conc[i]}')

        # Run simulation for the AP-SD model till steady state
        APsim.set_conc(drug_conc[i])
        log = APsim.simulate(save_signal=2, timestep=0.01,
                             log_var=base_log_key + qNet_current_list)

        # Check for EAD-like behaviour
        apd90 = APsim.APD90(log)
        if np.isnan(apd90).any():
            qNet = float('nan')
        else:
            qNet = APsim.qNet(log, period=2)
        qNet_SD_arr.append(qNet)

        # Run simulation for the AP-CS model till steady state
        reduction_scale = HillModel.simulate(Hill_coef, drug_conc[i])
        APsim.set_conc(0)
        APsim.set_CS_parameter(reduction_scale)
        log = APsim.simulate(save_signal=2, timestep=0.01,
                             log_var=base_log_key + qNet_current_list)
        APsim.set_CS_parameter(1)

        # Check for EAD-like behaviour
        apd90 = APsim.APD90(log)
        if np.isnan(apd90).any():
            qNet = float('nan')
        else:
            qNet = APsim.qNet(log, period=2)
        qNet_CS_arr.append(qNet)

    # Save qNet
    trapping_df['qNet'] = qNet_SD_arr
    conductance_df['qNet'] = qNet_CS_arr
    trapping_df.to_csv(os.path.join(data_dir, 'SD_APD_qNet_litHill.csv'))
    conductance_df.to_csv(os.path.join(data_dir, 'CS_APD_qNet_litHill.csv'))
