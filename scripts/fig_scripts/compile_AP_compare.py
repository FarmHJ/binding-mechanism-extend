# Compares the IKr, AP and APD90 of the SD model and the CS model
import matplotlib.pyplot as plt
import myokit
import numpy as np
import os
import pandas as pd
import sys

import modelling


# Define protocol
pulse_time = 1000
protocol_offset = 50
protocol = myokit.pacing.blocktrain(pulse_time, 0.5, offset=protocol_offset)

# Define constants
repeats = 1000
abs_tol = 1e-7
rel_tol = 1e-8

model_details = modelling.ModelDetails()
model_list = ['Grandi', 'TTP', 'Tomek', 'Lei']
drug_list = ['dofetilide', 'verapamil']
tuning_method = ['hERG_peak'] * 3 + ['AP_duration']

fig_dir = '../../figures/kinetics_comparison/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

drug_color = ['k', 'r']

fig = plt.figure()
for num, APmodel_name in enumerate(model_list):

    for drug_ind, drug in enumerate(drug_list):
        data_dir = \
            '../../simulation_data/kinetics_comparison/' + APmodel_name + '/' + \
            tuning_method[num] + '_match/' + drug + '/'

        SD_labelname = APmodel_name + '-SD model'
        CS_labelname = APmodel_name + '-CS model'

        # Load APD data
        APD_trapping_df = pd.read_csv(data_dir + 'SD_APD_fine.csv')
        APD_conductance_df = pd.read_csv(data_dir + 'CS_APD_fine.csv')

        # Identify EAD-like behaviour
        drug_conc = APD_trapping_df['drug concentration'].values.tolist()
        APD_trapping = APD_trapping_df['APD'].values.tolist()
        APD_conductance = APD_conductance_df['APD'].values.tolist()
        EAD_marker = [False if (i >= 1000 or j >= 1000) else True for (i, j)
                      in zip(APD_trapping, APD_conductance)]

        APD_trapping = np.array(APD_trapping)[np.array(EAD_marker)]
        APD_conductance = np.array(APD_conductance)[np.array(EAD_marker)]

        square_sum = 0
        count = 0
        for i in range(len(APD_trapping)):
            square_sum += (APD_trapping[i] - APD_conductance[i])**2
            count += 1
        RMSDiff_APD = np.sqrt(square_sum) / count

        plt.bar(num * 3 + drug_ind, RMSDiff_APD, color=drug_color[drug_ind])

        # SD_qNet = APD_trapping_df['qNet'].values.tolist()
        # CS_qNet = APD_conductance_df['qNet'].values.tolist()

        # SD_qNet = np.array(SD_qNet)[np.array(EAD_marker)]
        # CS_qNet = np.array(CS_qNet)[np.array(EAD_marker)]

        # square_sum = 0
        # count = 0
        # for i in range(len(SD_qNet)):
        #     square_sum += (SD_qNet[i] - CS_qNet[i])**2
        #     count += 1
        # RMSDiff_qNet = np.sqrt(square_sum) / count

        # plt.bar(num * 3 + 2, RMSDiff_qNet)

# # # Add panel letter
# # # fig.fig.set_size_inches(10, 5.5)
# # fig.fig.text(0.1, 0.905, '(A)', fontsize=11)
# # fig.fig.text(0.1, 0.455, '(B)', fontsize=11)
# # # fig.fig.text(0.49, 0.905, '(B)', fontsize=11)
# # fig.fig.text(0.53, 0.455, '(C)', fontsize=11)

# fig.savefig(fig_dir + "model_compare_all.pdf")
plt.savefig(fig_dir + "model_compare_all.pdf")
