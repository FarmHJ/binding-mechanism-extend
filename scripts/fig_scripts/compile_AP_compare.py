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

fig = modelling.figures.FigureStructure(figsize=(8, 3), gridspec=(1, 2),
                                        wspace=0.25)
plot = modelling.figures.FigurePlot()

APD_metric = {
    'dofetilide': [],
    'verapamil': []
}
qNet_metric = {
    'dofetilide': [],
    'verapamil': []
}
for num, APmodel_name in enumerate(model_list):
    for drug_ind, drug in enumerate(drug_list):
        data_dir = \
            '../../simulation_data/kinetics_comparison/' + APmodel_name + \
            '/' + tuning_method[num] + '_match/' + drug + '/'

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
        APD_metric[drug].append(RMSDiff_APD)

        # fig.axs[0][0].bar(num * 3 + drug_ind, RMSDiff_APD,
        #                   color=drug_color[drug_ind], label=drug)

        SD_qNet = APD_trapping_df['qNet'].values.tolist()
        CS_qNet = APD_conductance_df['qNet'].values.tolist()

        SD_qNet = np.array(SD_qNet)[np.array(EAD_marker)]
        CS_qNet = np.array(CS_qNet)[np.array(EAD_marker)]

        square_sum = 0
        count = 0
        for i in range(len(SD_qNet)):
            square_sum += (SD_qNet[i] - CS_qNet[i])**2
            count += 1
        RMSDiff_qNet = np.sqrt(square_sum) / count
        qNet_metric[drug].append(RMSDiff_qNet)

        # fig.axs[0][1].bar(num * 3 + drug_ind, RMSDiff_qNet,
        #                   color=drug_color[drug_ind], label=drug)


x = np.arange(len(model_list))
bar_width = 0.3

multiplier = 0
for drug, metric in APD_metric.items():
    offset = bar_width * multiplier
    fig.axs[0][0].bar(x + offset, metric, bar_width,
                      label=drug)
    multiplier += 1

multiplier = 0
for drug, metric in qNet_metric.items():
    offset = bar_width * multiplier
    fig.axs[0][1].bar(x + offset, metric, bar_width,
                      label=drug)
    multiplier += 1

# Adjust figure
fig.axs[0][0].set_title(r'APD$_{90}$')
fig.axs[0][1].set_title('qNet')
for i in range(2):
    fig.axs[0][i].set_ylabel('RMSD')
    fig.axs[0][i].set_xticks(x + bar_width / 2, model_list)
fig.axs[0][0].legend(handlelength=1)

# Add panel letter
fig.fig.text(0.1, 0.905, '(A)', fontsize=11)
fig.fig.text(0.52, 0.905, '(B)', fontsize=11)

fig.savefig(fig_dir + "model_compare_all.pdf")
