# Compares the IKr, AP and APD90 of the SD model and the CS model
import myokit
import numpy as np
import os
import pandas as pd

import modelling


# Define protocol
pulse_time = 1000
protocol_offset = 50
protocol = myokit.pacing.blocktrain(pulse_time, 0.5, offset=protocol_offset)

# Define constants
repeats = 1000
abs_tol = 1e-7
rel_tol = 1e-8

model_details = modelling.model_naming
model_list = model_details.APmodel_list
model_label = ['ORd-Li', 'Grandi-Li', 'ten Tusscher-Li', 'Tomek-Li', 'ORd-Lei']
drug_list = ['dofetilide', 'verapamil']
# tuning_method = ['hERG_peak'] * 4 + ['AP_duration']
tuning_method = ['AP_duration'] * 5

fig_dir = os.path.join(modelling.FIG_DIR, 'kinetics_comparison')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

drug_color = {
    'dofetilide': 'k',
    'verapamil': 'r'
}

fig = modelling.figures.FigureStructure(figsize=(8, 3), gridspec=(2, 2),
                                        height_ratios=[1, 3], wspace=0.25,
                                        hspace=0.15)
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
        data_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                                APmodel_name, f'{tuning_method[num]}_match',
                                drug)

        SD_labelname = f'{APmodel_name}-SD model'
        CS_labelname = f'{APmodel_name}-CS model'

        # Load APD data
        APD_trapping_df = pd.read_csv(os.path.join(data_dir, 'SD_APD_fine.csv'))
        APD_conductance_df = pd.read_csv(os.path.join(data_dir, 'CS_APD_fine.csv'))

        # Identify EAD-like behaviour
        drug_conc = APD_trapping_df['drug concentration'].values.tolist()
        APD_trapping = APD_trapping_df['APD'].values.tolist()
        APD_conductance = APD_conductance_df['APD'].values.tolist()
        EAD_marker = [False if (np.isnan((i, j)).any()) else True for (i, j)
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

        qNet_SD_df = pd.read_csv(os.path.join(data_dir, 'SD_qNet.csv'))
        qNet_CS_df = pd.read_csv(os.path.join(data_dir, 'CS_qNet.csv'))

        SD_qNet = qNet_SD_df['qNet'].values.tolist()
        CS_qNet = qNet_CS_df['qNet'].values.tolist()

        SD_qNet = np.array(SD_qNet)[np.array(EAD_marker)]
        CS_qNet = np.array(CS_qNet)[np.array(EAD_marker)]

        square_sum = 0
        count = 0
        for i in range(len(SD_qNet)):
            square_sum += (SD_qNet[i] - CS_qNet[i])**2
            count += 1
        RMSDiff_qNet = np.sqrt(square_sum) / count
        qNet_metric[drug].append(RMSDiff_qNet)

x = np.arange(len(model_list))
bar_width = 0.3

multiplier = 0
for drug, metric in APD_metric.items():
    offset = bar_width * multiplier
    fig.axs[0][0].bar(x + offset, metric, bar_width,
                      label=drug, color=drug_color[drug])
    fig.axs[1][0].bar(x + offset, metric, bar_width,
                      label=drug, color=drug_color[drug])
    multiplier += 1

multiplier = 0
for drug, metric in qNet_metric.items():
    offset = bar_width * multiplier
    fig.axs[0][1].bar(x + offset, metric, bar_width,
                      label=drug, color=drug_color[drug])
    fig.axs[1][1].bar(x + offset, metric, bar_width,
                      label=drug, color=drug_color[drug])
    multiplier += 1

y_bottom, y_top = fig.axs[0][0].get_ylim()
fig.axs[0][0].set_ylim(30, y_top)
fig.axs[1][0].set_ylim(y_bottom, 15)
fig.axs[0][0].spines['bottom'].set_visible(False)
fig.axs[1][0].spines['top'].set_visible(False)
# fig.axs[0][0].xaxis.tick_top()
fig.axs[0][0].tick_params(bottom=False, labelbottom=False)
fig.axs[1][0].xaxis.tick_bottom()

y_bottom, y_top = fig.axs[0][1].get_ylim()
fig.axs[0][1].set_ylim(0.017, y_top)
fig.axs[1][1].set_ylim(y_bottom, 0.012)
fig.axs[0][1].spines['bottom'].set_visible(False)
fig.axs[1][1].spines['top'].set_visible(False)
# fig.axs[0][1].xaxis.tick_top()
fig.axs[0][1].tick_params(bottom=False, labelbottom=False)
fig.axs[1][1].xaxis.tick_bottom()

# Add diagonal lines on yaxis
d = 0.015
for i in range(2):
    kwargs = dict(transform=fig.axs[0][i].transAxes, color='k', clip_on=False)
    # top-left diagonal
    fig.axs[0][i].plot((-d, +d), (-d, +d), **kwargs)
    # top-right diagonal
    fig.axs[0][i].plot((1 - d, 1 + d), (-d, +d), **kwargs)
    kwargs.update(transform=fig.axs[1][i].transAxes)
    # bottom-left diagonal
    fig.axs[1][i].plot((-d, +d), (1 - d / 2, 1 + d / 2), **kwargs)
    # bottom-right diagonal
    fig.axs[1][i].plot((1 - d, 1 + d), (1 - d / 2, 1 + d / 2), **kwargs)

# Adjust figure
fig.axs[0][0].set_title(r'APD$_{90}$')
fig.axs[0][1].set_title('qNet')
for i in range(2):
    fig.axs[1][i].set_ylabel('RMSD')
    fig.axs[1][i].set_xticks(x + bar_width / 2, model_label,
                             rotation=30, ha='right')
    # axs[1].set_xticklabels(axs[1].get_xticks(), rotation=45)
fig.axs[0][0].legend(handlelength=1)

# Add panel letter
fig.fig.text(0.1, 0.905, '(A)', fontsize=11)
fig.fig.text(0.52, 0.905, '(B)', fontsize=11)

fig.savefig(os.path.join(fig_dir, "model_compare_all.svg"), format='svg')
