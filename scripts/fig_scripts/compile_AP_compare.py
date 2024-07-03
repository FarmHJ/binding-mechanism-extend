# Compares the IKr, AP and APD90 of the SD model and the CS model
import numpy as np
import os
import pandas as pd

import modelling


# Define all AP models and drugs for comparison
model_details = modelling.model_naming
model_list = model_details.APmodel_list
model_label = ['ORd-Li', 'Grandi-Li', 'ten Tusscher-Li', 'Tomek-Li', 'ORd-Lei']
drug_list = ['dofetilide', 'verapamil']

# Define directory to save figure
fig_dir = os.path.join(modelling.FIG_DIR, 'kinetics_comparison')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Define color and label for each drug
drug_color = {
    'dofetilide': 'k',
    'verapamil': 'r'
}
drug_label = {
    'dofetilide': 'example drug T',
    'verapamil': 'example drug N'
}

APD_metric = {
    'dofetilide': [],
    'verapamil': []
}
qNet_metric = {
    'dofetilide': [],
    'verapamil': []
}

# To use functions defined within the ModelSimController class
FnClass = modelling.ModelComparison(modelling.ModelSimController('Grandi'))

model_drug_pair = [(m, d) for m in model_list for d in drug_list]
for APmodel_name, drug in model_drug_pair:

    data_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                            APmodel_name, 'AP_duration_match',
                            drug)

    SD_labelname = f'{APmodel_name}-SD model'
    CS_labelname = f'{APmodel_name}-CS model'

    # Load APD data
    trapping_df = pd.read_csv(os.path.join(data_dir, 'SD_APD_qNet.csv'))
    conductance_df = pd.read_csv(os.path.join(data_dir, 'CS_APD_qNet.csv'))

    # Identify EAD-like behaviour
    APD_trapping = trapping_df['APD'].values.tolist()
    APD_conductance = conductance_df['APD'].values.tolist()
    EAD_marker = [False if (np.isnan((i, j)).any()) else True for (i, j)
                  in zip(APD_trapping, APD_conductance)]

    APD_trapping = np.array(APD_trapping)[np.array(EAD_marker)]
    APD_conductance = np.array(APD_conductance)[np.array(EAD_marker)]

    # Calculate metric to quantify difference between the AP-SD model and
    # the AP-CS model
    # First, normalise APD90s of both the AP-SD model and the AP-CS model for
    # all drug concentrations
    # Then, compute the RMSD and MD
    # Finally, multiply the sign value of MD on to the RMSD
    APD_min = min(min(APD_trapping), min(APD_conductance))
    APD_max = max(max(APD_trapping), max(APD_conductance))
    APD_SD_norm = (APD_trapping - APD_min) / (APD_max - APD_min)
    APD_CS_norm = (APD_conductance - APD_min) / (APD_max - APD_min)
    FnClass.APD_trapping = APD_SD_norm
    FnClass.APD_conductance = APD_CS_norm
    getattr(FnClass, 'RMSE')()
    RMSD = FnClass.Error
    getattr(FnClass, 'ME')()
    MD = FnClass.Error
    APD_metric[drug].append(RMSD * MD / np.abs(MD))

    # Take the qNet data
    SD_qNet = trapping_df['qNet'].values.tolist()
    CS_qNet = conductance_df['qNet'].values.tolist()

    SD_qNet = np.array(SD_qNet)[np.array(EAD_marker)]
    CS_qNet = np.array(CS_qNet)[np.array(EAD_marker)]

    # Compute the metric for qNet
    qnet_min = min(min(SD_qNet), min(CS_qNet))
    qnet_max = max(max(SD_qNet), max(CS_qNet))
    qnet_SD_norm = (SD_qNet - qnet_min) / (qnet_max - qnet_min)
    qnet_CS_norm = (CS_qNet - qnet_min) / (qnet_max - qnet_min)
    FnClass.APD_trapping = qnet_SD_norm
    FnClass.APD_conductance = qnet_CS_norm

    getattr(FnClass, 'RMSE')()
    RMSD = FnClass.Error
    getattr(FnClass, 'ME')()
    MD = FnClass.Error
    qNet_metric[drug].append(RMSD * MD / np.abs(MD))

# Set up figure for bar plot
plot = modelling.figures.FigurePlot()
fig = modelling.figures.FigureStructure(figsize=(8, 3), gridspec=(1, 2),
                                        wspace=0.27)
x = np.arange(len(model_list))
bar_width = 0.3

# Plot metric of APD90 and qNet for all AP models
multiplier = 0
for drug, metric in APD_metric.items():
    offset = bar_width * multiplier
    fig.axs[0][0].bar(x + offset, metric, bar_width,
                      label=drug_label[drug], color=drug_color[drug])
    multiplier += 1

multiplier = 0
for drug, metric in qNet_metric.items():
    offset = bar_width * multiplier
    fig.axs[0][1].bar(x + offset, metric, bar_width,
                      label=drug_label[drug], color=drug_color[drug])
    multiplier += 1

# Adjust figures
fig.axs[0][0].set_ylabel(r'$\Delta \widetilde{\mathrm{APD}}_{90}$')
fig.axs[0][1].set_ylabel(r'$\Delta \widetilde{\mathrm{qNet}}$')

for i in range(2):
    fig.axs[0][i].set_xticks(x + bar_width / 2, model_label,
                             rotation=30, ha='right')
    fig.axs[0][i].spines[['right', 'top']].set_visible(False)

    fig.axs[0][i].axhline(0, 0, 1, color='k', lw=0.5)

fig.axs[0][0].legend(handlelength=1)

# Add panel letter
fig.fig.text(0.1, 0.905, '(A)', fontsize=11)
fig.fig.text(0.52, 0.905, '(B)', fontsize=11)

# Save figure
fig.savefig(os.path.join(fig_dir,
                         "model_compare_all_normsignedRMSD.pdf"))
