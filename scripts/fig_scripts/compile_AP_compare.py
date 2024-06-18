# Compares the IKr, AP and APD90 of the SD model and the CS model
import argparse
import numpy as np
import os
import pandas as pd

import modelling


parser = argparse.ArgumentParser(
    description="Comparison between AP-IKr-SD model and the AP-IKr-CS model")
parser.add_argument("--mode", default='cropped', choices=['cropped', 'full'],
                    help="Cropped or full image")
parser.add_argument("--metric", default='RMSE',
                    choices=['RMSE', 'ME', 'RMSE_ratio', 'ME_ratio'],
                    help="Choose metric to plot")
args = parser.parse_args()

model_details = modelling.model_naming
model_list = model_details.APmodel_list
model_label = ['ORd-Li', 'Grandi-Li', 'ten Tusscher-Li', 'Tomek-Li', 'ORd-Lei']
drug_list = ['dofetilide', 'verapamil']
tuning_method = ['AP_duration'] * len(model_list)

fig_dir = os.path.join(modelling.FIG_DIR, 'kinetics_comparison')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

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

FnClass = modelling.ModelComparison(modelling.ModelSimController('Grandi'))

for num, APmodel_name in enumerate(model_list):
    print('##################')
    print(APmodel_name)
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
        APD_trapping = APD_trapping_df['APD'].values.tolist()
        APD_conductance = APD_conductance_df['APD'].values.tolist()
        EAD_marker = [False if (np.isnan((i, j)).any()) else True for (i, j)
                      in zip(APD_trapping, APD_conductance)]

        APD_trapping = np.array(APD_trapping)[np.array(EAD_marker)]
        APD_conductance = np.array(APD_conductance)[np.array(EAD_marker)]

        APD_min = min(min(APD_trapping), min(APD_conductance))
        APD_max = max(max(APD_trapping), max(APD_conductance))
        APD_SD_norm = (APD_trapping - APD_min) / (APD_max - APD_min)
        APD_CS_norm = (APD_conductance - APD_min) / (APD_max - APD_min)
        # FnClass.APD_trapping = APD_trapping
        # FnClass.APD_conductance = APD_conductance
        FnClass.APD_trapping = APD_SD_norm
        FnClass.APD_conductance = APD_CS_norm
        getattr(FnClass, args.metric)()
        RMSD = FnClass.Error
        getattr(FnClass, 'ME')()
        MD = FnClass.Error
        # APD_metric[drug].append(FnClass.Error)
        APD_metric[drug].append(RMSD * MD / np.abs(MD))

        qNet_SD_df = pd.read_csv(os.path.join(data_dir, 'SD_qNet_EAD.csv'))
        qNet_CS_df = pd.read_csv(os.path.join(data_dir, 'CS_qNet_EAD.csv'))

        SD_qNet = qNet_SD_df['qNet'].values.tolist()
        CS_qNet = qNet_CS_df['qNet'].values.tolist()
        EAD_marker = [False if (np.isnan((i, j)).any()) else True for (i, j)
                      in zip(SD_qNet, CS_qNet)]

        SD_qNet = np.array(SD_qNet)[np.array(EAD_marker)]
        CS_qNet = np.array(CS_qNet)[np.array(EAD_marker)]

        # FnClass.APD_trapping = SD_qNet
        # FnClass.APD_conductance = CS_qNet
        qnet_min = min(min(SD_qNet), min(CS_qNet))
        qnet_max = max(max(SD_qNet), max(CS_qNet))
        qnet_SD_norm = (SD_qNet - qnet_min) / (qnet_max - qnet_min)
        qnet_CS_norm = (CS_qNet - qnet_min) / (qnet_max - qnet_min)
        FnClass.APD_trapping = qnet_SD_norm
        FnClass.APD_conductance = qnet_CS_norm

        getattr(FnClass, args.metric)()
        RMSD = FnClass.Error
        getattr(FnClass, 'ME')()
        MD = FnClass.Error
        # qNet_metric[drug].append(FnClass.Error)
        qNet_metric[drug].append(RMSD * MD / np.abs(MD))
        # print('qNet')
        # print(f'{drug}: {FnClass.Error}')

x = np.arange(len(model_list))
bar_width = 0.3

plot = modelling.figures.FigurePlot()
if args.mode == 'cropped':
    fig = modelling.figures.FigureStructure(figsize=(8, 3), gridspec=(2, 2),
                                            height_ratios=[1, 3], wspace=0.27,
                                            hspace=0.15)
elif args.mode == 'full':
    fig = modelling.figures.FigureStructure(figsize=(8, 3), gridspec=(1, 2),
                                            wspace=0.27)

multiplier = 0
for drug, metric in APD_metric.items():
    offset = bar_width * multiplier
    fig.axs[0][0].bar(x + offset, metric, bar_width,
                      label=drug_label[drug], color=drug_color[drug])
    if args.mode == 'cropped':
        fig.axs[1][0].bar(x + offset, metric, bar_width,
                          label=drug_label[drug], color=drug_color[drug])
    multiplier += 1

multiplier = 0
for drug, metric in qNet_metric.items():
    offset = bar_width * multiplier
    fig.axs[0][1].bar(x + offset, metric, bar_width,
                      label=drug_label[drug], color=drug_color[drug])
    if args.mode == 'cropped':
        fig.axs[1][1].bar(x + offset, metric, bar_width,
                          label=drug_label[drug], color=drug_color[drug])
    multiplier += 1

if args.mode == 'cropped':
    y_bottom, y_top = fig.axs[0][0].get_ylim()
    fig.axs[0][0].set_ylim(18, y_top)
    fig.axs[1][0].set_ylim(y_bottom, 14)
    fig.axs[0][0].spines['bottom'].set_visible(False)
    fig.axs[1][0].spines['top'].set_visible(False)
    fig.axs[0][0].tick_params(bottom=False, labelbottom=False)
    fig.axs[1][0].xaxis.tick_bottom()

    y_bottom, y_top = fig.axs[0][1].get_ylim()
    fig.axs[0][1].set_ylim(0.017, y_top)
    fig.axs[1][1].set_ylim(y_bottom, 0.014)
    fig.axs[0][1].spines['bottom'].set_visible(False)
    fig.axs[1][1].spines['top'].set_visible(False)
    fig.axs[0][1].tick_params(bottom=False, labelbottom=False)
    fig.axs[1][1].xaxis.tick_bottom()

    # Add diagonal lines on yaxis
    d = 0.015
    for i in range(2):
        kwargs = dict(transform=fig.axs[0][i].transAxes, color='k',
                      clip_on=False)
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
# fig.axs[0][0].set_title(r'APD$_{90}$')
# fig.axs[0][1].set_title('qNet')

fig.axs[0][0].set_ylabel(r'$\Delta \widetilde{\mathrm{APD}}_{90}$')
fig.axs[0][1].set_ylabel(r'$\Delta \widetilde{\mathrm{qNet}}$')

r = 0 if args.mode == 'full' else 1
for i in range(2):
    # fig.axs[r][i].set_ylabel(f'{args.metric}')
    fig.axs[r][i].set_xticks(x + bar_width / 2, model_label,
                             rotation=30, ha='right')
    fig.axs[r][i].spines[['right', 'top']].set_visible(False)

    if args.metric in ['ME', 'ME_ratio']:
        fig.axs[r][i].axhline(0, 0, 1, color='k', lw=0.5)

fig.axs[0][0].legend(handlelength=1)

# Add panel letter
fig.fig.text(0.1, 0.905, '(A)', fontsize=11)
fig.fig.text(0.52, 0.905, '(B)', fontsize=11)

fig.savefig(os.path.join(fig_dir,
                         f"model_compare_all_{args.mode}_normsigned{args.metric}_EAD2.pdf"))  # ,
            # format='svg')
