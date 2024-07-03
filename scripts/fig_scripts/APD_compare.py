# Compares the IKr, AP, APD90 and qNet of the SD model and the CS model
import argparse
import glob
import myokit
import numpy as np
import os
import pandas as pd

import modelling

# Define AP model, drug and tuning method
parser = argparse.ArgumentParser(
    description="Comparison between AP-IKr-SD model and the AP-IKr-CS model")
parser.add_argument("APmodel", help="Name of AP model")
parser.add_argument("drug", help="Drug")
parser.add_argument("-c", "--conc", default="conc", type=str,
                    choices=['dimless', 'conc'],
                    help='Choose to convert concentration to dimensionless')
args = parser.parse_args()

APmodel_name = args.APmodel
drug = args.drug

# Define directories to read data and save plotted figures
data_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                        APmodel_name, 'AP_duration_match', drug)
fig_dir = os.path.join(modelling.FIG_DIR, 'kinetics_comparison',
                       APmodel_name, 'AP_duration_match', drug)
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Set up structure of the figure
fig = modelling.figures.FigureStructure(figsize=(8, 5), gridspec=(2, 2),
                                        height_ratios=[1, 1], hspace=0.4,
                                        wspace=0.25, plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(2, 2), (1, 1), (1, 1)]
subgs = []
subgs.append(fig.gs[:2].subgridspec(*subgridspecs[0], wspace=0.08,
                                    hspace=0.08))
for i in range(1, 3):
    subgs.append(fig.gs[i + 1].subgridspec(*subgridspecs[i], wspace=0.08,
                                           hspace=0.08))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
         subgridspecs[k][1])] for i in range(subgridspecs[k][0])]
       for k in range(len(subgs))]

# Top panels
# Plot APs and the corresponding IKr of the AP-SD model and the
# AP-CS model stimulated to steady state
panel1 = axs[0]

# Read files name of action potential data
fname = f'AP_{args.conc}_'
SD_data_files = glob.glob(os.path.join(data_dir, f'SD_{fname}*.csv'))
CS_data_files = glob.glob(os.path.join(data_dir, f'CS_{fname}*.csv'))
key_len = len(fname) + 3
conc_label_SD = [os.path.basename(fname)[key_len:-4]
                 for fname in SD_data_files]
drug_conc_SD = [float(os.path.basename(fname)[key_len:-4])
                for fname in SD_data_files]
conc_label_CS = [os.path.basename(fname)[key_len:-4]
                 for fname in CS_data_files]
drug_conc_CS = [float(os.path.basename(fname)[key_len:-4])
                for fname in CS_data_files]

# Sort files in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc_SD), key=lambda x: x[1])]
sort_ind_CS = [i[0] for i in
               sorted(enumerate(drug_conc_CS), key=lambda x: x[1])]
drug_conc = sorted(drug_conc_SD)
trapping_data_files = [SD_data_files[i] for i in sort_ind]
conductance_data_files = [CS_data_files[i] for i in sort_ind_CS]
conc_label = [conc_label_SD[i] for i in sort_ind]
conc_label = [r"$10^{:d}$".format(int(np.log10(float(i)))) if float(i) >= 1e3
              else i for i in conc_label]
labels = [f'{i} nM' for i in conc_label]

# Load action potential data
trapping_AP_log = []
conductance_AP_log = []
for i in range(len(trapping_data_files)):
    trapping_AP_log.append(myokit.DataLog.load_csv(
        os.path.join(data_dir, trapping_data_files[i])))
    conductance_AP_log.append(myokit.DataLog.load_csv(
        os.path.join(data_dir, conductance_data_files[i])))

# Load APD90s and qNets
fname = f'APD_qNet_{args.conc}.csv'
trapping_df = pd.read_csv(os.path.join(data_dir, f'SD_{fname}'))
conductance_df = pd.read_csv(os.path.join(data_dir, f'CS_{fname}'))

APD_trapping = trapping_df['APD'].values.tolist()
APD_conductance = conductance_df['APD'].values.tolist()

# Initiate constants and variables
model_keys = modelling.model_naming.model_current_keys[APmodel_name]
current_key = model_keys['IKr']
Vm_key = model_keys['Vm']
plotting_pulse_time = 1000 * 2

SD_labelname = APmodel_name + '-SD model'
CS_labelname = APmodel_name + '-CS model'

# Plot AP and IKr at various drug concentrations
plot.add_multiple_continuous(panel1[0][0], trapping_AP_log,
                             Vm_key, labels=labels)
plot.add_multiple_continuous(panel1[1][0], trapping_AP_log,
                             current_key, labels=labels)
plot.add_multiple_continuous(panel1[0][1], conductance_AP_log,
                             Vm_key, labels=labels)
plot.add_multiple_continuous(panel1[1][1], conductance_AP_log,
                             current_key, labels=labels)
panel1[0][0].set_title(SD_labelname)
panel1[0][1].set_title(CS_labelname)
panel1[1][1].legend(handlelength=0.9, ncol=2, columnspacing=0.9,
                    loc=(1.04, 0.05))

# Adjust axes
fig.sharex(['Time (ms)'] * 2, [(0, plotting_pulse_time)] * 2,
           axs=panel1, subgridspec=subgridspecs[0])
fig.sharey(['Voltage (mV)', 'Current (A/F)'],
           axs=panel1, subgridspec=subgridspecs[0])

# Bottom left panel
# Plots the APD90 calculated for both models with drugs at various drug
# concentrations
panel2 = axs[1][0][0]

# Load APD data
# Identify EAD-like behaviour
drug_conc = trapping_df['drug concentration'].values.tolist()
EAD_marker = [1050 if (np.isnan((i, j)).any()) else None for (i, j)
              in zip(APD_trapping, APD_conductance)]

# Plot APD90 of both models
panel2.plot(drug_conc, APD_trapping, 'o', color='orange', label=SD_labelname)
panel2.plot(drug_conc, APD_conductance, '^', color='blue', label=CS_labelname,
            alpha=0.8)
panel2.scatter(drug_conc, EAD_marker, marker=(5, 2), color='k',
               label='EAD-like AP')
panel2.set_xscale("log", nonpositive='clip')
panel2.set_xlabel('Drug concentration (nM)')
panel2.set_ylabel(r'APD$_{90}$ (ms)')
panel2.legend(handlelength=1)

# Bottom right panel
# Plots the qNet calculated for both models
panel3 = axs[2][0][0]

# Load qNet data
SD_qNet = trapping_df['qNet'].values.tolist()
CS_qNet = conductance_df['qNet'].values.tolist()

# EAD_indicator = np.array([1 if i is None else None for i in EAD_marker])
SD_qNet = [None if np.isnan(APD_trapping[i]) else SD_qNet[i] for i in
           range(len(SD_qNet))]
CS_qNet = [None if np.isnan(APD_conductance[i]) else CS_qNet[i] for i in
           range(len(CS_qNet))]

# Plot APD90 of both models
panel3.plot(drug_conc, SD_qNet, 'o', color='orange', label=SD_labelname)
panel3.plot(drug_conc, CS_qNet, '^', color='blue', label=CS_labelname,
            alpha=0.8)
panel3.set_xlabel("Drug concentration (nM)")
panel3.set_xscale("log", nonpositive='clip')
l_lim, r_lim = panel2.get_xlim()
panel3.set_xlim(left=l_lim, right=r_lim)
panel3.set_ylabel('qNet (C/F)')
panel3.legend(handlelength=1)

# Add panel letter
fig.fig.text(0.1, 0.905, '(A)', fontsize=11)
fig.fig.text(0.1, 0.455, '(B)', fontsize=11)
fig.fig.text(0.53, 0.455, '(C)', fontsize=11)

# Save figure
fname = f'model_compare_{args.conc}.svg'
fig.savefig(os.path.join(fig_dir, fname), format='svg')
