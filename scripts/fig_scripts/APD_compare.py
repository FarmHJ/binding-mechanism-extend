# Compares the IKr, AP and APD90 of the SD model and the CS model
import matplotlib
import myokit
import numpy as np
import os
import pandas as pd
import sys

import modelling

# Define AP model, drug and tuning
APmodel_name = sys.argv[1]
drug = sys.argv[2]
protocol_name = 'Milnes'
tuning_method = sys.argv[3]

# Define directories to read data and save plotted figures
data_dir = \
    '../../simulation_data/kinetics_comparison/' + APmodel_name + '/' + \
    tuning_method + '_match/' + drug + '/'
fig_dir = \
    '../../figures/kinetics_comparison/' + APmodel_name + '/' + \
    tuning_method + '_match/' + drug + '/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Set up structure of the figure
fig = modelling.figures.FigureStructure(figsize=(8, 5), gridspec=(2, 2),
                                        height_ratios=[1, 1], hspace=0.4,
                                        wspace=0.25,
                                        # width_ratios=[2.5, 1.2],
                                        plot_in_subgrid=True)
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
# axs.append([[fig.fig.add_subplot(subgs[2][0, 0])]])

# Bottom left panel
# Plot action potentials and the corresponding IKr of the AP-SD model and the
# AP-CS model stimulated to steady state
panel2 = axs[0]

# Read files name of action potential data
SD_data_files = [f for f in os.listdir(data_dir) if
                 f.startswith('SD_AP_') and not
                 f.startswith('SD_AP_tran')]
CS_data_files = [f for f in os.listdir(data_dir) if
                 f.startswith('CS_AP_') and not
                 f.startswith('CS_AP_tran')]
conc_label_SD = [fname[6:-4] for fname in SD_data_files]
drug_conc_SD = [float(fname[6:-4]) for fname in SD_data_files]
conc_label_CS = [fname[6:-4] for fname in CS_data_files]
drug_conc_CS = [float(fname[6:-4]) for fname in CS_data_files]

# Sort in increasing order of drug concentration
sort_ind = [i[0] for i in sorted(enumerate(drug_conc_SD), key=lambda x:x[1])]
sort_ind_CS = [i[0] for i in
               sorted(enumerate(drug_conc_CS), key=lambda x:x[1])]
drug_conc = sorted(drug_conc_SD)
trapping_data_files = [SD_data_files[i] for i in sort_ind]
conductance_data_files = [CS_data_files[i] for i in sort_ind_CS]
conc_label = [conc_label_SD[i] for i in sort_ind]
conc_label = [r"$10^{:d}$".format(int(np.log10(float(i)))) if float(i) >= 1e3
              else i for i in conc_label]

# Initiate constants and variables
labels = [i + ' nM' for i in conc_label]
cmap = matplotlib.cm.get_cmap('viridis')

# Load action potential and APD data
trapping_AP_log = []
conductance_AP_log = []
for i in range(len(trapping_data_files)):
    trapping_AP_log.append(myokit.DataLog.load_csv(
        data_dir + trapping_data_files[i]))
    conductance_AP_log.append(myokit.DataLog.load_csv(
        data_dir + conductance_data_files[i]))

APD_trapping = pd.read_csv(data_dir + 'SD_APD_pulses2.csv')
APD_conductance = pd.read_csv(data_dir + 'CS_APD_pulses2.csv')

APD_trapping = [max(APD_trapping.loc[i].values.tolist()[1:-1]) for i in
                range(APD_trapping.shape[0])]
APD_conductance = [max(APD_conductance.loc[i].values.tolist()[1:-1]) for i in
                   range(APD_conductance.shape[0])]

# Initiate constants and variables
plotting_pulse_time = 1000 * 2

# Remove repeated signals at high concentrations
second_EAD_trap = [i for i, e in enumerate(APD_trapping)
                   if e == 1000][1:]
second_EAD_conduct = [i for i, e in enumerate(APD_conductance)
                      if e == 1000][1:]

if len(second_EAD_trap) <= len(second_EAD_conduct):
    chosen_conc_ind = second_EAD_trap
else:
    chosen_conc_ind = second_EAD_conduct

AP_trapping_plot = [e for i, e in enumerate(trapping_AP_log)
                    if i not in chosen_conc_ind]
AP_conductance_plot = [e for i, e in enumerate(conductance_AP_log)
                       if i not in chosen_conc_ind]

model_keys = modelling.ModelDetails().current_keys[APmodel_name]
current_key = model_keys['IKr']
Vm_key = model_keys['Vm']

SD_labelname = APmodel_name + '-SD model'
CS_labelname = APmodel_name + '-CS model'

# Plot AP and IKr at various drug concentrations
plot.add_multiple_continuous(panel2[0][0], AP_trapping_plot,
                             Vm_key, cmap=cmap,
                             labels=labels)
plot.add_multiple_continuous(panel2[1][0], AP_trapping_plot,
                             current_key, cmap=cmap, labels=labels)
plot.add_multiple_continuous(panel2[0][1], AP_conductance_plot,
                             Vm_key, cmap=cmap,
                             labels=labels)
plot.add_multiple_continuous(panel2[1][1], AP_conductance_plot,
                             current_key, cmap=cmap, labels=labels)
panel2[0][0].set_title(SD_labelname)
panel2[0][1].set_title(CS_labelname)
panel2[1][1].legend(handlelength=0.9, ncol=2, columnspacing=0.9,
                    loc=(1.04, 0.05))

# Adjust axes
fig.sharex(['Time (ms)'] * 2, [(0, plotting_pulse_time)] * 2,
           axs=panel2, subgridspec=subgridspecs[0])
fig.sharey(['Voltage (mV)', 'Current (A/F)'],
           axs=panel2, subgridspec=subgridspecs[0])

# Top right panel
# Plots the APD90 calculated for both models with drugs at various drug
# concentrations
panel3 = axs[1]

# Load APD data
APD_trapping_df = pd.read_csv(data_dir + 'SD_APD_fine.csv')
APD_conductance_df = pd.read_csv(data_dir + 'CS_APD_fine.csv')

# Identify EAD-like behaviour
drug_conc = APD_trapping_df['drug concentration'].values.tolist()
APD_trapping = APD_trapping_df['APD'].values.tolist()
APD_conductance = APD_conductance_df['APD'].values.tolist()
print(APD_trapping)
print(APD_conductance)
EAD_marker = [1050 if (np.isnan((i, j)).any()) else None for (i, j)
              in zip(APD_trapping, APD_conductance)]

# Plot APD90 of both models
panel3[0][0].plot(drug_conc, APD_trapping, 'o', color='orange',
                  label=SD_labelname)
panel3[0][0].plot(drug_conc, APD_conductance, '^', color='blue',
                  label=CS_labelname, alpha=0.8)
panel3[0][0].scatter(drug_conc, EAD_marker, marker=(5, 2),
                     color='k', label='EAD-like AP')
panel3[0][0].set_xscale("log", nonpositive='clip')
panel3[0][0].set_xlabel('Drug concentration (nM)')
panel3[0][0].set_ylabel(r'APD$_{90}$ (ms)')
panel3[0][0].legend(handlelength=1)

l_lim, r_lim = panel3[0][0].get_xlim()

# Bottom right panel
# Plots the qNet calculated for both models with drugs at a range of 0.5 to 25
# multiples of Cmax
panel4 = axs[2]

param_lib = modelling.BindingParameters()
Cmax = param_lib.binding_parameters[drug]['Cmax']

SD_qNet = APD_trapping_df['qNet'].values.tolist()
CS_qNet = APD_conductance_df['qNet'].values.tolist()

# EAD_indicator = np.array([1 if i is None else None for i in EAD_marker])
SD_qNet = [None if np.isnan(APD_trapping[i]) else SD_qNet[i] for i in
           range(len(SD_qNet))]
CS_qNet = [None if np.isnan(APD_conductance[i]) else CS_qNet[i] for i in
           range(len(CS_qNet))]

# Plot APD90 of both models
panel4[0][0].plot(drug_conc, SD_qNet, 'o', color='orange',
                  label=SD_labelname)
panel4[0][0].plot(drug_conc, CS_qNet, '^', color='blue',
                  label=CS_labelname, alpha=0.8)
panel4[0][0].set_xlabel("Drug concentration (nM)")
panel4[0][0].set_xscale("log", nonpositive='clip')
panel4[0][0].set_xlim(left=l_lim, right=r_lim)
panel4[0][0].set_ylabel('qNet (C/F)')
panel4[0][0].legend(handlelength=1)

# Add panel letter
# fig.fig.set_size_inches(10, 5.5)
fig.fig.text(0.1, 0.905, '(A)', fontsize=11)
fig.fig.text(0.1, 0.455, '(B)', fontsize=11)
# fig.fig.text(0.49, 0.905, '(B)', fontsize=11)
fig.fig.text(0.53, 0.455, '(C)', fontsize=11)

fig.savefig(fig_dir + "model_compare.pdf")
