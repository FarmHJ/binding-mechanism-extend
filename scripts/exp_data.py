import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

import modelling

dataset = modelling.DatasetLibrary()

drug_list = dataset.drug_list
protocol_list = dataset.protocol_list


fig = plt.figure(figsize=(6.5, 6))
gs = fig.add_gridspec(len(drug_list) + 1, len(protocol_list), height_ratios=[1] + ([3] * len(drug_list)),)
protocol_axs = [fig.add_subplot(gs[0, i]) for i in range(len(protocol_list))]
axs = [fig.add_subplot(gs[i + 1, j]) for i in range(len(drug_list)) for j in range(len(protocol_list))]

protocol_title = {"CIPA": "CiPA protocol", "Pharm": "Roche's protocol"}
compound_name = {"19": "Cisapride", "13": "Verapamil",
                 "Cisapride": "Cisapride", "Verapamil": "Verapamil"}
cell_choice = pd.DataFrame(data={"Cisapride": ["K22", "F17"], "Verapamil": ["H18", "N17"]})
cell_choice = cell_choice.rename(index={0: "CIPA", 1: "Pharm"})
drug_concentration = pd.DataFrame(data={"Cisapride": ["0.3uM", "1um"], "Verapamil": ["0.3uM", "0.3uM"]})
drug_concentration = drug_concentration.rename(index={0: "CIPA", 1: "Pharm"})

axs[0].sharex(axs[2])
axs[1].sharex(axs[3])

count = 0
for drug in drug_list:
    protocol_count = 0
    for protocol in protocol_list:
        cell_list = dataset.exp_data_list(protocol, drug)
        detail_list = dataset.detail_read(protocol, drug)
        if detail_list.index[0] == "Well ID":
            detail_list = detail_list.rename(index={"Well ID": "Parameter"})

        cell = cell_choice.loc[protocol, drug]
        cell_file_path = cell_list.loc[cell_list["cells"] == cell]["file_path"].reset_index(drop=True)[0]
        drug_conc = drug_concentration.loc[protocol, drug]
#         for cell in [cell_list["file_path"][0]]:
        for cell in [cell_file_path]:
            data = dataset.exp_data_read(cell)
            cell_name = cell_list.loc[cell_list["file_path"] == cell]["cells"].reset_index(drop=True)[0]

#             cell_name = cell[-8:-5]
            
            # get drug concentration info
            test = detail_list.loc[["Parameter", cell_name], :]
            detail = test.T.reset_index().rename(columns={"index" : "Sweep"})

            for i in range(len(detail.index)):
                sweep_string = detail.loc[i, "Sweep"]
                if len(sweep_string) == 9:
                    detail.loc[i, "Sweep"] = int(sweep_string[-3:])
                else:
                    detail.loc[i, "Sweep"] = int(sweep_string[-5:-2])

            detail = detail.rename(columns={cell_name : "values", "Parameter": cell_name})
            detail = detail.pivot(index='Sweep', columns=cell_name, values='values')
            
            compound_names = detail["Compound Name"].values.ravel()
            compound_names =  pd.unique(compound_names)
            print(compound_names)

            # get stimulus
            stimulus = data["Stimulus"] * 1000
            time = data["Sample Time (us)"] / 1000
            protocol_axs[protocol_count].plot(time, stimulus)
#             protocol_axs[protocol_count].set_title(protocol_title[protocol])
            protocol_axs[protocol_count].set_ylabel('Voltage (mV)')
            protocol_axs[protocol_count].set_xlabel('Time (ms)')
        
            compound_count1 = 0
            compound_count2 = 0
            for sweep in range(data.shape[1] - 4):
                signal = data.iloc[:, sweep + 3]
#                 conc = detail.loc[sweep + 1, "Concentration"]
                compound = detail.loc[sweep + 1, "Compound Name"]
                if compound == compound_names[0]:
                    compound_count1 += 1
                    label = 'before drug'
                    if compound_count1 >= 10:
                        axs[count].plot(time, signal / 1e-9, color='blue', label=label, alpha=0.5, zorder=-5)
                elif compound == compound_names[1]:
                    compound_count2 += 1
                    label = drug_conc + ' ' + compound_name[compound]
                    if compound_count2 >= 10:
                        axs[count].plot(time, signal / 1e-9, color='red', label=label, alpha=0.5, zorder=-5)
#                 else:
#                     label = 'after washout'
#                     axs[count].plot(time, signal, color='blue', label=compound)
#         axs[count].set_title(drug)
        axs[count].set_ylabel('Current (nA)')
        axs[count].set_xlabel('Time (ms)')
        axs[count].set_ylim(0, 1.1)
        axs[count].legend()
        axs[count].set_rasterization_zorder(0)
        legend_without_duplicate_labels(axs[count])
        count += 1
        protocol_count += 1
        
# axs[0].set_xticks([]) 
# axs[1].set_xticks([])
for ax in protocol_axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

fig.text(0.055, 0.975, '(a)', fontsize=10)
fig.text(0.555, 0.975, '(b)', fontsize=10)
fig.text(0.055, 0.805, '(c)', fontsize=10)
fig.text(0.555, 0.805, '(d)', fontsize=10)
fig.text(0.055, 0.405, '(e)', fontsize=10)
fig.text(0.555, 0.405, '(f)', fontsize=10)

# fig1.text(0.055, 0.975, '(a)', fontsize=10)
# fig1.text(0.055, 0.775, '(b)', fontsize=10)
# fig1.text(0.055, 0.425, '(c)', fontsize=10)

fig.set_tight_layout(True)
# plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("../figures/exp_data_visualisation.pdf", dpi=500)
plt.show()