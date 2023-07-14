import numpy as np
import pandas as pd

import modelling

dataset = modelling.DatasetLibrary()
trace_qc = modelling.QualityControl()

drug_list = dataset.drug_list
protocol_list = dataset.protocol_list

colors = ['blue', 'orange', 'green']

removed_cells = {
    'CIPA': {
        'dofetilide': [],
        'verapamil': [],
        'cisapride': []
    },
    'Pharm': {
        'dofetilide': [],
        'verapamil': [],
        'cisapride': []
    }
}

for protocol_count, protocol in enumerate(protocol_list):
    for drug_count, drug in enumerate(drug_list):
        cell_list = dataset.exp_data_list(protocol, drug)
        cell_list = cell_list.sort_values('drug_concentration')
        detail_list = dataset.detail_read(protocol, drug)

        for cell in cell_list['cells'].values:
            detail = dataset.cell_detail(detail_list, cell)
            cell_file_path = cell_list.loc[cell_list["cells"] == cell][
                "file_path"].values[0]
            data = dataset.exp_data_read(cell_file_path)

            # QC to choose data with experimental constants within threshold
            Rseal = [float(i) for i in detail.loc[:, 'Seal Resistance'].values]
            Rseries = [float(i) for i in
                       detail.loc[:, 'Series Resistance'].values]
            Cm = [float(i) for i in detail.loc[:, 'Capacitance'].values]
            QC_constants = trace_qc.qc_general(Rseal, Cm, Rseries)
            # if not QC_constants[2]:
            # if all(QC_constants):
            #     removed_cells[protocol][drug].append(cell)

            total_pulses = data.shape[1] - 4
            compound_change = []
            previous_compound = detail.loc[1, "Compound Name"]
            for sweep in range(total_pulses):
                signal = data.iloc[:, sweep + 3]
                compound = detail.loc[sweep + 1, "Compound Name"]
                if compound != previous_compound:
                    previous_compound = compound
                    compound_change.append(sweep + 1)

            # QC to make sure the traces are stable at the end of the block
            signal_stable = []
            for sweep in compound_change:
                valid_trace_ind = 1
                chosen_trace = []
                nan_trace = 2
                while nan_trace > 0:
                    temp_trace = data.iloc[:, sweep + 3 - valid_trace_ind]
                    if not any(np.isnan(np.array(temp_trace))):
                        chosen_trace.append(temp_trace)
                        nan_trace -= 1
                    valid_trace_ind += 1
                trace1 = chosen_trace[0] / 1e-9
                trace2 = chosen_trace[1] / 1e-9
                compound = detail.loc[sweep + 1 - 1, "Compound Name"]
                QC_stable = trace_qc.qc_stable(trace1, trace2)
                signal_stable.append(QC_stable)
            print(signal_stable)
            if all(signal_stable) and all(QC_constants):
                removed_cells[protocol][drug].append(cell)

fig_each_row = []
for prot in removed_cells.keys():
    for drug in removed_cells[prot].keys():
        cell_count = len(removed_cells[prot][drug])
        if cell_count > 0:
            fig_each_row.append(cell_count)

print('ready for figure')
gridspec = (len(fig_each_row), max(fig_each_row))
fig = modelling.figures.FigureStructure(
    figsize=(2 * max(fig_each_row), 2 * len(fig_each_row)),
    gridspec=gridspec, hspace=0.2, wspace=0.25,
    height_ratios=[1] * len(fig_each_row), plot_in_subgrid=True)
axs = [[fig.fig.add_subplot(fig.gs[i, j]) for j in
        range(fig_each_row[i])] for i in range(len(fig_each_row))]

row_index = 0
y_lb, y_ub = 0, 0.5
x_lb, x_ub = 10, 0
for prot in removed_cells.keys():
    for drug in removed_cells[prot].keys():

        cell_index = 0
        for cell in removed_cells[prot][drug]:
            cell_list = dataset.exp_data_list(prot, drug)
            detail_list = dataset.detail_read(prot, drug)
            detail = dataset.cell_detail(detail_list, cell)

            cell_file_path = cell_list.loc[cell_list["cells"] == cell][
                "file_path"].values[0]
            drug_conc = cell_list.loc[cell_list["cells"] == cell][
                "drug_concentration"].values[0]
            data = dataset.exp_data_read(cell_file_path)

            compound_names = detail["Compound Name"].values.ravel()
            compound_names = pd.unique(compound_names)

            time = data["Sample Time (us)"] / 1000
            compound_count = 0
            pulse_count = 0
            previous_compound = detail.loc[1, "Compound Name"]
            total_pulses = data.shape[1] - 4
            compound_change = []

            for sweep in range(total_pulses):
                signal = data.iloc[:, sweep + 3]
                compound = detail.loc[sweep + 1, "Compound Name"]
                if compound == previous_compound:
                    if pulse_count >= 10:
                        axs[row_index][cell_index].plot(
                            time, signal / 1e-9, color=colors[compound_count],
                            label=dataset.compound_function[
                                compound_names[compound_count]],
                            alpha=0.5, zorder=-1)
                    previous_compound = compound
                    pulse_count += 1
                elif compound != previous_compound:
                    pulse_count = 1
                    previous_compound = compound
                    compound_change.append(sweep + 1)
                    compound_count += 1

            drug_conc_str = "{0:.0e}".format(drug_conc / 1e-6)
            base, power = drug_conc_str.split("e")
            power = int(power)
            axs[row_index][cell_index].text(
                0.05, 0.95,
                base + r"$\times 10^{{{:d}}} \mu$".format(power) + 'M',
                fontsize=8, ha='left', va='top',
                transform=axs[row_index][cell_index].transAxes)

            y_bottom, y_top = axs[row_index][cell_index].get_ylim()
            if y_bottom < y_lb:
                y_lb = y_bottom
            if y_top > y_ub:
                y_ub = y_top
            x_left, x_right = axs[row_index][cell_index].get_xlim()
            if x_left < x_lb:
                x_lb = x_left
            if x_right > x_ub:
                x_ub = x_right
            cell_index += 1

        # if axs[row_index]:
        if removed_cells[prot][drug]:
            axs[row_index][0].text(0.05, 0.85, drug, fontsize=8,
                                   ha='left', va='top',
                                   transform=axs[row_index][0].transAxes)
            row_index += 1

for i in range(len(fig_each_row)):
    for j in range(fig_each_row[i]):
        axs[i][j].set_ylim(y_lb, y_ub)
        axs[i][j].set_rasterization_zorder(0)

        # Share x-axis with first column
        if j != 0:
            axs[i][j].sharex(axs[i][0])
            axs[i][j].sharey(axs[i][0])

        # # Label x-axis at the last row
        # if i != len(fig_each_row) - 1:
        #     axs[i][j].tick_params(labelbottom=False)
        # else:
        #     axs[i][j].set_xlabel('Time (ms)')

        # Label y-axis at the first column
        if j != 0:
            axs[i][j].tick_params(labelleft=False)
        else:
            axs[i][j].set_ylabel('Current (nA)')

free_panel_row = min([i for i in range(len(fig_each_row)) if
                      fig_each_row[i] < max(fig_each_row) and
                      fig_each_row[i] != 0])
free_panel_col = fig_each_row[free_panel_row]
legend_axs = fig.fig.add_subplot(
    fig.gs[free_panel_row, free_panel_col])
legend_axs.xaxis.set_visible(False)
legend_axs.yaxis.set_visible(False)
legend_axs.set_frame_on(False)
unique_label = fig.legend_without_duplicate_labels(axs[0][0])
legend_axs.legend(*zip(*unique_label), handlelength=1, loc='upper left')

min_thres, max_thres = trace_qc.Rseries_thres
# min_thres_str = "{0:.0e}".format(min_thres / 1e-12)
min_thres_str = "{:.0f}".format(min_thres / 1e6)
# min_base, min_power = min_thres_str.split("e")
# min_power = int(min_power)
# base + r"$\times 10^{{{:d}}} $".format(power) + 'pF'
max_thres_str = "{:.0f}".format(max_thres / 1e6)
# max_base, max_power = max_thres_str.split("e")
# fig.fig.suptitle('membrane capacitance: [' + str(int(min_thres)) + ', ' +
#                  str(int(max_thres)) + '] pF')
# fig.fig.suptitle('series resistance: [' + min_thres_str + ', ' +
#                  max_thres_str + '] ' + r"$M\Omega$")
fig.fig.suptitle('chosen_cells')

filename = "chosen_cells.pdf"
fig.savefig("../../figures/experimental_data/" + filename)
