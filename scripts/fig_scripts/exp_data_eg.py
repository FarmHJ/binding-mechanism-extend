import pandas as pd

import modelling


dataset = modelling.DatasetLibrary()

drug_list = dataset.drug_list
protocol_list = dataset.protocol_list

protocol_title = {"CIPA": "CiPA protocol", "Pharm": "Roche's protocol"}

fig = modelling.figures.FigureStructure(
    figsize=(8, 5), gridspec=(len(protocol_list), len(drug_list)),
    hspace=0.3, height_ratios=[1] * len(protocol_list), wspace=0.15)
plot = modelling.figures.FigurePlot()

count = 0
colors = ['blue', 'orange', 'green']

chosen_cells = {
    'CIPA': {
        'dofetilide': 'J07',
        'verapamil': 'A21',
        'cisapride': 'F19'
    },
    'Pharm': {
        'dofetilide': 'J19',
        'verapamil': 'N17',
        'cisapride': 'E14'
    }
}

chosen_drug_conc = {
    'CIPA': {
        'dofetilide': 0,
        'verapamil': 0,
        'cisapride': 0
    },
    'Pharm': {
        'dofetilide': 0,
        'verapamil': 0,
        'cisapride': 0
    }
}

for protocol_count, protocol in enumerate(protocol_list):

    y_lb, y_ub = 0, 0.5
    x_lb, x_ub = 10, 0

    for drug_count, drug in enumerate(drug_list):

        chosen_cell = chosen_cells[protocol][drug]
        cell_list = dataset.exp_data_list(protocol, drug)
        detail_list = dataset.detail_read(protocol, drug)
        if detail_list.index[0] == "Well ID":
            detail_list = detail_list.rename(index={"Well ID": "Parameter"})

        cell_file_path = cell_list.loc[cell_list["cells"] == chosen_cell][
            "file_path"].values[0]
        drug_conc = cell_list.loc[cell_list["cells"] == chosen_cell][
            "drug_concentration"].values[0]
        data = dataset.exp_data_read(cell_file_path)

        chosen_drug_conc[protocol][drug] = drug_conc

        # get drug concentration info
        test = detail_list.loc[["Parameter", chosen_cell], :]
        detail = test.T.reset_index().rename(columns={"index": "Sweep"})

        for i in range(len(detail.index)):
            sweep_string = detail.loc[i, "Sweep"]
            if len(sweep_string) == 9:
                detail.loc[i, "Sweep"] = int(sweep_string[-3:])
            else:
                detail.loc[i, "Sweep"] = int(sweep_string[-5:-2])

        detail = detail.rename(columns={chosen_cell: "values",
                                        "Parameter": chosen_cell})
        detail = detail.pivot(index='Sweep', columns=chosen_cell,
                              values='values')

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
                    compound_label = dataset.compound_function[
                        compound_names[compound_count]]
                    label = compound_label if compound_label in \
                        ['reference', 'positive block'] else 'drug'
                    fig.axs[protocol_count][drug_count].plot(
                        time, signal / 1e-9, color=colors[compound_count],
                        label=label, alpha=0.5, zorder=-1)
                previous_compound = compound
                pulse_count += 1
            elif compound != previous_compound:
                pulse_count = 1
                previous_compound = compound
                compound_change.append(sweep + 1)
                compound_count += 1

        y_bottom, y_top = fig.axs[protocol_count][drug_count].get_ylim()
        if y_bottom < y_lb:
            y_lb = y_bottom
        if y_top > y_ub:
            y_ub = y_top
        x_left, x_right = fig.axs[protocol_count][drug_count].get_xlim()
        if x_left < x_lb:
            x_lb = x_left
        if x_right > x_ub:
            x_ub = x_right

    for d, drug in enumerate(drug_list):
        fig.axs[protocol_count][d].set_ylim(y_lb, y_ub)
        fig.axs[protocol_count][d].set_rasterization_zorder(0)

        drug_conc_str = "{0:.0e}".format(chosen_drug_conc[protocol][drug] /
                                         1e-6)
        base, power = drug_conc_str.split("e")
        power = int(power)
        fig.axs[protocol_count][d].text(
            10, 0.95 * y_ub,
            base + r"$\times 10^{{{:d}}} \mu$".format(power) + 'M',
            fontsize=8, ha='left', va='top')

    # # get stimulus
    # stimulus = data["Stimulus"] * 1000
    # free_panel_row = min([i for i in range(len(unique_drug_concs)) if
    #                         cell_counts[unique_drug_concs[i]] <
    #                         max_cell_per_conc])
    # free_panel_col = cell_counts[unique_drug_concs[free_panel_row]]
    # protocol_axs = fig.fig.add_subplot(
    #     fig.gs[free_panel_row, free_panel_col])
    # protocol_axs.plot(time, stimulus, 'k')

    # protocol_axs.yaxis.tick_right()
    # protocol_axs.yaxis.set_label_position('right')
    # protocol_axs.set_ylabel('Voltage (mV)')
    # protocol_axs.sharex(axs[0][0])
    # protocol_axs.tick_params(labelbottom=False, labelleft=False)
    # protocol_axs.spines['top'].set_visible(False)
    # protocol_axs.spines['left'].set_visible(False)

for d in range(len(drug_list)):
    for i in range(len(protocol_list)):
        # Share x-axis with first column
        if d != 0:
            fig.axs[i][d].sharex(fig.axs[i][0])
            fig.axs[i][d].sharey(fig.axs[i][0])

        fig.axs[i][d].set_xlabel('Time (ms)')

        # Label y-axis at the first column
        if d != 0:
            fig.axs[i][d].tick_params(labelleft=False)
        else:
            fig.axs[i][d].set_ylabel('Current (nA)')

    fig.axs[0][d].set_title(drug_list[d])

unique_label = fig.legend_without_duplicate_labels(fig.axs[0][0])
fig.axs[0][0].legend(*zip(*unique_label), handlelength=1, loc=(0.03, 0.5))
fig.fig.text(0.1, 0.895, '(A)', fontsize=11)
fig.fig.text(0.1, 0.455, '(B)', fontsize=11)

filename = "exp_data_eg.svg"
fig.savefig("../../figures/experimental_data/" + filename)
