import pandas as pd

import modelling


dataset = modelling.DatasetLibrary()

drug_list = dataset.drug_list
protocol_list = dataset.protocol_list

variable_list = ['seal resistance', 'cell capacitance', 'series resistance']

fig = modelling.figures.FigureStructure(
    figsize=(8, 2), gridspec=(1, 3), wspace=0.35)
plot = modelling.figures.FigurePlot()

count = 0
colors = ['blue', 'orange', 'green']

# Experimental data examples
chosen_cells = {
    'seal resistance': {
        'protocol': 'CIPA',
        'drug': 'verapamil',
        'cell': 'D20'
    },
    'cell capacitance': {
        'protocol': 'CIPA',
        'drug': 'cisapride',
        'cell': 'D19'
    },
    'series resistance': {
        'protocol': 'Pharm',
        'drug': 'cisapride',
        'cell': 'D14'
    },
}

y_lb, y_ub = 0, 0.5
for variable_count, variable in enumerate(variable_list):

    chosen_cell = chosen_cells[variable]['cell']
    protocol = chosen_cells[variable]['protocol']
    drug = chosen_cells[variable]['drug']
    cell_list = dataset.exp_data_list(protocol, drug)
    detail_list = dataset.detail_read(protocol, drug)
    if detail_list.index[0] == "Well ID":
        detail_list = detail_list.rename(index={"Well ID": "Parameter"})

    cell_file_path = cell_list.loc[cell_list["cells"] == chosen_cell][
        "file_path"].values[0]
    data = dataset.exp_data_read(cell_file_path)

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
                fig.axs[0][variable_count].plot(
                    time, signal / 1e-9, color=colors[compound_count],
                    label=label, alpha=0.5, zorder=-1)
            previous_compound = compound
            pulse_count += 1
        elif compound != previous_compound:
            pulse_count = 1
            previous_compound = compound
            compound_change.append(sweep + 1)
            compound_count += 1

    y_bottom, y_top = fig.axs[0][variable_count].get_ylim()
    if y_bottom < y_lb:
        y_lb = y_bottom
    if y_top > y_ub:
        y_ub = y_top

for v in range(len(variable_list)):
    fig.axs[0][v].set_ylim(y_lb, y_ub)
    fig.axs[0][v].set_rasterization_zorder(0)

    fig.axs[0][v].set_xlabel('Time (ms)')
    fig.axs[0][v].set_ylabel('Current (nA)')

    fig.axs[0][v].set_title(variable_list[v])

unique_label = fig.legend_without_duplicate_labels(fig.axs[0][1])
fig.axs[0][1].legend(*zip(*unique_label), handlelength=1)  # , loc=(0.03, 0.5))
# fig.fig.text(0.1, 0.895, '(A)', fontsize=11)
# fig.fig.text(0.1, 0.455, '(B)', fontsize=11)

filename = "QC_eg.svg"
fig.savefig("../../figures/experimental_data/" + filename)
