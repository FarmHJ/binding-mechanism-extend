"""
Constant variables: seal resistance, membrane capacitance, series resistance
Variable names: Rseal, Cm, Rseries
Column name in data: 'Seal Resistance', 'Capacitance', 'Series Resistance'
Defined threshold: Rseal_thres, Cm_thres, Rseries_thres
y-axis boundaries: [1e10, -1], [1e-11, 1e-13], [1e7, 1e5]
"""

import numpy as np
import pandas as pd

import modelling

dataset = modelling.DatasetLibrary()

drug_list = dataset.drug_list
protocol_list = dataset.protocol_list

protocol_title = {"CIPA": "CiPA protocol", "Pharm": "Roche's protocol"}

count = 0
colors = ['blue', 'orange', 'green']
for drug_count, drug in enumerate(drug_list):
    protocol_count = 0
    for protocol in protocol_list:
        cell_list = dataset.exp_data_list(protocol, drug)
        cell_list = cell_list.sort_values('drug_concentration')
        detail_list = dataset.detail_read(protocol, drug)
        if detail_list.index[0] == "Well ID":
            detail_list = detail_list.rename(index={"Well ID": "Parameter"})

        drug_concs_list = cell_list["drug_concentration"].values
        unique_drug_concs = pd.unique(drug_concs_list)
        cell_counts = pd.Series(drug_concs_list).value_counts()
        max_cell_per_conc = max(cell_counts.values)

        # Set up structure of the figure
        gridspec = (len(unique_drug_concs), max_cell_per_conc)
        fig = modelling.figures.FigureStructure(
            figsize=(2 * max_cell_per_conc, 2 * len(unique_drug_concs)),
            gridspec=gridspec, hspace=0.2, wspace=0.25,
            height_ratios=[1] * len(unique_drug_concs), plot_in_subgrid=True)
        plot = modelling.figures.FigurePlot()

        axs = [[fig.fig.add_subplot(fig.gs[i, j]) for j in
                range(cell_counts[unique_drug_concs[i]])] for i in
               range(len(unique_drug_concs))]

        cell_count = 0
        y_lb, y_ub = 1e7, 1e5
        x_lb, x_ub = 10, 0

        for cell in cell_list['cells'].values:
            cell_file_path = cell_list.loc[cell_list["cells"] == cell][
                "file_path"].values[0]
            drug_conc = cell_list.loc[cell_list["cells"] == cell][
                "drug_concentration"].values[0]
            # data = dataset.exp_data_read(cell_file_path)

            # get drug concentration info
            test = detail_list.loc[["Parameter", cell], :]
            detail = test.T.reset_index().rename(columns={"index": "Sweep"})

            for i in range(len(detail.index)):
                sweep_string = detail.loc[i, "Sweep"]
                if len(sweep_string) == 9:
                    detail.loc[i, "Sweep"] = int(sweep_string[-3:])
                else:
                    detail.loc[i, "Sweep"] = int(sweep_string[-5:-2])

            detail = detail.rename(columns={cell: "values",
                                            "Parameter": cell})
            detail = detail.pivot(index='Sweep', columns=cell,
                                  values='values')

            compound_names = detail["Compound Name"].values.ravel()
            compound_names = pd.unique(compound_names)

            # time = data["Sample Time (us)"] / 1000
            drug_conc_index = np.where(unique_drug_concs == drug_conc)[0][0]
            cell_num_index = np.count_nonzero(
                drug_concs_list[:cell_count] == drug_conc)

            sweep_count = detail.shape[0]
            constants = [float(i) for i in
                         detail.loc[:, 'Series Resistance'].values]
            constants_thres = modelling.QualityControl().Rseries_thres

            compound_count = 0
            previous_compound = detail.loc[1, "Compound Name"]
            compound_change = []

            for sweep in range(sweep_count):
                compound = detail.loc[sweep + 1, "Compound Name"]
                if compound != previous_compound and not \
                        isinstance(compound, float):
                    compound_change.append(sweep + 1)
                    compound_count += 1
                previous_compound = compound

            axs[drug_conc_index][cell_num_index].plot(np.arange(sweep_count),
                                                      constants)
            constants_thres = modelling.QualityControl().Rseries_thres
            axs[drug_conc_index][cell_num_index].axhline(constants_thres[0], color='k')
            axs[drug_conc_index][cell_num_index].axhline(constants_thres[1], color='k')

            for c in range(compound_count):
                axs[drug_conc_index][cell_num_index].axvline(
                    compound_change[c], color='red')

            axs[drug_conc_index][cell_num_index].set_yscale('log', nonpositive='clip')
            y_bottom, y_top = axs[drug_conc_index][cell_num_index].get_ylim()
            # print("{:.2E}".format(y_bottom))
            if y_bottom < y_lb:
                y_lb = y_bottom
            if y_top > y_ub:
                y_ub = y_top
            x_left, x_right = axs[drug_conc_index][cell_num_index].get_xlim()
            if x_left < x_lb:
                x_lb = x_left
            if x_right > x_ub:
                x_ub = x_right
            cell_count += 1

        # print("{:.2E}".format(y_lb), "{:.2E}".format(y_ub))
        for i in range(len(unique_drug_concs)):
            for j in range(cell_counts[unique_drug_concs[i]]):
                axs[i][j].set_ylim(y_lb, y_ub)
                axs[i][j].set_xlim(x_lb, x_ub)
                axs[i][j].set_rasterization_zorder(0)

                # Label x-axis at the last row
                if i != len(unique_drug_concs) - 1:
                    axs[i][j].tick_params(labelbottom=False)
                else:
                    axs[i][j].set_xlabel('Sweep')

                # Label y-axis at the first column
                if j != 0:
                    axs[i][j].tick_params(labelleft=False)
                else:
                    axs[i][j].set_ylabel('Resistance (Ohm)')

            title_str_num = "{0:.0e}".format(unique_drug_concs[i] / 1e-6)
            base, power = title_str_num.split("e")
            power = int(power)
            axs[i][0].text(
                2, 0.8 * y_ub,
                base + r"$\times 10^{{{:d}}} \mu$".format(power) + 'M',
                fontsize=8, ha='left', va='top')

        # # get stimulus
        # stimulus = data["Stimulus"] * 1000
        # free_panel_row = min([i for i in range(len(unique_drug_concs)) if
        #                       cell_counts[unique_drug_concs[i]] <
        #                       max_cell_per_conc])
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

        filename = protocol + "_" + drug + "_Rseries.pdf"
        fig.savefig("../figures/experimental_data/" + filename)
