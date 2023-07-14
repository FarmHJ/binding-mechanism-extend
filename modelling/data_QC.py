import numpy as np
import pandas as pd

import modelling

class QualityControl(object):
    """
    Methods to remove unwanted experimental traces:
    General
    1. Check seal resistance, membrane capacitance and series resistance
    For Pharm protocol
    1. Traces look like the protocol - 2nd segment higher than 3rd segment
    (control)
    """
    QC_names = ['const_range.Rseal', 'const_range.Cm', 'const_range.Rseries']

    def __init__(self):
        super(QualityControl, self).__init__()

        # Define thresholds
        # Range of experimental constants
        self.Rseal_thres = [1e8, 1e12]
        self.Cm_thres = [1e-12, 1e-10]
        # self.Rseries_thres = [1e6, 2.5e7]
        self.Rseries_thres = [5e6, 2e7]

        # Stability of traces
        self.rmsd0c = 0.2

        # Stability of experimental constants
        self.Rseals_control_stab = 0.3
        self.Cms_control_stab = 0.3
        self.Rseriess_control_stab = 0.3

        self.Rseals_stab = 0.5
        self.Cms_stab = 0.5
        self.Rseriess_stab = 0.5

        self.drug_washin = 10  # pulses

        self.thresholds = {
            self.QC_names[0]: [self.Rseal_thres, r"$\Omega$"],
            self.QC_names[1]: [self.Cm_thres, 'F'],
            self.QC_names[2]: [self.Rseries_thres, r"$\Omega$"]
        }

    def threshold_used(self, qc_name):
        return self.thresholds[qc_name]

    def const_range(self, Rseal, Cm, Rseries):
        if any(np.array(Rseal) < self.Rseal_thres[0]) or \
                any(np.array(Rseal) > self.Rseal_thres[1]):
            # print('Rseal: ', Rseal)
            qc_Rseal = False
        else:
            qc_Rseal = True

        if any(np.array(Cm) < self.Cm_thres[0]) or \
                any(np.array(Cm) > self.Cm_thres[1]):
            # print('Cm: ', Cm)
            qc_Cm = False
        else:
            qc_Cm = True

        if any(np.array(Rseries) < self.Rseries_thres[0]) or \
                any(np.array(Rseries) > self.Rseries_thres[1]):
            # print('Rseries: ', Rseries)
            qc_Rseries = False
        else:
            qc_Rseries = True

        return [qc_Rseal, qc_Cm, qc_Rseries]

    def signal_stable(self, trace1, trace2):
        rmsd0_1 = np.sqrt(np.mean((trace1) ** 2))
        rmsd0_2 = np.sqrt(np.mean((trace2) ** 2))
        rmsdc = np.mean([rmsd0_1, rmsd0_2]) * self.rmsd0c

        rmsd_trace = np.sqrt(np.mean((trace1 - trace2) ** 2))
        if rmsd_trace > rmsdc or not \
                (np.isfinite(rmsd_trace) and np.isfinite(rmsdc)):
            print('rmsd: ', rmsd_trace)
            print('rmsd to zero: ', rmsdc)
            return False
        return True

    def const_stable(self, Rseals, Cms, Rseriess, compound_change_t):

        n_pulses = len(Rseals)
        # Rseal
        Rseals_control = Rseals[:compound_change_t[0]]
        if np.std(Rseals_control) / np.mean(Rseals_control) > \
                self.Rseals_control_stab or not \
                (np.isfinite(np.mean(Rseals_control)) and
                 np.isfinite(np.std(Rseals_control))):
            qc_Rs_control = False
        else:
            qc_Rs_control = True

        compound_change_start = compound_change_t[0]
        qc_Rs = True
        for i in compound_change_t[1:] + [n_pulses]:
            Rseals_comp = Rseals[compound_change_start + self.drug_washin:i]
            if np.std(Rseals_comp) / np.mean(Rseals_comp) > \
                self.Rseals_stab or not \
                (np.isfinite(np.mean(Rseals_comp)) and
                 np.isfinite(np.std(Rseals_comp))):
                print('check')
                qc_Rs = False
            else:
                qc_Rs = qc_Rs and True
            compound_change_start = i

        # Exp constants before and after change of compounds should be stable
        qc_Rs = qc_Rs and qc_Rs_control

        # Cm
        Cms_control = Cms[:compound_change_t[0]]
        if np.std(Cms_control) / np.mean(Cms_control) > \
                self.Cms_control_stab or not \
                (np.isfinite(np.mean(Cms_control)) and
                 np.isfinite(np.std(Cms_control))):
            qc_Cm_control = False
        else:
            qc_Cm_control = True

        compound_change_start = compound_change_t[0]
        qc_Cm = True
        for i in compound_change_t[1:] + [n_pulses]:
            Cms_comp = Cms[compound_change_start + self.drug_washin:i]
            if np.std(Cms_comp) / np.mean(Cms_comp) > self.Cms_stab or not \
                (np.isfinite(np.mean(Cms_comp)) and
                 np.isfinite(np.std(Cms_comp))):
                qc_Cm = False
            else:
                qc_Cm = qc_Cm and True
            compound_change_start = i

        # Exp constants before and after change of compounds should be stable
        qc_Cm = qc_Cm and qc_Cm_control

        # Rseries
        Rseriess_control = Rseriess[:compound_change_t[0]]
        if np.std(Rseriess_control) / np.mean(Rseriess_control) > \
                self.Rseriess_control_stab or not \
                (np.isfinite(np.mean(Rseriess_control)) and
                 np.isfinite(np.std(Rseriess_control))):
            qc_Rseries_control = False
        else:
            qc_Rseries_control = True

        compound_change_start = compound_change_t[0]
        qc_Rseries = True
        for i in compound_change_t[1:] + [n_pulses]:
            Rseriess_comp = Rseriess[compound_change_start + self.drug_washin:
                                     i]
            if np.std(Rseriess_comp) / np.mean(Rseriess_comp) > \
                self.Rseriess_stab or not \
                (np.isfinite(np.mean(Rseriess_comp)) and
                 np.isfinite(np.std(Rseriess_comp))):
                qc_Rseries = False
            else:
                qc_Rseries = qc_Rseries and True
            compound_change_start = i

        # Exp constants before and after change of compounds should be stable
        qc_Rseries = qc_Rseries and qc_Rseries_control

        return [qc_Rs, qc_Cm, qc_Rseries]

    # def positive_block_peak_m1(self, trace_before, trace_after, win=None):


def plot_qc_trace(removed_cells, qc_name, fig_title=None):
    '''
    Plot traces of cells that did not satisfy the QC criterion
    '''
    dataset = modelling.DatasetLibrary()

    fig_each_row = []
    for prot in removed_cells.keys():
        for drug in removed_cells[prot].keys():
            cell_count = len(removed_cells[prot][drug])
            if cell_count > 0:
                fig_each_row.append(cell_count)

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
    cell_count = 0
    colors = ['blue', 'orange', 'green']

    for prot in removed_cells.keys():
        for drug in removed_cells[prot].keys():
            cell_index = 0
            cell_list = dataset.exp_data_list(prot, drug)
            detail_list = dataset.detail_read(prot, drug)

            for cell in removed_cells[prot][drug]:
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
                                time, signal / 1e-9,
                                color=colors[compound_count],
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
                cell_count += 1

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

            axs[len(fig_each_row)][j].set_xlabel('Time (ms)')

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

    thres, units = modelling.QualityControl.threshold_used(qc_name)
    min_thres_str = "{0:.0e}".format(thres[0])
    # min_thres_str = "{:.0f}".format(min_thres / 1e6)
    min_base, min_power = min_thres_str.split("e")
    min_power = int(min_power)
    min_str = min_base + r"$\times 10^{{{:d}}} $".format(min_power)
    # max_thres_str = "{:.0f}".format(max_thres / 1e6)
    max_thres_str = "{0:.0e}".format(thres[1])
    max_base, max_power = max_thres_str.split("e")
    max_power = int(max_power)
    max_str = max_base + r"$\times 10^{{{:d}}} $".format(max_power)
    if fig_title is None:
        fig.fig.suptitle(qc_name + ': [' + min_str + ', ' +
                         max_str + ']' + units)
    else:
        fig.fig.suptitle(fig_title)
