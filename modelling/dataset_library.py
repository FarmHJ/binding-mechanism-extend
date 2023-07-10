import myokit
import os
import pandas as pd


class DatasetLibrary(object):
    """
    A data library class that reads the experimental data
    for different protocols and different drugs
    """
    def __init__(self):
        super(DatasetLibrary, self).__init__()

        self._directory = os.path.join(
            # os.path.dirname(os.path.dirname(os.path.dirname(
            #     os.path.dirname(os.path.dirname(os.path.dirname(
            #         __file__)))))),
            # "home/scratch/220122_exp_data")
            os.path.dirname(os.path.dirname(__file__)), "exp_data")

        self.protocol_list = ["CIPA", "Pharm"]
        self.drug_list = ["cisapride", "dofetilide", "verapamil"]
        self.protocol_title = {"CIPA": "CiPA protocol",
                               "Pharm": "Roche's protocol"}
        self.compound_name = {
            "cisapride": ["19", "Cisapride"],
            "dofetilide": ["110", "RO0319253-000-001"],
            "verapamil": ["13", "Verapamil"]}

        self.compound_function = {
            "NMDG60 0.2%DMSO": "reference",
            "NMDG60 1uM E4031 0.2%DMSO": "positive block",
            "19": "cisapride", "13": "verapamil", "110": 'dofetilide',
            "S Solution": "reference", "Cisapride": "cisapride",
            "Verapamil": "verapamil", "RO0319253-000-001": 'dofetilide'}

    def exp_data_list(self, protocol, drug):
        """
        Returns file directories for different cells of given
        protocol and drug.
        """

        # check protocol choice
        if protocol not in self.protocol_list:
            raise ValueError(
                "Choice of protocol must be either 'CIPA' or 'Pharm'")

        if drug not in self.drug_list:
            raise ValueError(
                "Choice of drug must be one of 'cisapride', 'dofetilide' \
                    and 'verapamil'")

        filepath = os.path.join(self._directory, protocol, drug)
        file_cells = os.listdir(filepath)

        # Find name of cells and their file path
        files = []
        cells = []
        for filename in file_cells:
            if filename.endswith(".csv") and not filename.startswith("OA"):
                cell_name = filename[-8:-5]
                cells.append(cell_name)
                files.append(os.path.join(filepath, filename))
            if filename.startswith("DMSO"):
                conc_file = os.path.join(filepath, filename)
                conc_info = pd.read_excel(conc_file)

        cell_name_path = pd.DataFrame(data={
            "cells": cells, "file_path": files,
            "drug_concentration": [0] * len(cells)})

        # Get drug concentration used for each cell
        for cell_num in cells:
            conc = conc_info.loc[conc_info["WELL ID"] == cell_num,
                                 "Concentration (Mol/L)"].iloc[0]
            cell_name_path.loc[cell_name_path["cells"] == cell_num,
                               "drug_concentration"] = conc

        return cell_name_path

    def exp_data_read(self, filepath):
        """
        Returns the dataframe of experimental data
        of given file path
        """
        data = pd.read_csv(filepath, header=2, sep=';')

        return data

    def detail_read(self, protocol, drug):
        """
        Returns parameters of the experiment.
        Including added compound name, its concentration,
        capacitance, seal resistance and series resistance.
        """
        # check protocol choice
        if protocol not in self.protocol_list:
            raise ValueError(
                "Choice of protocol must be either 'CIPA' or 'Pharm'")

        if drug not in self.drug_list:
            raise ValueError(
                "Choice of drug must be one of 'Cisapride', 'Dofetilide' \
                    and 'Verapamil'")

        filepath = os.path.join(self._directory, protocol, drug)
        file_cells = os.listdir(filepath)
        for filename in file_cells:
            if filename.endswith("_processed.csv"):
                file_dir = os.path.join(filepath, filename)
                data = pd.read_csv(file_dir, index_col=0)

        if data.index[0] == "Well ID":
            data = data.rename(index={"Well ID": "Parameter"})

        return data

    def cell_detail(self, detail_list, cell):
        """
        Returns parameters of the experiment for a specific cell.
        Including added compound name, its concentration,
        capacitance, seal resistance and series resistance.
        """
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

        return detail


class ProtocolLibrary(object):
    """
    A library class with known protocols
    """
    def __init__(self):
        super(ProtocolLibrary, self).__init__()

    def Milnes(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 800, period=t_max)
        protocol.schedule(-90, 800, 100, period=t_max)
        protocol.schedule(-80, 900, 100, period=t_max)
        protocol.schedule(-80, 11000, 14000, period=t_max)

        return protocol

    def current_impulse(self, t_max, offset=50):
        return myokit.pacing.blocktrain(t_max, 1, offset=offset)

    def validation3(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100)  # , period=t_max)
        protocol.schedule(-40, 100, 50)  # , period=t_max)
        protocol.schedule(20, 150, 500)  # , period=t_max)
        protocol.schedule(-40, 650, 500)  # , period=t_max)
        protocol.schedule(-80, 1150, 200)  # , period=t_max)

        return protocol

    def hERG_validation(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100)
        protocol.schedule(20, 100, 900)
        protocol.schedule(-80, 1000, 1000)

        return protocol

    def Pneg80(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 200, period=t_max)
        protocol.schedule(20, 200, 500, period=t_max)
        protocol.schedule(-50, 700, 200, period=t_max)
        protocol.schedule(-80, 900, 4500, period=t_max)

        return protocol

    def P0(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100, period=t_max)
        protocol.schedule(-60, 5100, 200, period=t_max)
        protocol.schedule(-80, 5300, 100, period=t_max)

        return protocol

    def P40(self, t_max):
        protocol = myokit.Protocol()
        protocol.schedule(-80, 0, 100, period=t_max)
        protocol.schedule(40, 100, 5000, period=t_max)
        protocol.schedule(-60, 5100, 200, period=t_max)
        protocol.schedule(-80, 5300, 100, period=t_max)

        return protocol
