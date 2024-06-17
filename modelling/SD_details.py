# Reference:
# Li Z, Dutta S, Sheng J, Tran PN, Wu W, Chang K, Mdluli T, Strauss DG,
# Colatsky T. Improving the In Silico Assessment of Proarrhythmia Risk by
# Combining hERG (Human Ether-Ã -go-go-Related Gene) Channel-Drug Binding
# Kinetics and Multichannel Pharmacology. Circ Arrhythm Electrophysiol.
# 2017 Feb;10(2):e004628. doi: 10.1161/CIRCEP.116.004628.

import myokit
import numpy as np
import os
import pandas as pd

import modelling


class BindingParameters(object):
    """
    To create a library of all the dynamic hERG-binding parameters for
    different drug compounds.
    (ref. Table 2)
    """

    def __init__(self, ikr_model='Li'):
        super(BindingParameters, self).__init__()

        self.ikr_model = ikr_model

    def load_SD_parameters(self, drug=None):
        param_file = os.path.join(modelling.PARAM_DIR,
                                  self.ikr_model + '-SD.csv')
        self.binding_params = pd.read_csv(param_file, index_col=0)
        if 'error' in self.binding_params.columns:
            self.binding_params = self.binding_params.drop(columns=['error'])

        if drug is not None:
            if drug not in drug_names:
                NameError("Choose a drug compound within the \
                        list in `drug_compounds`")

            return self.binding_params.loc[[drug]]

    def get_SD_parameters(self, drug):
        if drug not in drug_names:
            NameError("Choose a drug compound within the \
                      list in `drug_compounds`")

        return self.binding_params.loc[[drug]]

    def load_published_Hill_eq(self, drug, ikr_model='Li', channel=None):
        if drug not in drug_names:
            NameError("Choose a drug compound within the \
                      list in `drug_compounds`")

        param_file = os.path.join(modelling.PARAM_DIR, ikr_model + '-Hill.csv')
        self.Hill = pd.read_csv(param_file, index_col=0, header=[0, 1])
        channel_list = set([header[0] for header in self.Hill.columns])
        # see if can reduce the loading of csv file everytime

        if channel is None:
            return self.Hill.loc[[drug]]
        else:
            if channel not in channel_list:
                NameError("Choose an ion channel within the \
                          list in `channels`")
            return self.Hill.loc[[drug]][channel].values.tolist()[0]

    def load_Hill_eq(self, drug, ikr_model='Li'):
        if drug not in drug_names:
            NameError("Choose a drug compound within the \
                      list in `drug_compounds`")

        param_file = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                                  'Hill_curves', drug, ikr_model + '_Hill.csv')
        self.Hill = pd.read_csv(param_file, index_col=[0],
                                skipinitialspace=True)

        return self.Hill.T.values.tolist()[0]


# SD_Parameters
drug_names = ['dofetilide', 'bepridil', 'terfenadine',
              'cisapride', 'verapamil', 'ranolazine',
              'quinidine', 'sotalol', 'chlorpromazine',
              'ondansetron', 'diltiazem', 'mexiletine']
validation_drugs = ['vandetanib', 'ibutilide', 'azimilide', 'disopyramide',
                    'domperidone', 'droperidol', 'pimozide', 'clozapine',
                    'risperidone', 'astemizole', 'clarithromycin', 'tamoxifen',
                    'metoprolol', 'loratadine', 'nitrendipine', 'nifedipine']
channels = ['IKr', 'INaL', 'ICaL', 'INa', 'Ito', 'IK1', 'IKs']
SD_param_names = ['Kmax', 'Ku', 'halfmax', 'n', 'Vhalf']

drug_concentrations = {
    'dofetilide': {
        'coarse': [0, 0.1, 1, 10, 30, 100, 300, 500, 1000],
        'fine': 10.0**np.linspace(-1, 3, 20),
        'lit_default': [1, 3, 10, 30],
        'Cmax': [2, 4, 6, 8]
    },
    'verapamil': {
        'coarse': [0, 0.1, 1, 30, 300, 1000, 10000, 1e5],
        'fine': 10.0**np.linspace(-1, 5, 20),
        'lit_default': [30, 100, 300, 1000],
        'Cmax': [81, 2 * 81, 3 * 81, 4 * 81]
    },
    'bepridil': {
        'coarse': [0, 0.1, 1, 30, 100, 300, 1000, 10000],
        'fine': 10.0**np.linspace(-1, 5, 20),
        'lit_default': [10, 30, 100, 300]
    },
    'terfenadine': {
        'coarse': [0, 0.1, 1, 10, 30, 100, 300, 500, 1000, 10000],
        'fine': 10.0**np.linspace(-1, 5, 20),
        'lit_default': [3, 10, 30, 100]
    },
    'cisapride': {
        'coarse': [0, 0.1, 1, 10, 30, 100, 300, 500, 1000, 3000],
        'fine': 10.0**np.linspace(-1, 3, 20),
        'lit_default': [1, 10, 100, 300]
    },
    'ranolazine': {
        'coarse': [0, 1, 30, 300, 500, 1000, 10000, 1e5, 1e6],
        'fine': 10.0**np.linspace(1, 5.5, 20),
        'lit_default': [1000, 1e4, 3e4, 1e5]
    },
    'quinidine': {
        'coarse': [0, 1, 30, 300, 500, 1000, 3000, 10000, 1e5],
        'fine': 10.0**np.linspace(-1, 5, 20),
        'lit_default': [100, 300, 1000, 10000]
    },
    'sotalol': {
        'coarse': [0, 1, 30, 100, 300, 1000, 10000, 3e4, 1e5, 3e5,
                   1e6, 1e7],
        'fine': 10.0**np.linspace(-1, 7, 20),
        'lit_default': [1e4, 3e4, 1e5, 3e5]
    },
    'chlorpromazine': {
        'coarse': [0, 1, 30, 300, 500, 1000, 3000, 10000, 1e5],
        'fine': 10.0**np.linspace(-1, 4.5, 20),
        'lit_default': [100, 300, 1000, 3000]
    },
    'ondansetron': {
        'coarse': [0, 1, 30, 300, 500, 1000, 3000, 10000, 1e5, 3e5],
        'fine': 10.0**np.linspace(-1, 5.5, 20),
        'lit_default': [300, 1000, 3000, 1e4]
    },
    'diltiazem': {
        'coarse': [0, 1, 30, 100, 300, 1000, 3000, 10000, 3e4, 1e5,
                   1e6, 1e7],
        'fine': 10.0**np.linspace(-1, 6, 20),
        'lit_default': [3000, 1e4, 3e4, 1e5]
    },
    'mexiletine': {
        'coarse': [0, 1, 30, 100, 300, 1000, 10000, 3e4, 1e5, 1e6,
                   1e7],
        'fine': 10.0**np.linspace(-1, 7, 20),
        'lit_default': [1e4, 3e4, 1e5, 3e5]
    },
}


def Milnes(t_max):
    protocol = myokit.Protocol()
    protocol.schedule(-80, 0, 800, period=t_max)
    protocol.schedule(-90, 800, 100, period=t_max)
    protocol.schedule(-80, 900, 100, period=t_max)
    protocol.schedule(-80, 11000, 14000, period=t_max)

    return protocol
