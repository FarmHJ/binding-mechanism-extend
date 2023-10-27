import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

import modelling


class DataLibrary(object):
    """
    A data library class that reads the experimental data
    """
    def __init__(self):
        super(DataLibrary, self).__init__()

    def set_drug(self, drug):
        self.drug = drug

        self.data = pd.read_csv(os.path.join(modelling.DATA_DIR, self.drug + '.csv'))
        self.n_exp = max(self.data['exp'])
        self.n_sweeps = 10

        self.concs = np.sort(pd.unique(self.data['conc'].values))

    def get_conc_reps(self):
        conc_reps = {}
        for c in self.concs:
            reps = max(self.data.loc[self.data['conc'] == c, 'exp'])
            conc_reps.update({c: reps})

        return conc_reps

    def plot_signal(self):

        conc_reps = self.get_conc_reps()
        print(conc_reps)

        cmap = matplotlib.cm.get_cmap('viridis')
        norm = matplotlib.colors.Normalize(0, self.n_sweeps)

        fig = plt.figure(figsize=(2 * self.n_exp, len(self.concs)))
        gs = fig.add_gridspec(len(self.concs), self.n_exp,
                              hspace=0.2, wspace=0.25)
        axs = [[fig.add_subplot(gs[i, j]) for j in
                range(conc_reps[self.concs[i]])] for i in
               range(len(self.concs))]

        for c, conc in enumerate(self.concs):
            for n in range(conc_reps[conc]):
                subset_data = self.data.loc[(self.data['conc'] == conc)
                                            & (self.data['exp'] == n + 1)]
                for s in range(self.n_sweeps):
                    subset_data_sweep = subset_data.loc[
                        subset_data['sweep'] == s]
                    axs[c][n - 1].plot(subset_data_sweep['time'],
                                       subset_data_sweep['frac'],
                                       color=cmap(norm(s)))

        fig.savefig(os.path.join(modelling.FIG_DIR, 'checking_figures',
                                 'Li-data.pdf'), bbox_inches='tight')

