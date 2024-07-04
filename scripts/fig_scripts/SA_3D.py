import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os
import pandas as pd

import modelling


parser = argparse.ArgumentParser(
    description="Plot parameter space exploration outcome of an AP model")
parser.add_argument("APmodel", help="Name of AP model")
args = parser.parse_args()

APmodel_name = args.APmodel

# Define directory to save figure
fig_dir = os.path.join(modelling.FIG_DIR, 'parameter_exploration',
                       APmodel_name)
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)


def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"


# Set up figure structure
fig = modelling.figures.FigureStructure(figsize=(10, 3.5), gridspec=(1, 2),
                                        width_ratios=[2.5, 1], wspace=0.2,
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(1, 2), (1, 1)]
subgs = []
subgs.append(fig.gs[0].subgridspec(*subgridspecs[0], wspace=0.18))
subgs.append(fig.gs[1].subgridspec(*subgridspecs[1]))
SA_panel = [fig.fig.add_subplot(subgs[0][0, j], projection='3d')
            for j in range(2)]
hist_panel = fig.fig.add_subplot(subgs[1][0, 0], box_aspect=1)

# Load results of the chosen AP model
APsim = modelling.ModelSimController(APmodel_name)
APcompare = modelling.ModelComparison(APsim)
drug_axis_values, drug_error_space = APcompare.load_results(mode='drugs')
axis_values, error_space = APcompare.load_results()

# Define the similarity range of the AP model
data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA',
                        'parameter_n', APmodel_name)
overall_stats = pd.read_csv(os.path.join(data_dir, 'overall_stats.csv'),
                            index_col=[0])
error_thres = overall_stats['max'].values[0] * 1.2
error_similar = [True if np.abs(i) < error_thres else False
                 for i in error_space['error']]
print('Difference in APD where SD and CS model are similar: ', error_thres)
print('Total number of points within the range: ', np.sum(error_similar))

# Find the min and max RMSD of all results
cmin = min(min(drug_error_space['norm_error']),
           min(error_space['norm_error']))
cmax = max(max(drug_error_space['norm_error']),
           max(error_space['norm_error']))
print('min and max of APD difference in the parameter space:', cmin, cmax)

# Plot histogram of all RMSD within the parameter space
hist_panel.hist(error_space['error'], 20)
hist_panel.axvline(0, 0, 1, color='grey', ls='--', zorder=-1)
hist_panel.fill_between([-error_thres, error_thres], 0, 3300,
                        alpha=0.5, color='grey', zorder=-2)
hist_panel.spines[['right', 'top']].set_visible(False)
hist_panel.set_xlabel(r'$\Delta \mathrm{APD}_{90}$')

# Define color map
cmap = plt.get_cmap('viridis')
cmap_norm = matplotlib.colors.Normalize(cmin, cmax)
scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

# Plot RMSDs of normalised APD90
SA_panel[0].scatter(axis_values['Vhalf'], np.log10(axis_values['Kmax']),
                    np.log10(axis_values['Ku']),
                    c=scale_map.to_rgba(error_space['norm_error']),
                    marker='o', zorder=-10, alpha=0.5)
SA_panel[0].view_init(32, 55)
# SA_panel[0].view_init(32, -25)

# Plot points of 12 CiPA drugs and those with RMSD within the similarity
# range
xmin, xmax = min(axis_values['Vhalf']), max(axis_values['Vhalf'])
SA_panel[1].scatter(axis_values['Vhalf'][error_similar],
                    np.log10(axis_values['Kmax'][error_similar]),
                    np.log10(axis_values['Ku'][error_similar]),
                    c='dimgrey', s=3, marker='o', zorder=-10, alpha=0.2)
SA_panel[1].scatter(drug_axis_values['Vhalf'],
                    np.log10(drug_axis_values['Kmax']),
                    np.log10(drug_axis_values['Ku']),
                    c=scale_map.to_rgba(drug_error_space['norm_error']),
                    s=100, marker='^', zorder=-1)
SA_panel[1].scatter(xmin * np.ones(len(drug_axis_values['Vhalf'])),
                    np.log10(drug_axis_values['Kmax']),
                    np.log10(drug_axis_values['Ku']),
                    s=50, marker='o', zorder=-5, c='red')

for i in range(len(drug_axis_values['Vhalf'])):
    SA_panel[1].plot([xmin, drug_axis_values['Vhalf'][i]],
                     [np.log10(drug_axis_values['Kmax'][i]),
                      np.log10(drug_axis_values['Kmax'][i])],
                     zs=[np.log10(drug_axis_values['Ku'][i]),
                         np.log10(drug_axis_values['Ku'][i])],
                     color='red', linestyle='--', linewidth=0.7)
SA_panel[1].view_init(32, 55)
# SA_panel[1].view_init(32, -25)

# Adjust figure details
for i in range(2):
    SA_panel[i].set_xlabel(r"$V_\mathrm{half-trap}$")
    SA_panel[i].set_ylabel(r"$K_\mathrm{max}$")
    SA_panel[i].set_zlabel(r"$K_u$")

    SA_panel[i].set_xlim(min(axis_values['Vhalf']), max(axis_values['Vhalf']))
    SA_panel[i].set_ylim(min(np.log10(axis_values['Kmax'])),
                         max(np.log10(axis_values['Kmax'])))
    SA_panel[i].set_zlim(min(np.log10(axis_values['Ku'])),
                         max(np.log10(axis_values['Ku'])))

    SA_panel[i].xaxis.set_major_locator(mticker.MaxNLocator(nbins=6))
    SA_panel[i].yaxis.set_major_formatter(
        mticker.FuncFormatter(log_tick_formatter))
    SA_panel[i].yaxis.set_major_locator(mticker.MaxNLocator(nbins=8))
    SA_panel[i].zaxis.set_major_formatter(
        mticker.FuncFormatter(log_tick_formatter))
    SA_panel[i].zaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    SA_panel[i].set_rasterization_zorder(0)

cax = SA_panel[0].inset_axes([0.5, -0.15, 1, 0.03])
fig.fig.colorbar(scale_map, orientation='horizontal', ax=SA_panel, cax=cax,
                 label=r'$\Delta \widetilde{\mathrm{APD}}_{90}$')

# Add panel labels
fig.fig.text(0.075, 0.765, '(A)', fontsize=10)
fig.fig.text(0.375, 0.765, '(B)', fontsize=10)
fig.fig.text(0.655, 0.765, '(C)', fontsize=10)

# Save figure
fig.savefig(os.path.join(fig_dir, 'SA_3D_hist.pdf'))
