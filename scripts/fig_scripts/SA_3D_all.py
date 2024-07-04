import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os
import pandas as pd

import modelling


# Take AP modified models
model_details = modelling.model_naming
model_list = model_details.APmodel_list[1:-1]

# Define directory to save figure
fig_dir = os.path.join(modelling.FIG_DIR, 'parameter_exploration')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)


def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"


# Set up figure structure
fig = modelling.figures.FigureStructure(figsize=(10, 9), gridspec=(2, 1),
                                        height_ratios=[3, 1], hspace=0.25,
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()
plt.rcParams.update({'font.size': 9})
subgridspecs = [(2, 3), (1, 4)]
subgs = []
subgs.append(fig.gs[0].subgridspec(*subgridspecs[0], wspace=0.1,
                                   hspace=0.41, height_ratios=[1, 1]))
subgs.append(fig.gs[1].subgridspec(*subgridspecs[1], wspace=0.15))
SA_panel = [fig.fig.add_subplot(subgs[0][i, j], projection='3d')
            for j in range(3) for i in range(2)]
hist_panel = [fig.fig.add_subplot(subgs[1][0, j]) for j in range(4)]
cmap = plt.get_cmap('viridis')

################################################################
# Adapt results for ORd-Li model from previously published paper
# Load results
APsim = modelling.ModelSimController('ORd-Li')
APcompare = modelling.ModelComparison(APsim)
axis_values, error_space = APcompare.load_results()

# Load the max RMSD in varying parameter n for ORd-Li model
# Define the similarity range based on the max RMSD
data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA',
                        'parameter_n', 'ORd-Li')
overall_stats = pd.read_csv(os.path.join(data_dir, 'overall_stats.csv'),
                            index_col=[0])
error_thres = overall_stats['max'].values[0] * 1.2

# Plot histogram of all RMSD within the parameter space
max_count = 0
y, x, _ = hist_panel[0].hist(error_space['error'], 20)
hist_panel[0].axvline(0, 0, 1, color='grey', ls='--', zorder=-1)
hist_panel[0].fill_between([-error_thres, error_thres],
                           0, 3000, alpha=0.5, color='grey',
                           zorder=-2)
max_count = max(y.max(), max_count)
hist_panel[0].set_title('ORd-Li')
hist_panel[0].set_xlabel(r'$\Delta \mathrm{APD}_{90}$')

###################################
# Plot results for all AP-Li models
for n, APmodel_name in enumerate(model_list):

    APsim = modelling.ModelSimController(APmodel_name)
    APcompare = modelling.ModelComparison(APsim)

    # Load results of the parameter space and the 12 CiPA drugs
    drug_axis_values, drug_error_space = APcompare.load_results(mode='drugs')
    axis_values, error_space = APcompare.load_results()

    # Find the min and max RMSD of all results
    cmin = min(min(drug_error_space['norm_error']),
               min(error_space['norm_error']))
    cmax = max(max(drug_error_space['norm_error']),
               max(error_space['norm_error']))
    print('min and max of APD difference in the parameter space:', cmin, cmax)

    # Define the color map
    abs_max = max(np.abs(cmin), cmax)
    cmap_norm = matplotlib.colors.TwoSlopeNorm(vmin=abs_max * -1, vcenter=0,
                                               vmax=abs_max)
    scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

    # Plot RMSD of normalised APD90
    top = n * 2
    btm = n * 2 + 1
    SA_panel[top].scatter(axis_values['Vhalf'], np.log10(axis_values['Kmax']),
                          np.log10(axis_values['Ku']),
                          c=scale_map.to_rgba(error_space['norm_error']),
                          marker='o', zorder=-10, alpha=0.5)
    SA_panel[top].view_init(32, 55)
    SA_panel[top].set_rasterization_zorder(0)
    # SA_panel[top].view_init(32, -25)

    # Define the similarity range
    # Take points that lie within the similarity range
    data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA',
                            'parameter_n', APmodel_name)
    overall_stats = pd.read_csv(os.path.join(data_dir, 'overall_stats.csv'),
                                index_col=[0])
    error_thres = overall_stats['max'].values[0] * 1.2
    error_similar = [True if np.abs(i) < error_thres else False
                     for i in error_space['error']]
    print('Difference in APD where SD and CS model are similar: ', error_thres)
    print('Total number of points within the range: ', np.sum(error_similar))

    # Plot histogram of all RMSD within the parameter space
    y, x, _ = hist_panel[n + 1].hist(error_space['error'], 20)
    hist_panel[n + 1].axvline(0, 0, 1, color='grey', ls='--', zorder=-1)
    hist_panel[n + 1].fill_between([-error_thres, error_thres],
                                   0, 3000, alpha=0.5, color='grey',
                                   zorder=-2)
    max_count = max(y.max(), max_count)
    hist_panel[n + 1].set_title(
        model_details.AP_file_names[APmodel_name]['label'])

    # Plot points of 12 CiPA drugs and those with RMSD within the similarity
    # range
    xmin, xmax = min(axis_values['Vhalf']), max(axis_values['Vhalf'])
    SA_panel[btm].scatter(drug_axis_values['Vhalf'],
                          np.log10(drug_axis_values['Kmax']),
                          np.log10(drug_axis_values['Ku']),
                          c=scale_map.to_rgba(drug_error_space['norm_error']),
                          s=100, marker='^', zorder=-5)
    SA_panel[btm].scatter(axis_values['Vhalf'][error_similar],
                          np.log10(axis_values['Kmax'][error_similar]),
                          np.log10(axis_values['Ku'][error_similar]),
                          c='dimgrey', s=3, marker='o', zorder=-10, alpha=0.2)
    SA_panel[btm].scatter(xmin * np.ones(len(drug_axis_values['Vhalf'])),
                          np.log10(drug_axis_values['Kmax']),
                          np.log10(drug_axis_values['Ku']),
                          s=50, marker='o', zorder=-5, c='red')

    for i in range(len(drug_axis_values['Vhalf'])):
        SA_panel[btm].plot([xmin, drug_axis_values['Vhalf'][i]],
                           [np.log10(drug_axis_values['Kmax'][i]),
                            np.log10(drug_axis_values['Kmax'][i])],
                           zs=[np.log10(drug_axis_values['Ku'][i]),
                               np.log10(drug_axis_values['Ku'][i])],
                           color='red', linestyle='--', linewidth=0.7)
    # SA_panel[bottom].view_init(32, 80)
    # SA_panel[bottom].view_init(32, -25)
    SA_panel[btm].view_init(32, 55)

    # Adjust figure details
    for i in range(2):
        SA_panel[2 * n + i].set_xlabel(r"$V_\mathrm{half-trap}$")
        SA_panel[2 * n + i].set_ylabel(r"$K_\mathrm{max}$")
        SA_panel[2 * n + i].set_zlabel(r"$K_u$")

        SA_panel[2 * n + i].set_xlim(min(axis_values['Vhalf']),
                                     max(axis_values['Vhalf']))
        SA_panel[2 * n + i].set_ylim(min(np.log10(axis_values['Kmax'])),
                                     max(np.log10(axis_values['Kmax'])))
        SA_panel[2 * n + i].set_zlim(min(np.log10(axis_values['Ku'])),
                                     max(np.log10(axis_values['Ku'])))

        SA_panel[2 * n + i].xaxis.set_major_locator(
            mticker.MaxNLocator(nbins=6))
        SA_panel[2 * n + i].yaxis.set_major_formatter(
            mticker.FuncFormatter(log_tick_formatter))
        SA_panel[2 * n + i].yaxis.set_major_locator(
            mticker.MaxNLocator(nbins=6, integer=True))
        SA_panel[2 * n + i].zaxis.set_major_formatter(
            mticker.FuncFormatter(log_tick_formatter))
        SA_panel[2 * n + i].zaxis.set_major_locator(
            mticker.MaxNLocator(integer=True))

        SA_panel[2 * n + i].set_rasterization_zorder(0)

    cax = SA_panel[2 * n].inset_axes([0.08, -0.2, 0.8, 0.03])
    fig.fig.colorbar(scale_map, orientation='horizontal', ax=SA_panel, cax=cax,
                     label=r'$\Delta \widetilde{\mathrm{APD}}_{90}$')
    SA_panel[2 * n].set_title(
        model_details.AP_file_names[APmodel_name]['label'])

fig.fig.text(0.125, 0.88, '(A)', fontsize=10)
fig.fig.text(0.375, 0.88, '(B)', fontsize=10)
fig.fig.text(0.675, 0.88, '(C)', fontsize=10)
fig.fig.text(0.1, 0.29, '(D)', fontsize=10)

for i in range(3):
    hist_panel[i].set_ylim(0, max_count + 50)

hist_panel[0].spines[['right', 'top']].set_visible(False)
for i in range(1, 4):
    hist_panel[i].sharey(hist_panel[0])
    hist_panel[i].tick_params(labelleft=False)
    hist_panel[i].spines[['right', 'top']].set_visible(False)
    hist_panel[i].set_xlabel(r'$\Delta \mathrm{APD}_{90}$')

# Save figure
fig.savefig(os.path.join(fig_dir, 'SA_3D_hist.pdf'))
