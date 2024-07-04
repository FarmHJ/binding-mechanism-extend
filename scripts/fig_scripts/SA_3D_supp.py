import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os

import modelling


# Take all AP models
model_details = modelling.model_naming
model_list = model_details.APmodel_list

# Define directory to save figure
fig_dir = os.path.join(modelling.FIG_DIR, 'parameter_exploration')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)
plt.rcParams.update({'font.size': 9})


def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"


# Set up figure structure
fig = modelling.figures.FigureStructure(figsize=(10, 9), gridspec=(3, 1),
                                        height_ratios=[1.5, 1.5, 1],
                                        hspace=0.75, plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(1, 3), (1, 2), (1, 5)]
subgs = []
for i in range(3):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.15))
SA_panel = [fig.fig.add_subplot(subgs[0][0, j], projection='3d')
            for j in range(3)] + \
           [fig.fig.add_subplot(subgs[1][0, j], projection='3d')
            for j in range(2)]
hist_panel = [fig.fig.add_subplot(subgs[2][0, j]) for j in range(5)]

cmap = plt.get_cmap('viridis')
max_count = 0

for n, APmodel_name in enumerate(model_list):

    # Load results of the AP model
    APsim = modelling.ModelSimController(APmodel_name)
    APcompare = modelling.ModelComparison(APsim)
    drug_axis_values, drug_error_space = APcompare.load_results(mode='drugs')
    axis_values, error_space = APcompare.load_results()

    # Define the similarity range as 10ms
    error_thres = 10
    error_similar = [True if np.abs(i) < error_thres else False
                     for i in error_space['error']]

    # Find the min and max RMSD of all results
    cmin = min(min(drug_error_space['norm_error']),
               min(error_space['norm_error']))
    cmax = max(max(drug_error_space['norm_error']),
               max(error_space['norm_error']))
    print('min and max of APD difference in the parameter space:', cmin, cmax)
    print('Total number of points within the range: ', np.sum(error_similar))

    # Define color map
    abs_max = max(np.abs(cmin), cmax)
    cmap_norm = matplotlib.colors.TwoSlopeNorm(vmin=abs_max * -1, vcenter=0,
                                               vmax=abs_max)
    scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

    # Plot histogram of all RMSD within the parameter space
    y, x, _ = hist_panel[n].hist(error_space['error'], 20)
    hist_panel[n].axvline(0, 0, 1, color='grey', ls='--', zorder=-1)
    hist_panel[n].fill_between([-error_thres, error_thres],
                               0, 3300, alpha=0.5, color='grey',
                               zorder=-2)
    max_count = max(y.max(), max_count)
    hist_panel[n].set_title(
        model_details.AP_file_names[APmodel_name]['label'])

    # Plot points of 12 CiPA drugs and those with RMSD within the similarity
    # range
    xmin, xmax = min(axis_values['Vhalf']), max(axis_values['Vhalf'])
    SA_panel[n].scatter(axis_values['Vhalf'][error_similar],
                        np.log10(axis_values['Kmax'][error_similar]),
                        np.log10(axis_values['Ku'][error_similar]),
                        c='dimgrey', s=10, marker='o', zorder=-10, alpha=0.5)
    SA_panel[n].scatter(drug_axis_values['Vhalf'],
                        np.log10(drug_axis_values['Kmax']),
                        np.log10(drug_axis_values['Ku']),
                        c=scale_map.to_rgba(drug_error_space['norm_error']),
                        s=100, marker='^', zorder=1)
    SA_panel[n].scatter(xmin * np.ones(len(drug_axis_values['Vhalf'])),
                        np.log10(drug_axis_values['Kmax']),
                        np.log10(drug_axis_values['Ku']),
                        s=50, marker='o', zorder=1, c='red')

    for i in range(len(drug_axis_values['Vhalf'])):
        SA_panel[n].plot([xmin, drug_axis_values['Vhalf'][i]],
                         [np.log10(drug_axis_values['Kmax'][i]),
                          np.log10(drug_axis_values['Kmax'][i])],
                         zs=[np.log10(drug_axis_values['Ku'][i]),
                             np.log10(drug_axis_values['Ku'][i])],
                         color='red', linestyle='--', linewidth=0.7)
    # SA_panel[bottom].view_init(32, 80)
    # SA_panel[bottom].view_init(32, -25)
    SA_panel[n].view_init(32, 55)

    # Adjust figure
    SA_panel[n].set_xlabel(r"$V_\mathrm{half-trap}$")
    SA_panel[n].set_ylabel(r"$K_\mathrm{max}$")
    SA_panel[n].set_zlabel(r"$K_u$")

    SA_panel[n].set_xlim(min(axis_values['Vhalf']),
                         max(axis_values['Vhalf']))
    SA_panel[n].set_ylim(min(np.log10(axis_values['Kmax'])),
                         max(np.log10(axis_values['Kmax'])))
    SA_panel[n].set_zlim(min(np.log10(axis_values['Ku'])),
                         max(np.log10(axis_values['Ku'])))

    SA_panel[n].xaxis.set_major_locator(mticker.MaxNLocator(nbins=6))
    SA_panel[n].yaxis.set_major_formatter(
        mticker.FuncFormatter(log_tick_formatter))
    SA_panel[n].yaxis.set_major_locator(
        mticker.MaxNLocator(nbins=6, integer=True))
    SA_panel[n].zaxis.set_major_formatter(
        mticker.FuncFormatter(log_tick_formatter))
    SA_panel[n].zaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    SA_panel[n].set_rasterization_zorder(0)

    cax = SA_panel[n].inset_axes([0.08, -0.3, 0.8, 0.03])
    fig.fig.colorbar(scale_map, orientation='horizontal', ax=SA_panel, cax=cax,
                     label=r'$\Delta \widetilde{\mathrm{APD}}_{90}$')
    SA_panel[n].set_title(model_details.AP_file_names[APmodel_name]['label'])

    ##############################
    # Plot all RMSD in 2D layers
    # Set up figure structure
    fig_layer = plt.figure(figsize=(2 * 5, 2 * 4))
    gs = fig_layer.add_gridspec(4, 5, wspace=0.13, hspace=0.2)
    axs_layer = [[fig_layer.add_subplot(gs[i, j]) for j in range(5)]
                 for i in range(4)]

    # Get unique Ku values
    Ku_values = np.unique(axis_values['Ku'])
    for i, Ku in enumerate(Ku_values):
        # Get index where all Ku are as defined
        Ku_ind = axis_values['Ku'] == Ku

        # Get list of corresponding Vhalf, Kmax and RMSDs
        x_grid = axis_values['Vhalf'][Ku_ind]
        y_grid = np.log10(axis_values['Kmax'][Ku_ind])
        diff_list = np.array(error_space['norm_error'])[Ku_ind]

        # Plot the RMSD for all Vhalf and Kmax for a given Ku
        r, c = int(i / 5), i % 5
        im = axs_layer[r][c].scatter(x_grid, y_grid,
                                     c=scale_map.to_rgba(diff_list))
        axs_layer[r][c].yaxis.set_major_formatter(
            mticker.FuncFormatter(log_tick_formatter))

        # Adjust figures
        axs_layer[r][c].set_title(r"$K_u = $"'{:.2e}'.format(Ku),
                                  fontdict={'fontsize': 9})
        if c == 0:
            axs_layer[r][c].set_ylabel(r"$K_\mathrm{max}$")
        else:
            axs_layer[r][c].tick_params(labelleft=False)
        if r == 3:
            axs_layer[r][c].set_xlabel(r"$V_\mathrm{half-trap}$")
        else:
            axs_layer[r][c].tick_params(labelbottom=False)
    cax = axs_layer[0][4].inset_axes([1.05, -3, 0.1, 3.5])
    fig_layer.colorbar(scale_map, ax=axs_layer, cax=cax)

    # Save figure
    fig_layer.savefig(os.path.join(fig_dir, f'SA_3D_layer_{APmodel_name}.png'),
                      bbox_inches='tight', pad_inches=0.2)

# Adjust figure
hist_panel[0].set_ylim(0, max_count + 50)
hist_panel[0].set_xlabel(r'$\Delta \mathrm{APD}_{90}$')
hist_panel[0].spines[['right', 'top']].set_visible(False)
for i in range(1, 5):
    hist_panel[i].set_ylim(0, max_count + 50)
    hist_panel[i].sharey(hist_panel[0])
    hist_panel[i].tick_params(labelleft=False)
    hist_panel[i].spines[['right', 'top']].set_visible(False)
    hist_panel[i].set_xlabel(r'$\Delta \mathrm{APD}_{90}$')

# Save figure
fig.savefig(os.path.join(fig_dir, 'SA_3D_hist_absrange.pdf'))
