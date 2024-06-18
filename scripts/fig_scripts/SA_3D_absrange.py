import glob
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os
import pandas as pd

import modelling


model_details = modelling.model_naming
model_list = model_details.APmodel_list
model_label = ['ORd-Li', 'Grandi-Li', 'ten Tusscher-Li', 'Tomek-Li', 'ORd-Lei']

# Define directory to save figure
# data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA')
fig_dir = os.path.join(modelling.FIG_DIR, 'parameter_exploration')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)
plt.rcParams.update({'font.size': 9})


def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"


# Set up figure structure
fig = modelling.figures.FigureStructure(figsize=(10, 9), gridspec=(3, 1),
                                        height_ratios=[1.5, 1.5, 1], hspace=0.2,
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(1, 3), (1, 2), (1, 5)]
subgs = []
for i in range(3):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.15))
SA_panel = []
SA_panel.append([fig.fig.add_subplot(subgs[0][0, j], projection='3d')
                for j in range(3)])
SA_panel.append([fig.fig.add_subplot(subgs[1][0, j], projection='3d')
                 for j in range(2)])
hist_panel = [fig.fig.add_subplot(subgs[2][0, j]) for j in range(5)]

cmap = plt.get_cmap('viridis')
max_count = 0

for n, APmodel_name in enumerate(model_list):

    # Read simulated data for synthetic drugs
    data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA')
    filename = 'SA_alldrugs_' + APmodel_name + '.csv'
    df = pd.read_csv(os.path.join(data_dir, filename), header=[0, 1],
                     index_col=[0], skipinitialspace=True)

    Vhalf_list = df['param_values']['Vhalf'].values
    Kmax_list = df['param_values']['Kmax'].values
    Ku_list = df['param_values']['Ku'].values
    drug_list = df.index

    RMSError_drug = df['RMSE']['RMSE'].values
    MError_drug = df['ME']['ME'].values

    # Calculate signed RMSD
    Error_drug = np.array(RMSError_drug) * np.array(MError_drug) / \
        np.abs(np.array(MError_drug))

    # Read simulated data of virtual drugs in the parameter space
    data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration',
                            'SA_space', APmodel_name)
    file_prefix = 'SA_allparam'
    result_files = glob.glob(os.path.join(data_dir, file_prefix + '*.csv'))

    # Take 10ms difference in APD90 as significant
    error_range = 10

    # Load results and extract points where the RMSD value is within the
    # defined range
    first_iter = True
    for file in result_files:
        df = pd.read_csv(file, header=[0, 1], index_col=[0],
                         skipinitialspace=True)
        chosen_df = df.loc[df['RMSE']['RMSE'] < error_range]

        if first_iter:
            combined_df = df
            combined_chosen_df = chosen_df
            first_iter = False
        else:
            combined_chosen_df = pd.concat([combined_chosen_df, chosen_df])
            combined_df = pd.concat([combined_df, df])

    Vhalf_range = combined_df['param_values']['Vhalf'].values
    Kmax_range = combined_df['param_values']['Kmax'].values
    Ku_range = combined_df['param_values']['Ku'].values

    RMSError = combined_df['RMSE']['RMSE'].values
    MError = combined_df['ME']['ME'].values
    Error_space = RMSError * MError / np.abs(MError)

    error_key = ('error', 'signed_RMSE')
    combined_df.insert(len(combined_df.columns), error_key, Error_space)

    cmin = min(min(Error_drug), min(Error_space))
    cmax = max(max(Error_drug), max(Error_space))
    print('min and max of APD difference in the parameter space:', cmin, cmax)

    y, x, _ = hist_panel[n].hist(Error_space, 20)
    hist_panel[n].axvline(0, 0, 1, color='grey', ls='--', zorder=-1)
    hist_panel[n].fill_between([-error_range, error_range],
                               0, 3300, alpha=0.5, color='grey',
                               zorder=-2)
    max_count = max(y.max(), max_count)
    # hist_panel[n].set_title(model_label[n])

    # Plot points in the parameter space
    # cmap_norm = matplotlib.colors.Normalize(cmin, cmax)
    abs_max = max(np.abs(cmin), cmax)
    # cmap_norm = matplotlib.colors.Normalize(-abs_max, abs_max)
    cmap_norm = matplotlib.colors.TwoSlopeNorm(vmin=abs_max * -1, vcenter=0,
                                               vmax=abs_max)
    scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

    Vhalf_chosen = combined_chosen_df['param_values']['Vhalf'].values
    Kmax_chosen = combined_chosen_df['param_values']['Kmax'].values
    Ku_chosen = combined_chosen_df['param_values']['Ku'].values

    print('Total number of points within the range: ', len(Vhalf_chosen))

    # Plot points of all synthetic drugs and those with RMSD within the defined
    # range
    xmin, xmax = min(Vhalf_range), max(Vhalf_range)
    SA_panel[n].scatter(Vhalf_chosen, np.log10(Kmax_chosen),
                        np.log10(Ku_chosen), c='dimgrey', s=10, marker='o',
                        zorder=-10, alpha=0.5)
    SA_panel[n].scatter(Vhalf_list, np.log10(Kmax_list), np.log10(Ku_list),
                        c=scale_map.to_rgba(Error_drug), s=100, marker='^',
                        zorder=1)
    SA_panel[n].scatter(xmin * np.ones(len(Vhalf_list)), np.log10(Kmax_list),
                        np.log10(Ku_list), s=50, marker='o', zorder=1, c='red')

    for i in range(len(Vhalf_list)):
        SA_panel[n].plot([xmin, Vhalf_list[i]],
                         [np.log10(Kmax_list[i]), np.log10(Kmax_list[i])],
                         zs=[np.log10(Ku_list[i]), np.log10(Ku_list[i])],
                         color='red', linestyle='--', linewidth=0.7)
    # SA_panel[bottom].view_init(32, 80)
    # SA_panel[bottom].view_init(32, -25)
    SA_panel[n].view_init(32, 55)

    SA_panel[n].set_xlabel(r"$V_\mathrm{half-trap}$")
    SA_panel[n].set_ylabel(r"$K_\mathrm{max}$")
    SA_panel[n].set_zlabel(r"$K_u$")

    SA_panel[n].set_xlim(min(Vhalf_range), max(Vhalf_range))
    SA_panel[n].set_ylim(min(np.log10(Kmax_range)),
                         max(np.log10(Kmax_range)))
    SA_panel[n].set_zlim(min(np.log10(Ku_range)),
                         max(np.log10(Ku_range)))

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
    scale_map.set_array(Error_space)
    fig.fig.colorbar(scale_map, orientation='horizontal', ax=SA_panel, cax=cax)
    SA_panel[n].set_title(model_label[n])

# Save figure
for i in range(5):
    hist_panel[i].set_ylim(0, max_count + 50)

hist_panel[0].spines[['right', 'top']].set_visible(False)
for i in range(1, 5):
    hist_panel[i].sharey(hist_panel[0])
    hist_panel[i].tick_params(labelleft=False)
    hist_panel[i].spines[['right', 'top']].set_visible(False)
fig.savefig(os.path.join(fig_dir, 'SA_3D_hist_signedRMSD_absrange.pdf'),)
            # bbox_inches='tight', pad_inches=0.2)
