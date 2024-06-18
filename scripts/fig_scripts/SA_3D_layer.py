import glob
import itertools
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os
import pandas as pd
# import statistics

import modelling


model_details = modelling.model_naming
model_list = model_details.APmodel_list[1:]
# model_list = [model_details.APmodel_list[-1]]
model_label = ['Grandi-Li', 'ten Tusscher-Li', 'Tomek-Li', 'ORd-Lei']
# model_label = ['ORd-Lei']

# Define directory to save figure
# data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA')
fig_dir = os.path.join(modelling.FIG_DIR, 'parameter_exploration')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)
plt.rcParams.update({'font.size': 9})


def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"


cmap = plt.get_cmap('viridis')
# cmap = plt.get_cmap('RdBu_r')

APsim = modelling.ModelSimController('Grandi')
FnClass = modelling.ModelComparison(APsim)

# data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA')
# filename = 'SA_alldrugs_ORd-Li.csv'
# df = pd.read_csv(os.path.join(data_dir, filename), header=[0, 1],
#                  index_col=[0], skipinitialspace=True)

# RMSError_base_drug = df['RMSE']['RMSE'].values
# MError_base_drug = df['ME']['ME'].values

# data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration',
#                         'SA_space', 'ORd-Li')
# file_prefix = 'SA_allparam'
# result_files = glob.glob(os.path.join(data_dir, file_prefix + '*.csv'))

# # Load results and extract points where the RMSD value is within the defined
# # range
# first_iter = True
# for file in result_files:
#     df = pd.read_csv(file, header=[0, 1], index_col=[0],
#                      skipinitialspace=True)

#     if first_iter:
#         base_df = df
#         first_iter = False
#     else:
#         base_df = pd.concat([base_df, df])

# RMSError_base = base_df['RMSE']['RMSE'].values
# MError_base = base_df['ME']['ME'].values

for n, APmodel_name in enumerate(model_list):
    print(APmodel_name)
    # Set up figure structure
    fig = plt.figure(figsize=(2 * 5, 2 * 4))
    gs = fig.add_gridspec(4, 5, wspace=0.13, hspace=0.2)
    axs = [[fig.add_subplot(gs[i, j]) for j in range(5)] for i in range(4)]

    # Read simulated data of virtual drugs in the parameter space
    data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration',
                            'SA_space', APmodel_name)
    file_prefix = 'SA_allparam'
    result_files = glob.glob(os.path.join(data_dir, file_prefix + '*.csv'))

    # Load results and extract points where the RMSD value is within the defined
    # range
    first_iter = True
    for file in result_files:
        df = pd.read_csv(file, header=[0, 1], index_col=[0],
                         skipinitialspace=True)

        if first_iter:
            combined_df = df
            first_iter = False
        else:
            combined_df = pd.concat([combined_df, df])

    Vhalf_range = combined_df['param_values']['Vhalf'].values
    Kmax_range = combined_df['param_values']['Kmax'].values
    Ku_range = combined_df['param_values']['Ku'].values

    RMSError = combined_df['RMSE']['RMSE'].values
    MError = combined_df['ME']['ME'].values

    Error_space = RMSError * MError / np.abs(MError)
    error_key = ('error', 'signed_RMSE')
    combined_df.insert(len(combined_df.columns), error_key, Error_space)

    APD_SD_list = combined_df['APD_trapping']
    APD_CS_list = combined_df['APD_conductance']
    Error_norm_space = []
    for d in range(APD_SD_list.shape[0]):
        APD_SD = APD_SD_list.iloc[d, :].values
        APD_CS = APD_CS_list.iloc[d, :].values
        APD_min = min(min(APD_SD), min(APD_CS))
        APD_max = max(max(APD_SD), max(APD_CS))
        APD_SD_norm = (APD_SD - APD_min) / (APD_max - APD_min)
        APD_CS_norm = (APD_CS - APD_min) / (APD_max - APD_min)
        FnClass.APD_trapping = APD_SD_norm
        FnClass.APD_conductance = APD_CS_norm
        FnClass.RMSE()
        RMSD = FnClass.Error
        FnClass.ME()
        MD = FnClass.Error
        Error_norm_space.append(RMSD * MD / np.abs(MD))
    combined_df.insert(len(combined_df.columns), ('error', 'signed_norm_RMSE'),
                       Error_norm_space)

    cmin = min(Error_norm_space)
    cmax = max(Error_norm_space)
    print('min and max of APD difference in the parameter space:', cmin, cmax)
    print('shape of parameter space', combined_df.shape)

    # Plot points in the parameter space
    # cmap_norm = matplotlib.colors.Normalize(cmin, cmax)
    abs_max = max(np.abs(cmin), cmax)
    cmap_norm = matplotlib.colors.TwoSlopeNorm(vmin=abs_max * -1, vcenter=0,
                                               vmax=abs_max)
    scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

    third_dim = sorted(combined_df['param_values']['Ku'].unique())
    for i, Ku in enumerate(third_dim):

        r, c = int(i / 5), i % 5
        layer_df = combined_df[combined_df['param_values']['Ku'] == Ku]
        layer_df = layer_df.sort_values(by=[('param_values', 'Kmax'),
                                            ('param_values', 'Vhalf')])

        x_grid = layer_df['param_values']['Vhalf'].values
        y_grid = np.log10(layer_df['param_values']['Kmax'].values)
        diff_list = layer_df['error']['signed_norm_RMSE'].values
        # diff_grid = np.asarray(diff_list).reshape((20, 20))[::-1, :]

        xmin, xmax = np.min(x_grid), np.max(x_grid)
        ymin, ymax = np.min(y_grid), np.max(y_grid)
        extent = [xmin, xmax, ymin, ymax]

        # im = axs[r][c].imshow(diff_grid, extent=extent, aspect='auto',
        #                       norm=cmap_norm)
        im = axs[r][c].scatter(x_grid, y_grid, c=scale_map.to_rgba(diff_list))
        axs[r][c].yaxis.set_major_formatter(
            mticker.FuncFormatter(log_tick_formatter))

        axs[r][c].set_title(r"$K_u = $"'{:.2e}'.format(Ku),
                            fontdict={'fontsize': 9})
        if c == 0:
            axs[r][c].set_ylabel(r"$K_\mathrm{max}$")
        else:
            axs[r][c].tick_params(labelleft=False)
        if r == 3:
            axs[r][c].set_xlabel(r"$V_\mathrm{half-trap}$")
        else:
            axs[r][c].tick_params(labelbottom=False)

    cax = axs[0][4].inset_axes([1.05, -3, 0.1, 3.5])
    fig.colorbar(scale_map, ax=axs, cax=cax)

    # Save figure
    fig.savefig(os.path.join(fig_dir,
                             f'SA_3D_layer_{APmodel_name}_normsignedRMSD.png'),
                bbox_inches='tight', pad_inches=0.2)
    plt.close()
