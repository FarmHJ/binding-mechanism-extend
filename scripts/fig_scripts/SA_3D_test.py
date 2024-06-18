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
model_label = ['Grandi-Li', 'ten Tusscher-Li', 'Tomek-Li', 'ORd-Lei']

# Define directory to save figure
# data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA')
fig_dir = os.path.join(modelling.FIG_DIR, 'parameter_exploration')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)
plt.rcParams.update({'font.size': 9})


def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"


# Set up figure structure
fig = plt.figure(figsize=(6, 6))
gs = fig.add_gridspec(2, 2, wspace=0.15, hspace=0.5)
axs = [fig.add_subplot(gs[i, j], projection='3d') for j in range(2)
       for i in range(2)]

cmap = plt.get_cmap('rainbow')
# cmap = plt.get_cmap('RdBu_r')

fig_1D = plt.figure(figsize=(7, 2))
gs_1D = fig_1D.add_gridspec(1, 3)
axs_1D = [fig_1D.add_subplot(gs_1D[0, j]) for j in range(3)]

fig_hist = plt.figure(figsize=(9, 2))
gs_hist = fig_hist.add_gridspec(1, 4, wspace=0.15)
axs_hist = [fig_hist.add_subplot(gs_hist[0, j]) for j in range(4)]
max_count = 0

APsim = modelling.ModelSimController('Grandi')
FnClass = modelling.ModelComparison(APsim)

data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA')
filename = 'SA_alldrugs_ORd-Li.csv'
df = pd.read_csv(os.path.join(data_dir, filename), header=[0, 1],
                 index_col=[0], skipinitialspace=True)

RMSError_base_drug = df['RMSE']['RMSE'].values
MError_base_drug = df['ME']['ME'].values

data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration',
                        'SA_space', 'ORd-Li')
file_prefix = 'SA_allparam'
result_files = glob.glob(os.path.join(data_dir, file_prefix + '*.csv'))

# Load results and extract points where the RMSD value is within the defined
# range
first_iter = True
for file in result_files:
    df = pd.read_csv(file, header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    if first_iter:
        base_df = df
        first_iter = False
    else:
        base_df = pd.concat([base_df, df])

RMSError_base = base_df['RMSE']['RMSE'].values
MError_base = base_df['ME']['ME'].values
print(len(RMSError_base))

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
    # Error_drug = np.array(RMSError_drug) * np.array(MError_drug) / \
    #     np.abs(np.array(MError_drug))
    # Error_drug = (np.array(RMSError_drug) - np.array(RMSError_base_drug)) \
    #     / np.array(RMSError_base_drug)

    APD_SD_list = df['APD_trapping']
    APD_CS_list = df['APD_conductance']
    Error_drug = []
    for d in range(APD_SD_list.shape[0]):
        APD_SD = APD_SD_list.iloc[d, :].values
        APD_CS = APD_CS_list.iloc[d, :].values
        APD_min = min(min(APD_SD), min(APD_CS))
        APD_max = max(max(APD_SD), max(APD_CS))
        APD_SD_norm = (APD_SD - APD_min) / (APD_max - APD_min)
        APD_CS_norm = (APD_CS - APD_min) / (APD_max - APD_min)
        FnClass.APD_trapping = APD_SD_norm
        FnClass.APD_conductance = APD_CS_norm
        FnClass.ME()
        Error_drug.append(FnClass.Error)

    # Read simulated data of virtual drugs in the parameter space
    data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration',
                            'SA_space', APmodel_name)
    file_prefix = 'SA_allparam'
    result_files = glob.glob(os.path.join(data_dir, file_prefix + '*.csv'))

    # Define the range where the RMSD between the APD90s of the two models are
    # small
    data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA',
                            'parameter_n', APmodel_name)
    overall_stats = pd.read_csv(os.path.join(data_dir, 'overall_stats.csv'),
                                index_col=[0])
    error_range = overall_stats['max'].values[0] * 1.2
    print('Difference in APD where SD and CS model are similar: ', error_range)

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

    # Error_space = RMSError * MError / np.abs(MError)
    # print(len(RMSError))
    # print(APmodel_name)
    # Error_space = (RMSError - RMSError_base) / RMSError_base

    APD_SD_list = combined_df['APD_trapping']
    APD_CS_list = combined_df['APD_conductance']
    Error_space = []
    for d in range(APD_SD_list.shape[0]):
        APD_SD = APD_SD_list.iloc[d, :].values
        APD_CS = APD_CS_list.iloc[d, :].values
        APD_min = min(min(APD_SD), min(APD_CS))
        APD_max = max(max(APD_SD), max(APD_CS))
        APD_SD_norm = (APD_SD - APD_min) / (APD_max - APD_min)
        APD_CS_norm = (APD_CS - APD_min) / (APD_max - APD_min)
        FnClass.APD_trapping = APD_SD_norm
        FnClass.APD_conductance = APD_CS_norm
        FnClass.ME()
        Error_space.append(FnClass.Error)

    error_key = ('error', 'signed_RMSE')
    combined_df.insert(len(combined_df.columns), error_key, Error_space)

    cmin = min(min(Error_drug), min(Error_space))
    cmax = max(max(Error_drug), max(Error_space))
    print('min and max of APD difference in the parameter space:', cmin, cmax)

    y, x, _ = axs_hist[n].hist(Error_space, 20)
    axs_hist[n].axvline(0, 0, 1, color='grey', zorder=-1)
    max_count = max(y.max(), max_count)
    axs_hist[n].set_title(APmodel_name)

    # Plot points in the parameter space
    cmap_norm = matplotlib.colors.Normalize(cmin, cmax)
    # abs_max = max(np.abs(cmin), cmax)
    # cmap_norm = matplotlib.colors.Normalize(-abs_max, abs_max)
    # cmap_norm = matplotlib.colors.TwoSlopeNorm(vmin=abs_max * -1, vcenter=0,
    #                                            vmax=abs_max)
    scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

    axs[n].scatter(Vhalf_range, np.log10(Kmax_range), np.log10(Ku_range),
                   c=scale_map.to_rgba(Error_space), s=5,
                   marker='o', zorder=-10, alpha=0.5)
    # axs[n].view_init(32, 55)
    axs[n].view_init(32, -25)

    Vhalf_values = set(Vhalf_range)
    Kmax_values = set(Kmax_range)
    Ku_values = set(Ku_range)
    param_values = {
        'Vhalf': Vhalf_values,
        'Kmax': Kmax_values,
        'Ku': Ku_values
    }
    param_list = ['Vhalf', 'Kmax', 'Ku']
    param_unit = ['mV', '1', r'ms$^{-1}$']
    for c, p_interest in enumerate(param_list):
        max_diff = 0
        param_xy = [i for i in param_list if i != p_interest]
        p1, p2 = param_xy[0], param_xy[1]
        diff_dict = {}
        for p1_value, p2_value in itertools.product(param_values[p1], param_values[p2]):
            temp = combined_df[(combined_df[('param_values', p1)] == p1_value)
                               & (combined_df[('param_values', p2)] == p2_value)]
            diff = max(temp[error_key]) - min(temp[error_key])
            diff_dict.update({(p1_value, p2_value): diff})

        # Get max point
        chosen_pt = max(diff_dict, key=diff_dict.get)
        max_diff = diff_dict[chosen_pt]
        chosen_id = combined_df[(combined_df[('param_values', p1)] == chosen_pt[0]) &
                                (combined_df[('param_values', p2)] == chosen_pt[1])].index
        print('chosen id: ', chosen_id)
        print('parameter of interest: ', p_interest)
        print('max diff: ', max_diff)
        print('param point: ', (p1, p2))
        print('max point: ', chosen_pt)
        print('param id: ', chosen_id.values)
        dir = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration',
                           'check', APmodel_name)
        np.savetxt(os.path.join(dir, f'chosen_id_{p_interest}.txt'),
                   chosen_id.values, fmt='%d')

        temp = combined_df[(combined_df[('param_values', p1)] == chosen_pt[0])
                           & (combined_df[('param_values', p2)] == chosen_pt[1])]
        p_interest_range = temp[('param_values', p_interest)]
        p_interest_error = temp[error_key]

        error_norm = (p_interest_error - cmin) / (cmax - cmin)
        axs_1D[c].plot(p_interest_range, error_norm, 'o', label=model_label[n])
        # axs_1D[c].set_title(p_interest)
        axs_1D[c].set_xlabel(f'{p_interest} ({param_unit[c]})')

        if p_interest == 'Vhalf':
            step = np.mean(p_interest_range[1:] - p_interest_range[:-1])
            p_interest_range = list(p_interest_range.values)
            p_interest_range += [p_interest_range[0] - step, p_interest_range[-1] + step]
        else:
            step = np.mean(np.log10(p_interest_range[1:]) - np.log10(p_interest_range[:-1]))
            p_interest_range = list(p_interest_range.values)
            p_interest_range += [np.exp(np.log(p_interest_range[0]) - step), np.exp(np.log(p_interest_range[-1]) + step)]

    # Adjust figure details
    for i in range(2):
        axs[n].set_xlabel(r"$V_\mathrm{half-trap}$")
        axs[n].set_ylabel(r"$K_\mathrm{max}$")
        axs[n].set_zlabel(r"$K_u$")

        axs[n].set_xlim(min(Vhalf_range), max(Vhalf_range))
        axs[n].set_ylim(min(np.log10(Kmax_range)), max(np.log10(Kmax_range)))
        axs[n].set_zlim(min(np.log10(Ku_range)), max(np.log10(Ku_range)))

        axs[n].xaxis.set_major_locator(mticker.MaxNLocator(nbins=6))
        axs[n].yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
        axs[n].yaxis.set_major_locator(mticker.MaxNLocator(nbins=6, integer=True))
        axs[n].zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
        axs[n].zaxis.set_major_locator(mticker.MaxNLocator(integer=True))

        axs[n].set_rasterization_zorder(0)

    cax = axs[n].inset_axes([0.08, -0.2, 0.8, 0.03])
    scale_map.set_array(Error_space)
    fig.colorbar(scale_map, orientation='horizontal', ax=axs, cax=cax)
    axs[n].set_title(model_label[n])

# Save figure
fig.savefig(os.path.join(fig_dir, 'SA_3D_normMD.png'),
            bbox_inches='tight', pad_inches=0.2)

axs_1D[1].set_xscale('log')
axs_1D[2].set_xscale('log')
axs_1D[1].legend(handlelength=1, handletextpad=0.6)
axs_1D[0].set_ylabel('Normalised max difference')
fig_1D.text(0.1, 0.925, '(A)', fontsize=10)
fig_1D.text(0.37, 0.925, '(B)', fontsize=10)
fig_1D.text(0.645, 0.925, '(C)', fontsize=10)
for i in range(3):
    axs_1D[i].set_ylim(-0.05, 1.05)
    axs_hist[i].set_ylim(0, max_count + 50)
fig_1D.savefig(os.path.join(fig_dir, 'SA_1D_normMD.pdf'),
               bbox_inches='tight')

for i in range(1, 4):
    axs_hist[i].sharey(axs_hist[0])
    axs_hist[i].tick_params(labelleft=False)
fig_hist.savefig(os.path.join(fig_dir, 'SA_hist_normMD.pdf'),
                 bbox_inches='tight')
