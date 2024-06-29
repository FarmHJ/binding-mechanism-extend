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
model_list = model_details.APmodel_list[1:-1]
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
fig = modelling.figures.FigureStructure(figsize=(10, 9), gridspec=(2, 1),
                                        height_ratios=[3, 1], hspace=0.25,
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(2, 3), (1, 4)]
subgs = []
subgs.append(fig.gs[0].subgridspec(*subgridspecs[0], wspace=0.1,
                                   hspace=0.41, height_ratios=[1, 1]))
subgs.append(fig.gs[1].subgridspec(*subgridspecs[1], wspace=0.15))
SA_panel = [fig.fig.add_subplot(subgs[0][i, j], projection='3d') for j in range(3)
            for i in range(2)]
hist_panel = [fig.fig.add_subplot(subgs[1][0, j]) for j in range(4)]

cmap = plt.get_cmap('viridis')
# cmap = plt.get_cmap('RdBu_r')

fig_1D = plt.figure(figsize=(7, 2))
gs_1D = fig_1D.add_gridspec(1, 3)
axs_1D = [fig_1D.add_subplot(gs_1D[0, j]) for j in range(3)]

# fig_hist = plt.figure(figsize=(9, 2))
# gs_hist = fig_hist.add_gridspec(1, 4, wspace=0.15)
# axs_hist = [fig_hist.add_subplot(gs_hist[0, j]) for j in range(4)]
max_count = 0

APsim = modelling.ModelSimController('Grandi')
FnClass = modelling.ModelComparison(APsim)

# Read simulated data of virtual drugs in the parameter space
data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration',
                        'SA_space', 'ORd-Li')
file_prefix = 'SA_allparam'
result_files = glob.glob(os.path.join(data_dir, file_prefix + '*.csv'))

# Load results
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

data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA',
                        'parameter_n', 'ORd-Li')
overall_stats = pd.read_csv(os.path.join(data_dir, 'overall_stats.csv'),
                            index_col=[0])
error_range = overall_stats['max'].values[0] * 1.2

# # Use 10% difference in APD90 as similarity range threshold
# APsim = modelling.ModelSimController('ORd-Li')
# log = APsim.simulate()
# apd90 = APsim.APD90(log)
# error_range = apd90 * 0.1

# # Use 10ms different in APD90 as similarity range threshold
# error_range = 10

y, x, _ = hist_panel[0].hist(Error_space, 20)
hist_panel[0].axvline(0, 0, 1, color='grey', ls='--', zorder=-1)
hist_panel[0].fill_between([-error_range, error_range],
                           0, 3000, alpha=0.5, color='grey',
                           zorder=-2)
max_count = max(y.max(), max_count)
hist_panel[0].set_title('ORd-Li')
hist_panel[0].set_xlabel(r'$\Delta \mathrm{APD}_{90}$')

for n, APmodel_name in enumerate(model_list):

    # Read simulated data for synthetic drugs
    data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA')
    filename = f'SA_alldrugs_{APmodel_name}.csv'
    df = pd.read_csv(os.path.join(data_dir, filename), header=[0, 1],
                     index_col=[0], skipinitialspace=True)

    Vhalf_list = df['param_values']['Vhalf'].values
    Kmax_list = df['param_values']['Kmax'].values
    Ku_list = df['param_values']['Ku'].values
    drug_list = df.index

    RMSError_drug = df['RMSE']['RMSE'].values
    MError_drug = df['ME']['ME'].values

    APD_SD_list = df['APD_trapping']
    APD_CS_list = df['APD_conductance']
    Error_norm_drug = []
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
        Error_norm_drug.append(RMSD * MD / np.abs(MD))

    # Calculate signed RMSD
    Error_drug = np.array(RMSError_drug) * np.array(MError_drug) / \
        np.abs(np.array(MError_drug))
    df.insert(len(df.columns), ('error', 'signed_RMSE'), Error_drug)

    # Read simulated data of virtual drugs in the parameter space
    data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration',
                            'SA_space', APmodel_name)
    file_prefix = 'SA_allparam'
    result_files = glob.glob(os.path.join(data_dir, file_prefix + '*.csv'))

    ## Define the range where the RMSD between the APD90s of the two models are
    ## small
    data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA',
                            'parameter_n', APmodel_name)
    overall_stats = pd.read_csv(os.path.join(data_dir, 'overall_stats.csv'),
                                index_col=[0])
    error_range = overall_stats['max'].values[0] * 1.2
    print('Difference in APD where SD and CS model are similar: ', error_range)
    ## Take 10% increase in APD90 as significant
    # APsim = modelling.ModelSimController(APmodel_name)
    # APsim.set_ikr_rescale_method('AP_duration')
    # log = APsim.simulate()
    # apd90 = APsim.APD90(log)
    # error_range = apd90 * 0.1
    # Take 10ms difference in APD90 as significant
    # error_range = 10

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

    error_key = ('error', 'signed_RMSE')
    combined_df.insert(len(combined_df.columns), error_key, Error_space)

    cmin = min(min(Error_norm_drug), min(Error_norm_space))
    cmax = max(max(Error_norm_drug), max(Error_norm_space))
    print('min and max of APD difference in the parameter space:', cmin, cmax)

    # hist_panel = axs[0][1]
    y, x, _ = hist_panel[n + 1].hist(Error_space, 20)
    hist_panel[n + 1].axvline(0, 0, 1, color='grey', ls='--', zorder=-1)
    hist_panel[n + 1].fill_between([-error_range, error_range],
                                   0, 3000, alpha=0.5, color='grey',
                                   zorder=-2)
    max_count = max(y.max(), max_count)
    hist_panel[n + 1].set_title(model_label[n])

    # Plot points in the parameter space
    # cmap_norm = matplotlib.colors.Normalize(cmin, cmax)
    abs_max = max(np.abs(cmin), cmax)
    # cmap_norm = matplotlib.colors.Normalize(-abs_max, abs_max)
    cmap_norm = matplotlib.colors.TwoSlopeNorm(vmin=abs_max * -1, vcenter=0,
                                               vmax=abs_max)
    scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

    # SA_panel = axs[0][0]
    top = n * 2
    bottom = n * 2 + 1
    SA_panel[top].scatter(Vhalf_range, np.log10(Kmax_range), np.log10(Ku_range),
                          c=scale_map.to_rgba(Error_norm_space),  # s=5,
                          marker='o', zorder=-10, alpha=0.5)
    # axs[top].scatter([V_Ku] * len(Ku_error), [np.log10(K_Ku)] * len(Ku_error),
    #                  np.log10(Ku_line),
    #                  c=['k'] * len(Ku_error), s=15,
    #                  marker='o', zorder=-10, alpha=0.8)
    # axs[top].scatter(V_line, [np.log10(Km_V)] * len(V_error),
    #                  [np.log10(Ku_V)] * len(V_error),
    #                  c=['k'] * len(V_error), s=15,
    #                  marker='o', zorder=-10, alpha=0.8)
    SA_panel[top].view_init(32, 55)
    SA_panel[top].set_rasterization_zorder(0)
    # SA_panel[top].view_init(32, -25)

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
            # if diff > max_diff:
            #     max_diff = diff
            #     chosen_pt1 = p1_value
            #     chosen_pt2 = p2_value
            #     p_interest_range = temp[('param_values', p_interest)]
            #     p_interest_error = temp[error_key]

        # Get max point
        chosen_pt = max(diff_dict, key=diff_dict.get)
        max_diff = diff_dict[chosen_pt]
        # if APmodel_name == 'Grandi':
        #     chosen_id = combined_df[(combined_df[('param_values', p1)] == chosen_pt[0]) &
        #                             (combined_df[('param_values', p2)] == chosen_pt[1])].loc[:, ('param_id', 'param_id')]
        # else:
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

        # # Get median point
        # median_diff = statistics.median(diff_dict.values())
        # closest = 10
        # for i in diff_dict.keys():
        #     diff = np.abs(diff_dict[i] - median_diff)
        #     if diff < closest:
        #         closest = diff
        #         chosen_pt1, chosen_pt2 = i
        # print('median diff: ', median_diff)
        # print('median point: ', (chosen_pt1, chosen_pt2))

        temp = combined_df[(combined_df[('param_values', p1)] == chosen_pt[0]) & (combined_df[('param_values', p2)] == chosen_pt[1])]
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
        # line_plot = {p_interest: p_interest_range,
        #              p1: [chosen_pt[0]] * 22,
        #              p2: [chosen_pt[1]] * 22}
        # SA_panel[bottom].scatter(line_plot['Vhalf'], np.log10(line_plot['Kmax']),
        #                          np.log10(line_plot['Ku']), c=['k'] * 22, s=15,
        #                          marker='o', zorder=-10, alpha=0.8)

    Vhalf_chosen = combined_chosen_df['param_values']['Vhalf'].values
    Kmax_chosen = combined_chosen_df['param_values']['Kmax'].values
    Ku_chosen = combined_chosen_df['param_values']['Ku'].values

    print('Total number of points within the range: ', len(Vhalf_chosen))

    # Plot points of all synthetic drugs and those with RMSD within the defined
    # range
    xmin, xmax = min(Vhalf_range), max(Vhalf_range)
    SA_panel[bottom].scatter(Vhalf_list, np.log10(Kmax_list),
                             np.log10(Ku_list),
                             c=scale_map.to_rgba(Error_norm_drug),
                             s=100, marker='^', zorder=-5)
    SA_panel[bottom].scatter(Vhalf_chosen, np.log10(Kmax_chosen),
                             np.log10(Ku_chosen), c='dimgrey',
                             s=3, marker='o', zorder=-10,
                             alpha=0.2)
    SA_panel[bottom].scatter(xmin * np.ones(len(Vhalf_list)),
                             np.log10(Kmax_list), np.log10(Ku_list),
                             s=50, marker='o', zorder=-5, c='red')

    for i in range(len(Vhalf_list)):
        SA_panel[bottom].plot([xmin, Vhalf_list[i]],
                              [np.log10(Kmax_list[i]), np.log10(Kmax_list[i])],
                              zs=[np.log10(Ku_list[i]), np.log10(Ku_list[i])],
                              color='red', linestyle='--', linewidth=0.7)
    # SA_panel[bottom].view_init(32, 80)
    # SA_panel[bottom].view_init(32, -25)
    SA_panel[bottom].view_init(32, 55)

    # Adjust figure details
    for i in range(2):
        SA_panel[2 * n + i].set_xlabel(r"$V_\mathrm{half-trap}$")
        SA_panel[2 * n + i].set_ylabel(r"$K_\mathrm{max}$")
        SA_panel[2 * n + i].set_zlabel(r"$K_u$")

        SA_panel[2 * n + i].set_xlim(min(Vhalf_range), max(Vhalf_range))
        SA_panel[2 * n + i].set_ylim(min(np.log10(Kmax_range)),
                                     max(np.log10(Kmax_range)))
        SA_panel[2 * n + i].set_zlim(min(np.log10(Ku_range)),
                                     max(np.log10(Ku_range)))

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
    scale_map.set_array(Error_space)
    fig.fig.colorbar(scale_map, orientation='horizontal', ax=SA_panel, cax=cax,
                     label=r'$\Delta \widetilde{\mathrm{APD}}_{90}$')
    SA_panel[2 * n].set_title(model_label[n])

# Save figure
fig.fig.text(0.125, 0.88, '(A)', fontsize=10)
fig.fig.text(0.375, 0.88, '(B)', fontsize=10)
fig.fig.text(0.675, 0.88, '(C)', fontsize=10)
fig.fig.text(0.1, 0.29, '(D)', fontsize=10)

axs_1D[1].set_xscale('log')
axs_1D[2].set_xscale('log')
axs_1D[1].legend()
axs_1D[0].set_ylabel('Normalised max difference')
fig_1D.text(0.1, 0.925, '(A)', fontsize=10)
fig_1D.text(0.37, 0.925, '(B)', fontsize=10)
fig_1D.text(0.645, 0.925, '(C)', fontsize=10)
for i in range(3):
    axs_1D[i].set_ylim(-0.05, 1.05)
    hist_panel[i].set_ylim(0, max_count + 50)
# fig_1D.savefig(os.path.join(fig_dir, 'SA_1D_new.pdf'),
#                bbox_inches='tight')

hist_panel[0].spines[['right', 'top']].set_visible(False)
for i in range(1, 4):
    hist_panel[i].sharey(hist_panel[0])
    hist_panel[i].tick_params(labelleft=False)
    hist_panel[i].spines[['right', 'top']].set_visible(False)
    hist_panel[i].set_xlabel(r'$\Delta \mathrm{APD}_{90}$')
fig.savefig(os.path.join(fig_dir, 'SA_3D_hist_normsignedRMSD.pdf'),)
            # bbox_inches='tight', pad_inches=0.2)
