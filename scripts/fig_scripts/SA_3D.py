import argparse
import glob
import itertools
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

# Read simulated data for synthetic drugs
data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_SA')
filename = f'SA_alldrugs_{APmodel_name}.csv'
df = pd.read_csv(os.path.join(data_dir, filename), header=[0, 1],
                 index_col=[0], skipinitialspace=True)

APsim = modelling.ModelSimController('Grandi')
FnClass = modelling.ModelComparison(APsim)

Vhalf_list = df['param_values']['Vhalf'].values
Kmax_list = df['param_values']['Kmax'].values
Ku_list = df['param_values']['Ku'].values
drug_list = df.index

RMSError_drug = df['RMSE']['RMSE'].values
MError_drug = df['ME']['ME'].values

# Calculate signed RMSD
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
    chosen_df = df.loc[df['RMSE']['RMSE'] < error_range]
    # chosen_df = df.loc[df['ME']['ME'] > 0]

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

# Remove points where there is numerical issue in the simulation
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
Error_space = RMSError * MError / np.abs(MError)
error_key = ('error', 'signed_RMSE')
combined_df.insert(len(combined_df.columns), error_key, Error_norm_space)

cmin = min(min(Error_norm_drug), min(Error_norm_space))
cmax = max(max(Error_norm_drug), max(Error_norm_space))
print('min and max of APD difference in the parameter space:', cmin, cmax)

Vhalf_chosen = combined_chosen_df['param_values']['Vhalf'].values
Kmax_chosen = combined_chosen_df['param_values']['Kmax'].values
Ku_chosen = combined_chosen_df['param_values']['Ku'].values

print('Total number of points within the range: ', len(Vhalf_chosen))


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

hist_panel.hist(Error_space, 20)
hist_panel.axvline(0, 0, 1, color='grey', ls='--', zorder=-1)
hist_panel.fill_between([-error_range, error_range], 0, 3300,
                        alpha=0.5, color='grey', zorder=-2)
hist_panel.spines[['right', 'top']].set_visible(False)
hist_panel.set_xlabel(r'$\Delta \mathrm{APD}_{90}$')

cmap = plt.get_cmap('viridis')
cmap_norm = matplotlib.colors.Normalize(cmin, cmax)
# abs_max = max(np.abs(cmin), cmax)
# cmap_norm = matplotlib.colors.Normalize(-abs_max, abs_max)
# cmap = plt.get_cmap('RdBu_r')
# cmap_norm = matplotlib.colors.TwoSlopeNorm(vmin=cmin, vcenter=0, vmax=cmax)

# cmap_norm = matplotlib.colors.TwoSlopeNorm(vmin=abs_max * -1, vcenter=0,
#                                            vmax=abs_max)
scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

# Plot points in the parameter space
SA_panel[0].scatter(Vhalf_range, np.log10(Kmax_range), np.log10(Ku_range),
                    c=scale_map.to_rgba(Error_norm_space),
                    s=5, marker='o', zorder=-10, alpha=0.5)
SA_panel[0].view_init(32, 55)
# SA_panel[0].view_init(32, -25)

# Plot points of all synthetic drugs and those with RMSD within the defined
# range
xmin, xmax = min(Vhalf_range), max(Vhalf_range)
SA_panel[1].scatter(Vhalf_chosen, np.log10(Kmax_chosen), np.log10(Ku_chosen),
                    c='dimgrey', s=3, marker='o', zorder=-10, alpha=0.2)
SA_panel[1].scatter(Vhalf_list, np.log10(Kmax_list), np.log10(Ku_list),
                    c=scale_map.to_rgba(Error_norm_drug),
                    s=100, marker='^', zorder=-1)
SA_panel[1].scatter(xmin * np.ones(len(Vhalf_list)), np.log10(Kmax_list),
                    np.log10(Ku_list), s=50, marker='o', zorder=-1, c='red')

for i in range(len(Vhalf_list)):
    SA_panel[1].plot([xmin, Vhalf_list[i]],
                     [np.log10(Kmax_list[i]), np.log10(Kmax_list[i])],
                     zs=[np.log10(Ku_list[i]), np.log10(Ku_list[i])],
                     color='red', linestyle='--', linewidth=0.7)
SA_panel[1].view_init(32, 55)
# SA_panel[1].view_init(32, -25)

# Adjust figure details
for i in range(2):
    SA_panel[i].set_xlabel(r"$V_\mathrm{half-trap}$")
    SA_panel[i].set_ylabel(r"$K_\mathrm{max}$")
    SA_panel[i].set_zlabel(r"$K_u$")

    SA_panel[i].set_xlim(min(Vhalf_range), max(Vhalf_range))
    SA_panel[i].set_ylim(min(np.log10(Kmax_range)), max(np.log10(Kmax_range)))
    SA_panel[i].set_zlim(min(np.log10(Ku_range)), max(np.log10(Ku_range)))

    SA_panel[i].xaxis.set_major_locator(mticker.MaxNLocator(nbins=6))
    SA_panel[i].yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
    SA_panel[i].yaxis.set_major_locator(mticker.MaxNLocator(nbins=8))
    SA_panel[i].zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
    SA_panel[i].zaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    SA_panel[i].set_rasterization_zorder(0)

cax = SA_panel[0].inset_axes([0.5, -0.15, 1, 0.03])
scale_map.set_array(Error_norm_space)
fig.fig.colorbar(scale_map, orientation='horizontal', ax=SA_panel, cax=cax,
                 label=r'$\Delta \widetilde{\mathrm{APD}}_{90}$')

# Add panel labels
fig.fig.text(0.075, 0.765, '(A)', fontsize=10)
fig.fig.text(0.375, 0.765, '(B)', fontsize=10)
fig.fig.text(0.655, 0.765, '(C)', fontsize=10)

# Save figure
# fig.savefig(os.path.join(fig_dir, 'SA_3D_anchor0.png'),)
fig.savefig(os.path.join(fig_dir, 'SA_3D_hist_normsignedRMSD.pdf'),)
            # bbox_inches='tight', pad_inches=0.2)

# fig_1D = plt.figure(figsize=(7, 2))
# gs_1D = fig_1D.add_gridspec(1, 3)
# axs_1D = [fig_1D.add_subplot(gs_1D[0, j]) for j in range(3)]

# Vhalf_values = set(Vhalf_range)
# Kmax_values = set(Kmax_range)
# Ku_values = set(Ku_range)
# param_values = {
#     'Vhalf': Vhalf_values,
#     'Kmax': Kmax_values,
#     'Ku': Ku_values
# }
# param_list = ['Vhalf', 'Kmax', 'Ku']
# param_unit = ['mV', '1', r'ms$^{-1}$']
# for c, p_interest in enumerate(param_list):
#     max_diff = 0
#     param_xy = [i for i in param_list if i != p_interest]
#     p1, p2 = param_xy[0], param_xy[1]
#     diff_dict = {}
#     for p1_value, p2_value in itertools.product(param_values[p1],
#                                                 param_values[p2]):

#         temp = combined_df[(combined_df[('param_values', p1)] == p1_value) &
#                            (combined_df[('param_values', p2)] == p2_value)]
#         diff = max(temp[error_key]) - min(temp[error_key])
#         diff_dict.update({(p1_value, p2_value): diff})

#     # Get max point
#     chosen_pt = max(diff_dict, key=diff_dict.get)
#     max_diff = diff_dict[chosen_pt]
#     chosen_id = combined_df[(combined_df[('param_values', p1)] == chosen_pt[0]) &
#                             (combined_df[('param_values', p2)] == chosen_pt[1])].index
#     chosen_df = combined_df[(combined_df[('param_values', p1)] == chosen_pt[0]) &
#                             (combined_df[('param_values', p2)] == chosen_pt[1])]
#     print('parameter of interest: ', p_interest)
#     print('max diff: ', max_diff)
#     print('param point: ', (p1, p2))
#     print('max point: ', chosen_pt)
#     print('param id: ', chosen_id.values)
#     dir = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration',
#                        'check', APmodel_name)
#     np.savetxt(os.path.join(dir, f'chosen_id_{p_interest}.txt'),
#                chosen_id.values, fmt='%d')

#     temp = combined_df[(combined_df[('param_values', p1)] == chosen_pt[0]) &
#                        (combined_df[('param_values', p2)] == chosen_pt[1])]
#     p_interest_range = temp[('param_values', p_interest)]
#     # p_interest_error = temp[error_key]
#     p_interest_error = temp[('ME', 'ME')]
#     # print(p_interest_range.columns)
#     # print('param values: ', p_interest_range.sort_values())

#     error_norm = (p_interest_error - cmin) / (cmax - cmin)
#     print(p_interest_error)
#     axs_1D[c].plot(p_interest_range, p_interest_error, 'o', label=p_interest)
#     # axs_1D[c].set_title(p_interest)
#     axs_1D[c].set_xlabel(f'{p_interest} ({param_unit[c]})')

#     if p_interest == 'Vhalf':
#         step = np.mean(p_interest_range[1:] - p_interest_range[:-1])
#         p_interest_range = list(p_interest_range.values)
#         p_interest_range += [p_interest_range[0] - step,
#                              p_interest_range[-1] + step]
#     else:
#         step = np.mean(np.log10(p_interest_range[1:]) -
#                        np.log10(p_interest_range[:-1]))
#         p_interest_range = list(p_interest_range.values)
#         p_interest_range += [np.exp(np.log(p_interest_range[0]) - step),
#                              np.exp(np.log(p_interest_range[-1]) + step)]

# axs_1D[1].set_xscale('log')
# axs_1D[2].set_xscale('log')
# axs_1D[1].legend()
# axs_1D[0].set_ylabel('Normalised max difference')
# fig_1D.text(0.1, 0.925, '(A)', fontsize=10)
# fig_1D.text(0.37, 0.925, '(B)', fontsize=10)
# fig_1D.text(0.645, 0.925, '(C)', fontsize=10)
# # for i in range(3):
# #     axs_1D[i].set_ylim(-0.05, 1.05)
# fig_1D.savefig(os.path.join(fig_dir, 'SA_1D_normsignedRMSD.pdf'),
#                bbox_inches='tight')
