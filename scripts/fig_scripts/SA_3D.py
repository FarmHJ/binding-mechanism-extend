import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os
import pandas as pd
import sys

import modelling

APmodel_name = sys.argv[1]

# Define directory to save figure
fig_dir = '../../figures/parameter_exploration/' + APmodel_name + '/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Get list of synthetic drugs' names
param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

# Get the name of parameters
# SA_model = modelling.SensitivityAnalysis()
# param_names = SA_model.param_names

# starting_param_df = pd.DataFrame([1] * 5, index=param_names).T
# ComparisonController = modelling.ModelComparison(starting_param_df)

# # Read simulated data for synthetic drugs
root_dir = '../../simulation_data/'
data_dir = root_dir + 'parameter_SA/'
filename = 'SA_alldrugs_' + APmodel_name + '.csv'
df = pd.read_csv(data_dir + filename, header=[0, 1], index_col=[0],
                 skipinitialspace=True)

Vhalf_list = df['param_values']['Vhalf'].values
Kmax_list = df['param_values']['Kmax'].values
Ku_list = df['param_values']['Ku'].values
drug_list = df['drug']['drug'].values

RMSError_drug = df['RMSE']['RMSE'].values
MError_drug = df['ME']['ME'].values

# Calculate signed RMSD
Error_drug = np.array(RMSError_drug) * np.array(MError_drug) / \
    np.abs(np.array(MError_drug))

# Read simulated data of virtual drugs in the parameter space
data_dir = root_dir + 'parameter_space_exploration/SA_space/' + \
    APmodel_name + '/'
file_prefix = 'SA_allparam'
result_files = [data_dir + f for f in os.listdir(data_dir)
                if f.startswith(file_prefix)]

# Define the range where the RMSD between the APD90s of the two models are
# small
data_dir = root_dir + 'parameter_SA/APD90diff_N/' + APmodel_name + '/'
overall_stats = pd.read_csv(data_dir + 'overall_stats.csv', index_col=[0])
error_range = overall_stats['max'].values[0] * 1.2
print(error_range)

# Load results and extract points where the RMSD value is within the defined
# range
first_iter = True
for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)
    # chosen_df = df.loc[df['RMSE']['RMSE'] < error_range]
    chosen_df = df.loc[df['ME']['ME'] > 0]

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
nan_ind = [i for i in range(len(RMSError)) if np.isnan(RMSError[i]) or
           np.isnan(MError[i])]
Error_space = RMSError * MError / np.abs(MError)

cmin = min(min(Error_drug), min(Error_space))
cmax = max(max(Error_drug), max(Error_space))
# cmin = min(Error_space)
# cmax = max(Error_space)
print('min and max:', cmin, cmax)

Vhalf_range = [Vhalf_range[i] for i in range(len(Vhalf_range))
               if i not in nan_ind]
Kmax_range = [Kmax_range[i] for i in range(len(Kmax_range))
              if i not in nan_ind]
Ku_range = [Ku_range[i] for i in range(len(Ku_range))
            if i not in nan_ind]
Error_space = [Error_space[i] for i in range(len(Error_space))
               if i not in nan_ind]

Vhalf_chosen = combined_chosen_df['param_values']['Vhalf'].values
Kmax_chosen = combined_chosen_df['param_values']['Kmax'].values
Ku_chosen = combined_chosen_df['param_values']['Ku'].values

print(type(Vhalf_chosen))
print(len(Vhalf_chosen))
# # Read simulated data of virtual drugs in the parameter space around the
# # surface where the RMSD value is small
# data_dir = data_dir + '../SA_curve/'
# file_prefix = 'SA_curve'
# result_files2 = [data_dir + f for f in os.listdir(data_dir)
#                  if f.startswith(file_prefix)]

# min_Ku = min(Ku_range)

# # Load results where the RMSD value is within the defined range
# first_iter = True
# for file in result_files2:
#     df = pd.read_csv(file,
#                      header=[0, 1], index_col=[0],
#                      skipinitialspace=True)
#     chosen_df = df.loc[df['RMSE']['RMSE'] < error_range]

#     if first_iter:
#         curve_chosen_df = chosen_df
#         first_iter = False
#     else:
#         curve_chosen_df = pd.concat([curve_chosen_df, chosen_df])

# Vhalf_curve = curve_chosen_df['param_values']['Vhalf'].values
# Kmax_curve = curve_chosen_df['param_values']['Kmax'].values
# Ku_curve = curve_chosen_df['param_values']['Ku'].values


def log_tick_formatter(val, pos=None):
    return f"$10^{{{int(val)}}}$"


# Set up figure structure
fig = plt.figure(figsize=(10, 5))

gs = fig.add_gridspec(1, 2, wspace=0.1)
axs = [fig.add_subplot(gs[0, j], projection='3d') for j in range(2)]
# gs = fig.add_gridspec(1, 1)
# axs = fig.add_subplot(gs[0, 0], projection='3d')

cmap = plt.get_cmap('rainbow')
cmap_norm = matplotlib.colors.Normalize(cmin, cmax)
scale_map = matplotlib.cm.ScalarMappable(norm=cmap_norm, cmap=cmap)

# Plot points in the parameter space
axs[0].scatter(Vhalf_range, np.log10(Kmax_range), np.log10(Ku_range),
               c=scale_map.to_rgba(Error_space),
               s=5, marker='o', zorder=-10, alpha=0.5)
axs[0].view_init(32, 55)

# Plot points of all synthetic drugs and those with RMSD within the defined
# range
xmin, xmax = min(Vhalf_range), max(Vhalf_range)
axs[1].scatter(Vhalf_chosen, np.log10(Kmax_chosen), np.log10(Ku_chosen),
               c='dimgrey', s=10, marker='o', zorder=-10, alpha=0.5)
# axs[1].scatter(Vhalf_curve, np.log10(Kmax_curve), np.log10(Ku_curve),
#                c='dimgrey', s=10, marker='o', zorder=-10, alpha=0.5)
axs[1].scatter(Vhalf_list, np.log10(Kmax_list), np.log10(Ku_list),
               c=scale_map.to_rgba(Error_drug),
               s=100, marker='^', zorder=-1)
axs[1].scatter(xmin * np.ones(len(Vhalf_list)), np.log10(Kmax_list),
               np.log10(Ku_list), s=50, marker='o', zorder=-1, c='red')

for i in range(len(Vhalf_list)):
    axs[1].plot([xmin, Vhalf_list[i]],
                [np.log10(Kmax_list[i]), np.log10(Kmax_list[i])],
                zs=[np.log10(Ku_list[i]), np.log10(Ku_list[i])],
                color='red', linestyle='--', linewidth=0.7)
axs[1].view_init(32, 55)

# Adjust figure details
for i in range(2):
    axs[i].set_xlabel(r"$V_\mathrm{half-trap}$")
    axs[i].set_ylabel(r"$K_\mathrm{max}$")
    axs[i].set_zlabel(r"$K_u$")

    axs[i].set_xlim(min(Vhalf_range), max(Vhalf_range))
    axs[i].set_ylim(min(np.log10(Kmax_range)), max(np.log10(Kmax_range)))
    # axs[i].set_zlim(min(np.log10(Ku_curve)), max(np.log10(Ku_curve)))
    axs[i].set_zlim(min(np.log10(Ku_range)), max(np.log10(Ku_range)))

    axs[i].xaxis.set_major_locator(mticker.MaxNLocator(nbins=6))
    axs[i].yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
    axs[i].yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    axs[i].zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
    axs[i].zaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    axs[i].set_rasterization_zorder(0)

# cax = axs.inset_axes([0.01, -0.08, 1, 0.03])
cax = axs[0].inset_axes([0.5, -0.08, 1, 0.03])
scale_map.set_array(Error_space)
fig.colorbar(scale_map, orientation='horizontal', ax=axs, cax=cax)

# Add panel labels
fig.text(0.075, 0.75, '(A)', fontsize=11)
fig.text(0.5, 0.75, '(B)', fontsize=11)

# Save figure
# plt.subplots_adjust(hspace=0)
plt.savefig(fig_dir + 'SA_3D_positive.png', bbox_inches='tight', pad_inches=0.2)

# #
# # Plot previous figures at different angles (for supplementary materials)
# #
# # Set up figure structure
# fig = plt.figure(figsize=(10, 5))

# gs = fig.add_gridspec(1, 2, wspace=0.1)
# axs = [fig.add_subplot(gs[0, j], projection='3d') for j in range(2)]

# axs[0].scatter(Vhalf_range, np.log10(Kmax_range), np.log10(Ku_range),
#                c=scale_map.to_rgba(Error_space),
#                s=5, marker='o', zorder=-10, alpha=0.5)
# axs[0].view_init(32, 55)

# # Plot points of all synthetic drugs and those with RMSD within the defined
# # range
# xmin, xmax = min(Vhalf_range), max(Vhalf_range)
# axs[1].scatter(Vhalf_chosen, np.log10(Kmax_chosen), np.log10(Ku_chosen),
#                c='dimgrey', s=10, marker='o', zorder=-10, alpha=0.5)
# axs[1].scatter(Vhalf_curve, np.log10(Kmax_curve), np.log10(Ku_curve),
#                c='dimgrey', s=10, marker='o', zorder=-10, alpha=0.5)
# axs[1].scatter(Vhalf_list, np.log10(Kmax_list), np.log10(Ku_list),
#                c=scale_map.to_rgba(Error_drug),
#                s=100, marker='^', zorder=-1)
# axs[1].scatter(xmin * np.ones(len(Vhalf_list)), np.log10(Kmax_list),
#                np.log10(Ku_list), s=50, marker='o', zorder=-1, c='red')

# for i in range(len(Vhalf_list)):
#     axs[1].plot([xmin, Vhalf_list[i]],
#                 [np.log10(Kmax_list[i]), np.log10(Kmax_list[i])],
#                 zs=[np.log10(Ku_list[i]), np.log10(Ku_list[i])],
#                 color='red', linestyle='--', linewidth=0.7)

# # Initiate different viewing angles
# axs[0].view_init(10, 45)
# axs[1].view_init(10, 45)

# # Adjust figure details
# for i in range(2):
#     axs[i].set_xlabel(r"$V_\mathrm{half-trap}$")
#     axs[i].set_ylabel(r"$K_\mathrm{max}$")
#     axs[i].set_zlabel(r"$K_u$")

#     axs[i].set_xlim(min(Vhalf_range), max(Vhalf_range))
#     axs[i].set_ylim(min(np.log10(Kmax_range)), max(np.log10(Kmax_range)))
#     axs[i].set_zlim(min(np.log10(Ku_curve)), max(np.log10(Ku_curve)))

#     axs[i].xaxis.set_major_locator(mticker.MaxNLocator(nbins=6))
#     axs[i].yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
#     axs[i].yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
#     axs[i].zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
#     axs[i].zaxis.set_major_locator(mticker.MaxNLocator(integer=True))

#     axs[i].set_rasterization_zorder(0)

# cax = axs[0].inset_axes([0.5, -0.08, 1, 0.03])
# scale_map.set_array(Error_space)
# fig.colorbar(scale_map, orientation='horizontal', ax=axs, cax=cax)

# fig.text(0.075, 0.75, '(A)', fontsize=11)
# fig.text(0.5, 0.75, '(B)', fontsize=11)

# # Save figure
# plt.subplots_adjust(hspace=0)
# fig_dir = '../../figures/supp_mat/'
# if not os.path.isdir(fig_dir):
#     os.makedirs(fig_dir)
# plt.savefig(fig_dir + 'FigS_SA_3D.png', bbox_inches='tight')
