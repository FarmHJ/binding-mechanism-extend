import numpy as np
import os
import pandas as pd
import plotly.graph_objects as go

import modelling

APmodel_name = 'Grandi'

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

# Load results and extract points where the RMSD value is within the defined
# range
first_iter = True
for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
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

# Remove points where there is numerical issue in the simulation
nan_ind = [i for i in range(len(RMSError)) if np.isnan(RMSError[i]) or
           np.isnan(MError[i])]
Error_space = RMSError * MError / np.abs(MError)

cmin = min(min(Error_drug), min(Error_space))
cmax = max(max(Error_drug), max(Error_space))

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

fig = go.Figure()

hovertext = np.empty(shape=(13, 3, 1), dtype='object')
hovertext[:, 0] = np.array(drug_list).reshape(-1, 1)
hovertext[:, 1] = np.array(RMSError_drug).reshape(-1, 1)

fig.add_trace(
    go.Scatter3d(
        x=Vhalf_list,
        y=Kmax_list,
        z=Ku_list,
        mode='markers',
        marker_symbol='diamond',
        name='',
        customdata=hovertext,
        hovertemplate='<b>%{customdata[0]}</b> <br>RMSD = %{customdata[1]:.2e}',
        # <br>MD = %{customdata[2]:.2e}',
        marker=dict(
            color=RMSError_drug,
            colorscale='Portland',
            colorbar=dict(thickness=20),
            cmin=cmin,
            cmax=cmax
        )
    )
)

# bin_arr = np.arange(cmin, cmax, 30)
bin_arr = np.linspace(cmin, cmax, 20)
bins = [(bin_arr[i], bin_arr[i + 1]) for i in range(len(bin_arr) - 1)]

for lb, ub in bins:
    chosen_ind = [i for i, e in enumerate(Error_space) if e < ub and e > lb]
    Vhalf_chosen = np.array([Vhalf_range[i] for i in chosen_ind])
    Kmax_chosen = np.array([Kmax_range[i] for i in chosen_ind])
    Ku_chosen = np.array([Ku_range[i] for i in chosen_ind])
    RMSE_chosen = np.array([Error_space[i] for i in chosen_ind])
    # paramid_chosen = np.array([param_id[i] for i in chosen_ind])
#         MAE_chosen = np.array([MAError[i] for i in chosen_ind])

    Vhalf_bg = np.array([Vhalf_range[i] for i in range(len(Error_space))
                         if i not in chosen_ind])
    Kmax_bg = np.array([Kmax_range[i] for i in range(len(Error_space))
                        if i not in chosen_ind])
    Ku_bg = np.array([Ku_range[i] for i in range(len(Error_space))
                      if i not in chosen_ind])
    # paramid_bg = np.array([param_id[i] for i in range(len(Error_space))
    #                        if i not in chosen_ind])

    # hovertext = np.empty(shape=(len(paramid_chosen),3,1), dtype='object')
    hovertext = np.empty(shape=(len(RMSE_chosen), 3, 1), dtype='object')
    # hovertext[:,0] = np.array(paramid_chosen).reshape(-1,1)
    hovertext[:, 0] = np.array(RMSE_chosen).reshape(-1, 1)
#         hovertext[:,2] = np.array(MAE_chosen).reshape(-1,1)

    fig.add_trace(
        go.Scatter3d(
            visible=True,
            x=Vhalf_chosen,
            y=Kmax_chosen,
            z=Ku_chosen,
            mode='markers',
            name='',
            customdata=hovertext,
            hovertemplate='<b>id: %{customdata[0]}</b> <br>RMSD = %{customdata[1]}',  # <br>MD = %{customdata[2]}',
            marker=dict(
                color=RMSE_chosen,
                colorscale='Portland',
                opacity=0.7,
                size=3,
                colorbar=dict(thickness=20),
                cmin=cmin,
                cmax=cmax
            )
        )
    )

sets = []
for i in range(int((len(fig.data) - 1))):
    param_set = dict(
        method="update",
        args=[{"visible": [False] * len(fig.data)},
              {"title": "APD90 difference at range " + "%d" % bins[i][0] +
               " to " + "%d" % bins[i][1]}],
        label="(" + "%d" % bins[i][0] + ", " + "%d" % bins[i][1] + ")"
    )
    param_set["args"][0]["visible"][0] = True
    param_set["args"][0]["visible"][i + 1] = True
#     param_set["args"][0]["visible"][2 * i + 2] = True
    sets.append(param_set)

sliders = [dict(
    active=5,
    currentvalue={"prefix": "Bins at "},
    pad={"t": 5},
    steps=sets
)]

fig.update_layout(sliders=sliders,
                  scene=dict(
                      xaxis_title='Vhalf',
                      yaxis_title='Kmax',
                      zaxis_title='Ku',
                      xaxis=dict(range=[min(Vhalf_range), max(Vhalf_range)]),
                      yaxis=dict(dtick=1,
                                 type='log',
                                 range=[np.log10(min(Kmax_range)),
                                        np.log10(max(Kmax_range))]),
                      zaxis=dict(dtick=1,
                                 type='log',
                                 range=[np.log10(min(Ku_range)),
                                        np.log10(max(Ku_range))])),
                  scene_aspectmode='manual',
                  scene_aspectratio=dict(x=1, y=1.2, z=1))
#                   margin=dict(r=20, l=10, b=10, t=10))

fig.show()
