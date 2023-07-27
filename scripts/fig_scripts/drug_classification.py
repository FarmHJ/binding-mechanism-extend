import numpy as np
import os
import pandas as pd

import modelling


fig_dir = '../../figures/parameter_exploration/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)
root_dir = '../../simulation_data/'

# Get list of synthetic drugs' names
param_lib = modelling.BindingParameters()
drug_list = param_lib.drug_compounds

model_details = modelling.ModelDetails()
model_name_list = model_details.APmodels[1:]
drug_similar_count_arr = []
for APmodel_name in model_name_list:

    # Read simulated data for synthetic drugs
    data_dir = root_dir + 'parameter_SA/'
    filename = 'SA_alldrugs_' + APmodel_name + '.csv'
    df = pd.read_csv(data_dir + filename, header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    drug_list = df['drug']['drug'].values
    RMSError_drug = df['RMSE']['RMSE'].values

    # Define the range where the RMSD between the APD90s of the two models are
    # small
    data_dir = root_dir + 'parameter_SA/APD90diff_N/' + APmodel_name + '/'
    overall_stats = pd.read_csv(data_dir + 'overall_stats.csv', index_col=[0])
    error_range = overall_stats['max'].values[0] * 1.1

    drug_similar_count = np.sum(RMSError_drug < error_range)
    drug_similar_count_arr.append(drug_similar_count)

fig = modelling.figures.FigureStructure(figsize=(4, 3))
plot = modelling.figures.FigurePlot()

model_number_arr = np.arange(len(drug_similar_count_arr))
fig.axs[0][0].bar(model_number_arr,
                  drug_similar_count_arr, 0.8)
fig.axs[0][0].set_xticks(model_number_arr, model_name_list)
fig.axs[0][0].set_ylabel('Number of drugs')
fig.axs[0][0].set_title("Drug within similarity range")
fig.savefig(fig_dir + "models_similar_drugs.pdf")
