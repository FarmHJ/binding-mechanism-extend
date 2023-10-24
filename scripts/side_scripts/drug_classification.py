import numpy as np
import os
import pandas as pd

import modelling


fig_dir = os.path.join(modelling.FIG_DIR, 'parameter_exploration')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

model_list = modelling.model_naming.APmodel_list[1:]
drug_similar_count_arr = []
for APmodel_name in model_list:

    # Read simulated data for synthetic drugs
    fpath = os.path.join(modelling.RESULT_DIR, 'parameter_SA',
                         'SA_alldrugs_' + APmodel_name + '.csv')
    df = pd.read_csv(fpath, header=[0, 1], index_col=[0],
                     skipinitialspace=True)

    drug_list = df['drug']['drug'].values
    RMSError_drug = df['RMSE']['RMSE'].values

    # Define the range where the RMSD between the APD90s of the two models are
    # small
    fpath = os.path.join(modelling.RESULT_DIR, 'parameter_SA',
                         'APD90diff_N', APmodel_name, 'overall_stats.csv')
    overall_stats = pd.read_csv(fpath, index_col=[0])
    error_range = overall_stats['max'].values[0] * 1.2

    drug_similar_count = np.sum(RMSError_drug < error_range)
    drug_similar_count_arr.append(drug_similar_count)

fig = modelling.figures.FigureStructure(figsize=(4, 3))
plot = modelling.figures.FigurePlot()

model_number_arr = np.arange(len(drug_similar_count_arr))
fig.axs[0][0].bar(model_number_arr, drug_similar_count_arr, 0.8)
fig.axs[0][0].set_xticks(model_number_arr, model_list)
fig.axs[0][0].set_ylabel('Number of drugs')
fig.axs[0][0].set_title("Drug within similarity range")
fig.savefig(os.path.join(fig_dir, "models_similar_drugs.pdf"))
