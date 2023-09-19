import os
import pandas as pd

data_filepath = '../simulation_data/parameter_space_exploration/SA_space_IKr/'
file_prefix = 'SA_allparam'
result_files = [data_filepath + f for f in os.listdir(data_filepath)
                if f.startswith(file_prefix)]

first_iter = True
for file in result_files:
    df = pd.read_csv(file,
                     header=[0, 1], index_col=[0],
                     skipinitialspace=True)
    df = df[['param_id', 'drug_conc_Hill', 'Hill_curve', 'param_values']]

    if first_iter:
        combined_df = df
        first_iter = False
    else:
        combined_df = pd.concat([combined_df, df])

# print(combined_df.head())
combined_df.to_csv(data_filepath + 'Hill_curves.csv')