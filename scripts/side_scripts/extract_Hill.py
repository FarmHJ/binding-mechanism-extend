import glob
import os
import pandas as pd

import modelling

data_fpath = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration',
                          'SA_space_IKr')
result_files = glob.glob(os.path.join(data_fpath, 'SA_allparam*'))

first_iter = True
for file in result_files:
    df = pd.read_csv(file, header=[0, 1], index_col=[0],
                     skipinitialspace=True)
    df = df[['param_id', 'drug_conc_Hill', 'Hill_curve', 'param_values']]

    if first_iter:
        combined_df = df
        first_iter = False
    else:
        combined_df = pd.concat([combined_df, df])

combined_df.to_csv(os.path.join(data_fpath, 'Hill_curves.csv'))
