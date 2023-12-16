import os
import pandas as pd

import modelling


drug_list = modelling.SD_details.drug_names
results_dir = os.path.join(modelling.PARAM_DIR, 'Lei-SD-inference')

compile_param = pd.DataFrame()

for drug in drug_list:

    fname = os.path.join(results_dir, f'{drug}-Milnes.csv')
    scores_dict = pd.read_csv(fname, index_col=0)
    scores_dict = scores_dict.sort_values(by=['error'])
    # scores_dict.to_csv(fname)

    p = scores_dict.iloc[[0], :]
    p = p.reset_index(drop=True)
    p = p.rename(index={0: drug})

    compile_param = pd.concat([compile_param, p])

compile_param.to_csv(os.path.join(modelling.PARAM_DIR, 'Farm-SD.csv'))
