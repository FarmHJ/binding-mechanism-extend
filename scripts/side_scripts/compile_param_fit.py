import os
import pandas as pd

import modelling


# Define directory where optimised parameters are saved
results_dir = os.path.join(modelling.PARAM_DIR, 'Lei-SD-inference')

compile_param = pd.DataFrame()
drug_list = modelling.SD_details.drug_names
for drug in drug_list:

    # Load inference results of given drug
    fname = os.path.join(results_dir, f'{drug}.csv')
    scores_dict = pd.read_csv(fname, index_col=0)
    scores_dict = scores_dict.sort_values(by=['error'])

    # Take best performing parameter set
    p = scores_dict.iloc[[0], :]
    p = p.reset_index(drop=True)
    p = p.rename(index={0: drug})

    # Compile best performing parameter into a file
    compile_param = pd.concat([compile_param, p])

# Save file
compile_param.to_csv(os.path.join(modelling.PARAM_DIR, 'Lei-SD.csv'))
