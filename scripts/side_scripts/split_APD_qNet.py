import os
import pandas as pd

import modelling


for drug in ['dofetilide', 'verapamil']:
    data_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                            'ORd-Li', 'hERG_peak_match', drug)

    for model in ['SD', 'CS']:
        APD_df = pd.read_csv(os.path.join(data_dir, model + '_APD_fine.csv'))
        qNet_df = APD_df[['drug concentration', 'qNet']]
        qNet_df.to_csv(os.path.join(data_dir, model + '_qNet.csv'))
