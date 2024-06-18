import argparse
import glob
import os
import pandas as pd

import modelling

parser = argparse.ArgumentParser(
    description="Combine all simulated output")
parser.add_argument("APmodel", help="Name of AP model")
args = parser.parse_args()

data_dir = os.path.join(modelling.RESULT_DIR, 'parameter_space_exploration')
file_prefix = 'SA_paramid_'
results_dir = os.path.join(data_dir, 'SA_space', args.APmodel)

result_files = glob.glob(os.path.join(results_dir,
                                      f'{file_prefix}*.csv'))
id_list = []
for paths in result_files:
    id_num = int(os.path.basename(paths)[len(file_prefix):-4])
    id_list.append(id_num)
id_list = sorted(id_list)

filenum_list = list(set([int(i / 1000) for i in id_list]))

for fnum in filenum_list:

    sub_id_list = [i for i in id_list if int(i / 1000) == fnum]

    filename = f'SA_allparam_{fnum}.csv'
    filepath = os.path.join(results_dir, filename)

    if os.path.exists(filepath):
        combined_df = pd.read_csv(filepath, header=[0, 1], index_col=[0],
                                  skipinitialspace=True)
    else:
        fname_id = f'{file_prefix}{sub_id_list[0]}.csv'
        combined_df = pd.read_csv(os.path.join(results_dir, fname_id),
                                  header=[0, 1], index_col=[0],
                                  skipinitialspace=True)
        sub_id_list = sub_id_list[1:]
    for i in sub_id_list:
        output = pd.read_csv(os.path.join(results_dir,
                                          f'{file_prefix}{i}.csv'),
                             header=[0, 1], index_col=[0],
                             skipinitialspace=True)
        combined_df = pd.concat([combined_df, output])
    combined_df.to_csv(filepath)
