import glob
import numpy as np
import os
import pandas as pd

import modelling


###############################
# Li-SD parameters to csv file
###############################
# Save parameters of SD model to csv file
# binding_params = modelling.BindingParameters().binding_parameters

# df = pd.DataFrame.from_dict(binding_params).T

# dir = '../../parameters/'
# df.to_csv(dir + 'Li-SD.csv')

######################################################
# Hill coefficients from Li et al. (2017) to csv file
######################################################
# # Save Hill coefficients (from Li et al. 2017) to csv file
# Hill_curves = modelling.BindingParameters().Hill_curve

# Hill_param_names = ['Hill_coef', 'IC50']

# drug = list(Hill_curves.keys())[0]
# channels = Hill_curves[drug].keys()

# index_dict = {}
# for c in channels:
#     index_dict[c] = Hill_param_names

# all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
# index = pd.MultiIndex.from_tuples(all_index)

# combined_df = pd.DataFrame(index=index).T
# for d in Hill_curves.keys():
#     Hill_coefs = []
#     for c in channels:
#         Hill_coefs += list(Hill_curves[d][c].values())
#     df = pd.DataFrame(Hill_coefs, index=index).T
#     df = df.rename(index={0: d})
#     combined_df = pd.concat([combined_df, df])

# dir = '../../parameters/'
# combined_df.to_csv(dir + 'Li-Hill.csv')

################################
# Lei-SD parameters to csv file
################################
# param_dir = os.path.join(modelling.MAIN_DIR, '..', 'Lei-SD-fits',
#                          'Milnes-data-fits')

# drug_param_dict = {}
# for fpath in glob.glob(os.path.join(
#         param_dir, 'Drug-*-lei-model12-protocol-Milnes-fit-RMSE.txt')):

#     fname = os.path.basename(fpath)
#     starter = 0
#     occ_count = 0
#     pos_ind = []
#     while occ_count < 3:
#         for i in range(len(fname)):
#             j = fname.find('-', starter)
#             if j != -1:
#                 if occ_count == 0:
#                     pos_ind.append(j + 1)
#                 else:
#                     pos_ind.append(j)
#                 starter = j + 1
#                 occ_count += 1

#     drug_name = fname[pos_ind[0]:pos_ind[1]]

#     drug_params = np.loadtxt(fpath, unpack=True)
#     param_name = ['Kmax', 'Ku', 'halfmax', 'n', 'Vhalf']
#     param_dict = {}
#     for i in range(len(drug_params)):
#         param_dict[param_name[i]] = drug_params[i]
#     # param_df = pd.DataFrame(drug_params,
#     #                         index=['Kmax', 'Ku', 'halfmax', 'n', 'Vhalf']).T
#     param_dict['Vhalf'] = -1 * param_dict['Vhalf']
#     drug_param_dict[drug_name] = param_dict

# df = pd.DataFrame.from_dict(drug_param_dict).T

# dir = '../../parameters/'
# df.to_csv(dir + 'Lei-SD.csv')

###########################################################
# Change style of saved Hill coefficients from Li-SD model
###########################################################
for drug in ['dofetilide', 'verapamil']:
    Hill_coef_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                                'Hill_curves', drug)
    Hill_coef_df = pd.read_csv(os.path.join(Hill_coef_dir, 'Hill_curves.csv'),
                               index_col=[0], skipinitialspace=True).T
    Hill_coef_df.to_csv(os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                                     'Hill_curves', drug, 'Lei_Hill.csv'))