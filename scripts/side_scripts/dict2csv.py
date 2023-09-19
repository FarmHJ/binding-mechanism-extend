import pandas as pd

import modelling


# Save parameters of SD model to csv file
# binding_params = modelling.BindingParameters().binding_parameters

# df = pd.DataFrame.from_dict(binding_params).T

# dir = '../../parameters/'
# df.to_csv(dir + 'Li-SD.csv')

# Save Hill coefficients (from Li et. al. 2017) to csv file
Hill_curves = modelling.BindingParameters().Hill_curve

Hill_param_names = ['Hill_coef', 'IC50']

drug = list(Hill_curves.keys())[0]
channels = Hill_curves[drug].keys()

index_dict = {}
for c in channels:
    index_dict[c] = Hill_param_names

all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
index = pd.MultiIndex.from_tuples(all_index)

combined_df = pd.DataFrame(index=index).T
for d in Hill_curves.keys():
    Hill_coefs = []
    for c in channels:
        Hill_coefs += list(Hill_curves[d][c].values())
    df = pd.DataFrame(Hill_coefs, index=index).T
    df = df.rename(index={0: d})
    combined_df = pd.concat([combined_df, df])

dir = '../../parameters/'
combined_df.to_csv(dir + 'Li-Hill.csv')
