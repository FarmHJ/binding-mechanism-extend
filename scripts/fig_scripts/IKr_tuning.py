import myokit
import os
import pandas as pd

import modelling


# Define protocol
pulse_time = 1000
protocol_offset = 50
protocol = myokit.pacing.blocktrain(pulse_time, 0.5, offset=protocol_offset)

# Define constants
repeats = 1000
abs_tol = 1e-7
rel_tol = 1e-8

fig_dir = '../../figures/kinetics_comparison/'
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)
data_dir = '../../simulation_data/'

# Set up figure for reversal potential, AP and current contribution
plot = modelling.figures.FigurePlot()
fig = modelling.figures.FigureStructure(figsize=(9, 5), gridspec=(2, 2),
                                        height_ratios=[1] * 2, hspace=0.37,
                                        wspace=0.1, plot_in_subgrid=True)

subgridspecs = [(2, 1)] * 4
subgs = []
for i in range(4):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], hspace=0.1))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for k in
    range(len(subgs))]

# Figure parameters
model_details = modelling.ModelDetails()
current_list = model_details.current_list
plotting_pulse_time = 800
current_colours = model_details.current_colours

model_list = ['Grandi', 'TTP', 'Tomek-Cl', 'Lei']

# Figure y-axes limits
AP_bottom_list = []
AP_top_list = []
IKr_bottom_list = []
IKr_top_list = []

for num, APmodel_name in enumerate(model_list):
    # Load base model
    model_filenames = model_details.file_names[APmodel_name]
    APmodel = '../../' + model_filenames['AP_path']
    APmodel, _, x = myokit.load(APmodel)
    AP_model = modelling.Simulation(APmodel, base_constant=None)
    AP_model.protocol = protocol

    model_keys = model_details.current_keys[APmodel_name]
    Vm_key = model_keys['Vm']
    IKr_key = model_keys['IKr']
    current_head_key = IKr_key[:IKr_key.index('.')]

    # Simulate AP
    base_log = AP_model.model_simulation(repeats, abs_tol=abs_tol,
                                         rel_tol=rel_tol)

    # Plot AP and hERG
    panel = axs[num]
    panel[0][0].plot(base_log.time(), base_log[Vm_key], 'k--',
                     label='base AP model')
    panel[1][0].plot(base_log.time(), base_log[IKr_key], 'k--')

    # Load Grandi-SD model
    APmodel = '../../' + model_filenames['AP_SD_path']
    APmodel, _, x = myokit.load(APmodel)
    AP_model = modelling.Simulation(APmodel, current_head_key=current_head_key)
    AP_model.protocol = protocol

    log = AP_model.model_simulation(repeats,
                                    conductance_name='tune.ikr_rescale',
                                    conductance_value=1,
                                    abs_tol=abs_tol, rel_tol=rel_tol)

    plot.add_single(panel[0][0], log, Vm_key, color='k',
                    label=r'$I_\mathrm{Kr}$ replaced')
    plot.add_single(panel[1][0], log, IKr_key, color='k')

    # Load IKr scale
    scaling_df_filepath = '../../simulation_data/' + APmodel_name + \
        '_conductance_scale.csv'
    scaling_df = pd.read_csv(scaling_df_filepath, index_col=[0])
    if APmodel_name == 'Lei':
        scaling_factor = scaling_df.loc['AP_duration']['conductance scale']
    else:
        scaling_factor = scaling_df.loc['hERG_peak']['conductance scale']
    log = AP_model.model_simulation(repeats,
                                    conductance_name='tune.ikr_rescale',
                                    conductance_value=scaling_factor,
                                    abs_tol=abs_tol, rel_tol=rel_tol)

    plot.add_single(panel[0][0], log, Vm_key, color='r',
                    label=r'$I_\mathrm{Kr}$ replaced $+  I_\mathrm{Kr}$ tuned')
    plot.add_single(panel[1][0], log, IKr_key, color='r')
    panel[0][0].set_title(model_filenames['label'])

    AP_y_bottom, AP_y_top = panel[0][0].get_ylim()
    IKr_y_bottom, IKr_y_top = panel[1][0].get_ylim()

    AP_bottom_list.append(AP_y_bottom)
    AP_top_list.append(AP_y_top)
    IKr_bottom_list.append(IKr_y_bottom)
    IKr_top_list.append(IKr_y_top)

    fig.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
               axs=panel, subgridspec=subgridspecs[num])

# Adjust figures
AP_y_min = min(AP_bottom_list)
AP_y_max = max(AP_top_list)
IKr_y_min = min(IKr_bottom_list)
IKr_y_max = max(IKr_top_list)

axs[0][0][0].legend()
for i in range(4):
    axs[i][0][0].set_ylim(AP_y_min, AP_y_max)
    axs[i][1][0].set_ylim(IKr_y_min, IKr_y_max)
    if i == 0 or i == 2:
        axs[i][0][0].set_ylabel('AP')
        axs[i][1][0].set_ylabel(r"$I_\mathrm{Kr}$")
    else:
        axs[i][0][0].set_yticklabels([])
        axs[i][1][0].set_yticklabels([])

fig.savefig(fig_dir + 'AP_tune_IKr_test.svg')
