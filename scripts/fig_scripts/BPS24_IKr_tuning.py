import os

import modelling


fig_dir = os.path.join(modelling.FIG_DIR, 'kinetics_comparison')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Set up figure for reversal potential, AP and current contribution
plot = modelling.figures.FigurePlot()
fig = modelling.figures.FigureStructure(figsize=(9, 2.5), gridspec=(1, 3),
                                        wspace=0.1, plot_in_subgrid=True)

subgridspecs = [(2, 1)] * 3
fig.subgrid(subgridspecs, hspace=0.1)

# Figure parameters
model_details = modelling.model_naming
current_list = model_details.current_list
plotting_pulse_time = 450

model_list = model_details.APmodel_list[1:-1]
model_label = ['Grandi-Li', 'ten Tusscher-Li', 'Tomek-Li']

# Figure y-axes limits
AP_bottom_list = []
AP_top_list = []
IKr_bottom_list = []
IKr_top_list = []

for num, APmodel_name in enumerate(model_list):
    print('Simulating AP for', APmodel_name)
    # Load base model
    APsim = modelling.ModelSimController(APmodel_name, ikr_modified=False)

    model_keys = model_details.model_current_keys[APmodel_name]

    # Simulate AP
    base_log = APsim.simulate(prepace=5)

    # Plot AP and hERG
    panel = fig.axs[num]
    panel[0][0].plot(base_log.time(), base_log[APsim.Vm_key], 'k--',
                     label='AP model')
    panel[1][0].plot(base_log.time(), base_log[APsim.ikr_key], 'k--',
                     label='AP model')

    # Load AP-IKr model
    APsim = modelling.ModelSimController(APmodel_name)
    log = APsim.simulate(prepace=5)

    panel[0][0].plot(log.time(), log[APsim.Vm_key], 'k',
                     label='AP-Li model')
    panel[1][0].plot(log.time(), log[APsim.ikr_key], 'k',
                     label='AP-Li model')

    APsim.set_ikr_rescale_method('AP_duration')
    log = APsim.simulate(prepace=5)

    panel[0][0].plot(log.time(), log[APsim.Vm_key], 'r',
                     label='calibrated \nAP-Li model')
    panel[1][0].plot(log.time(), log[APsim.ikr_key], 'r',
                     label='calibrated \nAP-Li model')
    panel[0][0].set_title(model_label[num])

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

fig.axs[0][1][0].legend(handlelength=1.5)
for i in range(3):
    fig.axs[i][0][0].set_ylim(AP_y_min, AP_y_max)
    fig.axs[i][1][0].set_ylim(IKr_y_min, IKr_y_max)
    if i == 0:
        fig.axs[i][0][0].set_ylabel('Voltage (mV)')
        fig.axs[i][1][0].set_ylabel("Current (A/F)")
    else:
        fig.axs[i][0][0].set_yticklabels([])
        fig.axs[i][1][0].set_yticklabels([])

fig.savefig(os.path.join(fig_dir, 'AP_tune_IKr_BPS24.pdf'))
