import os

import modelling


# Set up figure for reversal potential, AP and current contribution
plot = modelling.figures.FigurePlot()
fig_AP = modelling.figures.FigureStructure(figsize=(9, 3), gridspec=(1, 4),
                                           wspace=0.1, plot_in_subgrid=True)
subgridspecs = [(2, 1)] * 4
fig_AP.subgrid(subgridspecs, hspace=0.1)

fig_IKr = modelling.figures.FigureStructure(figsize=(9, 1.5), gridspec=(1, 4),
                                            wspace=0.1)

prepace = 1000
# Figure parameters
model_details = modelling.model_naming
model_list = model_details.APmodel_list[1:]
model_label = model_details.AP_file_names
plotting_pulse_time = 600

APsim = modelling.ModelSimController('ORd-Li')
ORdLi_log = APsim.simulate(prepace=prepace)

# Figure y-axes limits
AP_bottom_list = []
AP_top_list = []
IKr_bottom_list = []
IKr_top_list = []

for num, APmodel_name in enumerate(model_list):
    print('Simulating AP for', APmodel_name)
    # Load base model
    base_APsim = modelling.ModelSimController(APmodel_name, ikr_modified=False)
    model_keys = model_details.model_current_keys[APmodel_name]

    # Simulate AP
    base_log = base_APsim.simulate(prepace=prepace)

    # Plot AP and hERG
    panel = fig_AP.axs[num]
    panel[0][0].plot(base_log.time(), base_log[base_APsim.Vm_key], 'k--')
    panel[1][0].plot(base_log.time(), base_log[base_APsim.ikr_key], 'k--',
                     label='base AP model')

    # Load AP-IKr model
    APsim = modelling.ModelSimController(APmodel_name)
    log = APsim.simulate(prepace=prepace)

    panel[0][0].plot(log.time(), log[APsim.Vm_key], 'k')
    panel[1][0].plot(log.time(), log[APsim.ikr_key], 'k',
                     label=r'$I_\mathrm{Kr}$ replaced')

    APsim.set_ikr_rescale_method('AP_duration')
    log = APsim.simulate(prepace=prepace)

    panel[0][0].plot(log.time(), log[APsim.Vm_key], 'r')
    panel[1][0].plot(log.time(), log[APsim.ikr_key], 'r',
                     label=r'$I_\mathrm{Kr}$ replaced'
                     r'$+  I_\mathrm{Kr}$ tuned')
    panel[0][0].set_title(model_label[APmodel_name]['label'])

    AP_y_bottom, AP_y_top = panel[0][0].get_ylim()
    IKr_y_bottom, IKr_y_top = panel[1][0].get_ylim()
    panel[0][0].spines[['right', 'top']].set_visible(False)
    panel[1][0].spines[['right', 'top']].set_visible(False)

    AP_bottom_list.append(AP_y_bottom)
    AP_top_list.append(AP_y_top)
    IKr_bottom_list.append(IKr_y_bottom)
    IKr_top_list.append(IKr_y_top)

    fig_AP.sharex(['Time (ms)'], [(0, plotting_pulse_time)],
                  axs=panel, subgridspec=subgridspecs[num])

    # Plot IKrs
    for i in range(4):
        if i == num:
            fig_IKr.axs[0][i].plot(log.time(), log[APsim.ikr_key], 'k',
                                   zorder=5)
            fig_IKr.axs[0][i].set_title(model_label[APmodel_name]['label'])
        else:
            fig_IKr.axs[0][i].plot(log.time(), log[APsim.ikr_key], '#cccccc',
                                   alpha=0.9)
            fig_IKr.axs[0][i].plot(ORdLi_log.time(), ORdLi_log['ikr.IKr'],
                                   '#cccccc', linestyle='--', alpha=0.9)

# Adjust figures
AP_y_min = min(AP_bottom_list)
AP_y_max = max(AP_top_list)
IKr_y_min = min(IKr_bottom_list)
IKr_y_max = max(IKr_top_list)

fig_AP.axs[0][1][0].legend(handlelength=1.5, loc='upper left',
                           handletextpad=0.5, bbox_to_anchor=(0, 1.05))
for i in range(4):
    fig_AP.axs[i][0][0].set_ylim(AP_y_min, AP_y_max)
    fig_AP.axs[i][1][0].set_ylim(IKr_y_min, IKr_y_max)
    if i == 0:
        fig_AP.axs[i][0][0].set_ylabel('Voltage (mV)')
        fig_AP.axs[i][1][0].set_ylabel("Current (A/F)")
    else:
        fig_AP.axs[i][0][0].set_yticklabels([])
        fig_AP.axs[i][1][0].set_yticklabels([])

fig_AP.savefig(os.path.join(modelling.FIG_DIR, 'AP_tune_IKr.pdf'))

fig_IKr.sharex(['Time (ms)'] * 4, [(0, plotting_pulse_time)] * 4)
fig_IKr.sharey(["Current (A/F)"])
fig_IKr.savefig(os.path.join(modelling.FIG_DIR, 'IKr_magnitude.pdf'))
