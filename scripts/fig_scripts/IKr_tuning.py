import os

import modelling


fig_dir = os.path.join(modelling.FIG_DIR, 'kinetics_comparison')
if not os.path.isdir(fig_dir):
    os.makedirs(fig_dir)

# Set up figure for reversal potential, AP and current contribution
plot = modelling.figures.FigurePlot()
fig = modelling.figures.FigureStructure(figsize=(9, 5), gridspec=(2, 2),
                                        height_ratios=[1] * 2, hspace=0.4,
                                        wspace=0.1, plot_in_subgrid=True)

subgridspecs = [(2, 1)] * 4
subgs = []
for i in range(4):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], hspace=0.1))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for k in
    range(len(subgs))]

# Figure parameters
current_list = modelling.model_naming.current_list
plotting_pulse_time = 800

model_list = ['Grandi', 'TTP', 'Tomek-Cl', 'ORd-Lei']
model_label = ['Grandi-Li', 'ten Tusscher-Li', 'Tomek-Li', 'ORd-Lei']

# Figure y-axes limits
AP_bottom_list = []
AP_top_list = []
IKr_bottom_list = []
IKr_top_list = []

for num, APmodel_name in enumerate(model_list):
    print(APmodel_name)
    # Load base model
    APsim = modelling.ModelSimController(APmodel_name, ikr_modified=False)

    model_keys = modelling.model_naming.model_current_keys[APmodel_name]

    # Simulate AP
    base_log = APsim.simulate()

    # Plot AP and hERG
    panel = axs[num]
    panel[0][0].plot(base_log.time(), base_log[APsim.Vm_key], 'k--',
                     label='base AP model')
    panel[1][0].plot(base_log.time(), base_log[APsim.ikr_key], 'k--')

    # Load AP-IKr model
    APsim = modelling.ModelSimController(APmodel_name)
    log = APsim.simulate()

    panel[0][0].plot(log.time(), log[APsim.Vm_key], 'k',
                     label=r'$I_\mathrm{Kr}$ replaced')
    panel[1][0].plot(log.time(), log[APsim.ikr_key], 'k')

    APsim.set_ikr_rescale_method('AP_duration')
    log = APsim.simulate()

    panel[0][0].plot(log.time(), log[APsim.Vm_key], 'r',
                     label=r'$I_\mathrm{Kr}$ replaced '
                     r'$+  I_\mathrm{Kr}$ tuned')
    panel[1][0].plot(log.time(), log[APsim.ikr_key], 'r')
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

axs[0][0][0].legend()
for i in range(4):
    axs[i][0][0].set_ylim(AP_y_min, AP_y_max)
    axs[i][1][0].set_ylim(IKr_y_min, IKr_y_max)
    if i == 0 or i == 2:
        axs[i][0][0].set_ylabel('Voltage (mV)')
        axs[i][1][0].set_ylabel("Current (A/F)")
    else:
        axs[i][0][0].set_yticklabels([])
        axs[i][1][0].set_yticklabels([])

fig.savefig(os.path.join(fig_dir, 'AP_tune_IKr_test.svg'))

#########################################
# IKr magnitude comparison
#########################################
model_list = modelling.model_naming.APmodel_list[:4]
model_label = ['ORd-Li', 'Grandi-Li', 'ten Tusscher-Li', 'Tomek-Li']

# Set up figure for reversal potential, AP and current contribution
plot = modelling.figures.FigurePlot()
fig = modelling.figures.FigureStructure(figsize=(7, 2.5), gridspec=(2, 2),
                                        height_ratios=[1] * 2, hspace=0.3,
                                        wspace=0.1)

# Figure y-axes limits
IKr_bottom_list = []
IKr_top_list = []

for num, APmodel_name in enumerate(model_list):
    # Load base model
    APsim = modelling.ModelSimController(APmodel_name)

    # Load IKr scale
    if APmodel_name != 'ORd-Li':
        APsim.set_ikr_rescale_method('hERG_peak')
    log = APsim.simulate()

    for i in range(4):
        if i == num:
            r, c = int(num / 2), num % 2
            fig.axs[r][c].plot(log.time(), log[APsim.ikr_key], 'k', zorder=5)
            fig.axs[r][c].set_title(model_label[num])
        else:
            r, c = int(i / 2), i % 2
            fig.axs[r][c].plot(log.time(), log[APsim.ikr_key], '#cccccc',
                               alpha=0.9)

    IKr_y_bottom, IKr_y_top = fig.axs[r][c].get_ylim()
    IKr_bottom_list.append(IKr_y_bottom)
    IKr_top_list.append(IKr_y_top)

    fig.sharex(['Time (ms)'] * 2, [(0, plotting_pulse_time)] * 2)

# Adjust figures
IKr_y_min = min(IKr_bottom_list)
IKr_y_max = max(IKr_top_list)

for i in range(4):
    r, c = int(i / 2), i % 2
    fig.axs[r][c].set_ylim(IKr_y_min, IKr_y_max)
    if i == 0 or i == 2:
        fig.axs[r][c].set_ylabel("Current (A/F)")
    else:
        fig.axs[r][c].set_yticklabels([])

fig.savefig(os.path.join(fig_dir, 'IKr_magnitude.svg'))
