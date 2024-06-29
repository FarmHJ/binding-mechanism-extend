import matplotlib
import matplotlib.pyplot as plt
import myokit
import numpy as np
import os

import modelling

fig = modelling.figures.FigureStructure(figsize=(10, 7), gridspec=(3, 2),
                                        height_ratios=[2, 3, 3],
                                        hspace=0.53, wspace=0.3,
                                        plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

ikr_gridspec = (2, 1)
state_gridspec = [(1, 2), (1, 2)]
subgs = []
subgs.append(fig.gs[2:4].subgridspec(*ikr_gridspec, hspace=0.08,
                                     height_ratios=[1, 1]))
for i in range(2):
    subgs.append(fig.gs[i + 4].subgridspec(*state_gridspec[i], wspace=0.07))
ikr_axs = [fig.fig.add_subplot(subgs[0][i, 0]) for i in range(2)]
state_axs = [[[fig.fig.add_subplot(subgs[k + 1][i, j]) for j in range(
    state_gridspec[k][1])] for i in range(state_gridspec[k][0])] for
    k in range(2)]

data_dir = os.path.join(modelling.RESULT_DIR, 'background')
fig_dir = os.path.join(modelling.FIG_DIR, 'background')
APmodels = ['ORd-Li', 'ORd-Lei']
current_keys = modelling.model_naming.model_current_keys

# Plot state occupancies of the hERG channel
# color_seq = ['#ff0000', '#ff5555', '#666666', '#999999', '#808080',
#              '#333333', '#4d4d4d']
color_seq = ['#00df00', '#50c850', '#666666', '#999999', '#808080',
             '#333333', '#4d4d4d']
model_color = plt.get_cmap('Dark2')

# 80% gray: 333333
# 70% gray: 4d4d4d
# 60% gray: 666666
# 50% gray: 808080
# 40% gray: 999999

li_states = ['O', 'IO', 'IC1', 'C1', 'IC2', 'C2']
# state 'y3' is open state
# lei_states = ['y3', 'y2', 'y1', 'y4']
lei_states = ['O', 'I', 'CI', 'C']

##################
# Plot AP and IKr
##################
for m, APmodel in enumerate(APmodels):
    print('AP of model: ', APmodel)
    # Plot AP-clamp protocol
    prot_log = myokit.DataLog.load_csv(os.path.join(data_dir,
                                                    f'APclamp_{APmodel}.csv'))
    ikr_axs[0].plot(prot_log.time(), prot_log[current_keys[APmodel]['Vm']],
                    color=model_color(m), label=APmodel)

    m_row, m_col = int(m / 2), m % 2

    # Plot IKr - Li model
    li_log = myokit.DataLog.load_csv(os.path.join(
        data_dir, f'{APmodel}_APclamp_Li_current.csv'))
    ikr_axs[1].plot(li_log.time(), li_log['ikr.IKr'], color=model_color(m))

    # Compute area under the curve of both O and IO states
    auc_li = np.trapz(np.array(li_log['ikr.O']) + np.array(li_log['ikr.IO']),
                      x=li_log.time())
    print('Li: ', auc_li)

    # Plot state occupancy - Li model
    state_panel = state_axs[0][m_row][m_col]
    state_panel.stackplot(
        li_log.time(), *[li_log['ikr.' + s] for s in li_states],
        labels=[s for s in li_states], colors=color_seq[:6], zorder=-10,
        edgecolor='k')
    state_panel.set_ylim(bottom=0, top=1)
    state_panel.text(0.95, 0.95, APmodel, fontsize=8, ha='right', va='top',
                     transform=state_panel.transAxes)
    state_panel.set_rasterization_zorder(0)

    # Plot IKr - Lei model
    lei_log = myokit.DataLog.load_csv(os.path.join(
        data_dir, f'{APmodel}_APclamp_Lei_current.csv'))
    ikr_axs[1].plot(lei_log.time(), lei_log['ikr.IKr'], '--',
                    color=model_color(m))

    # Compute area under the curve of both O and I states
    auc_lei = np.trapz(np.array(lei_log['ikr.O']) + np.array(lei_log['ikr.I']),
                       x=lei_log.time())
    print('Lei: ', auc_lei)

    # Plot state occupancy - Lei model
    state_panel = state_axs[1][m_row][m_col]
    state_panel.stackplot(
        lei_log.time(), *[lei_log['ikr.' + s] for s in lei_states],
        labels=[s for s in lei_states], colors=color_seq[:4], zorder=-10,
        edgecolor='k')
    state_panel.set_ylim(bottom=0, top=1)
    state_panel.text(0.95, 0.95, APmodel, fontsize=8, ha='right', va='top',
                     transform=state_panel.transAxes)
    state_panel.set_rasterization_zorder(0)

state_axs[0][0][0].legend(ncol=6, loc='lower left', handlelength=1,
                          bbox_to_anchor=(0, 1.0),
                          columnspacing=1, labelspacing=0.3)
state_axs[1][0][0].legend(ncol=4, loc='lower left', handlelength=1,
                          bbox_to_anchor=(0, 1.0),
                          columnspacing=1, labelspacing=0.3)
# Adjust axes
pulse_time = 1000
ikr_axs[1].set_xlim((0, 600))
ikr_axs[0].set_ylabel('Voltage \n (ms)')
ikr_axs[1].set_ylabel('Current \n (A/F)')
ikr_axs[0].sharex(ikr_axs[1])
ikr_axs[0].tick_params(labelbottom=False)
ikr_axs[1].set_xlabel('Time (ms)')

lines = []
for i in range(2):
    ikr_axs[i].set_box_aspect(0.15)
    ikr_axs[i].spines[['right', 'top']].set_visible(False)
    lines.append(matplotlib.lines.Line2D([0], [0], color=model_color(i), lw=5))
lines.append(matplotlib.lines.Line2D([0], [0], color='k', linestyle='-'))
lines.append(matplotlib.lines.Line2D([0], [0], color='k', linestyle='--'))

ikr_axs[0].legend(lines, APmodels + ['Li', 'Lei'], loc='upper left',
                  handlelength=1.5, bbox_to_anchor=(1.01, 0.5), ncol=1)

for i in range(2):
    fig.sharex(['Time (ms)'] * 2, [(0, pulse_time)] * 2,
               axs=state_axs[i], subgridspec=state_gridspec[i])
    fig.sharey(['State\noccupancy'],
               axs=state_axs[i], subgridspec=state_gridspec[i])

fig.fig.text(0.11, 0.925, '(A)', fontsize=9)
fig.fig.text(0.56, 0.925, '(B)', fontsize=9)
fig.fig.text(0.23, 0.66, '(C)', fontsize=9)
fig.fig.text(0.09, 0.35, '(D)', fontsize=9)
fig.fig.text(0.53, 0.35, '(E)', fontsize=9)
fig.savefig(os.path.join(fig_dir, "hERGmodels_APs_LiLei.svg"), format='svg')
