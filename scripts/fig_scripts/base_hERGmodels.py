import myokit
import os

import modelling

fig = modelling.figures.FigureStructure(figsize=(12, 6), gridspec=(2, 2),
                                        height_ratios=[1, 2], hspace=0.3,
                                        wspace=0.3, plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(3, 1), (3, 1)]
subgs = []
for i in range(2):
    height_ratio = [1] * 2 + [2]
    subgs.append(fig.gs[i + 2].subgridspec(*subgridspecs[i], wspace=0.1,
                                           hspace=0.05,
                                           height_ratios=height_ratio))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(len(subgs))]

data_dir = os.path.join(modelling.RESULT_DIR, 'background')
fig_dir = os.path.join(modelling.FIG_DIR, 'background')

model_list = ['Li', 'Lei']
# Load steady state IKr data for drug free, addition of a dofetilide-like drug
# and a verapamil-like drug conditions (AP clamp protocol)
cipa_log = myokit.DataLog.load_csv(
    os.path.join(data_dir, model_list[0] + '_APclamp_current.csv'))
lei_log = myokit.DataLog.load_csv(
    os.path.join(data_dir, model_list[1] + '_APclamp_current.csv'))
log_all = [cipa_log, lei_log]
pulse_time = 1000

# Plot state occupancies of the hERG channel
color_seq = ['#ff0000', '#ff5555', '#999999', '#666666', '#333333',
             '#808080', '#4d4d4d']

# 80% gray: 333333
# 70% gray: 4d4d4d
# 60% gray: 666666
# 50% gray: 808080
# 40% gray: 999999

panel_cipa = axs[0]
panel_lei = axs[1]

cipa_states = ['O', 'IO', 'IC1', 'C1', 'IC2', 'C2']
# state 'y3' is open state
# lei_states = ['y3', 'y2', 'y1', 'y4']
lei_states = ['O', 'I', 'CI', 'C']

# Plot state occupancies of the hERG channel
panel_cipa[2][0].stackplot(cipa_log.time(),
                           *[cipa_log['ikr.' + s] for s in cipa_states],
                           labels=[s for s in cipa_states],
                           colors=color_seq[:6], zorder=-10,
                           edgecolor='k')
panel_cipa[2][0].legend(ncol=3, loc='upper right', handlelength=1,
                        columnspacing=1, labelspacing=0.3)
panel_cipa[2][0].set_ylim(bottom=0, top=1)
panel_cipa[2][0].set_rasterization_zorder(0)
panel_lei[2][0].stackplot(lei_log.time(),
                          *[lei_log['ikr.' + s] for s in lei_states],
                          labels=[s for s in lei_states],
                          colors=color_seq[:4], zorder=-10,
                          edgecolor='k')
panel_lei[2][0].legend(ncol=2, loc='upper right', handlelength=1,
                       columnspacing=1, labelspacing=0.3)
panel_lei[2][0].set_ylim(bottom=0, top=1)
panel_lei[2][0].set_rasterization_zorder(0)

plot.add_single(panel_cipa[0][0], log_all[0], 'membrane.V')
plot.add_single(panel_lei[0][0], log_all[1], 'membrane.V')
plot.add_single(panel_cipa[1][0], log_all[0], 'ikr.IKr')
plot.add_single(panel_lei[1][0], log_all[1], 'ikr.IKr')

# Adjust axes
fig.sharex(['Time (ms)'], [(0, pulse_time)],
           axs=panel_cipa, subgridspec=subgridspecs[0])
fig.sharey(['Voltage\n(mV)', 'Current\n(A/F)', 'State\noccupancy'],
           axs=panel_cipa, subgridspec=subgridspecs[0])
fig.sharex(['Time (ms)'], [(0, pulse_time)],
           axs=panel_lei, subgridspec=subgridspecs[1])
fig.sharey(['Voltage\n(mV)', 'Current\n(A/F)', 'State\noccupancy'],
           axs=panel_lei, subgridspec=subgridspecs[1])

fig.fig.text(0.0855, 0.925, '(A)', fontsize=11)
fig.fig.text(0.52, 0.925, '(B)', fontsize=11)
fig.fig.text(0.0855, 0.56, '(C)', fontsize=11)
fig.fig.text(0.52, 0.56, '(D)', fontsize=11)
fig.savefig(os.path.join(fig_dir, "hERGmodels.svg"), format='svg')
