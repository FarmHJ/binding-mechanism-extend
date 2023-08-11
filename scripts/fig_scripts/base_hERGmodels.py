import myokit

import modelling

fig = modelling.figures.FigureStructure(figsize=(12, 8), gridspec=(2, 2),
                                        height_ratios=[2, 3], hspace=0.3,
                                        wspace=0.3, plot_in_subgrid=True)
plot = modelling.figures.FigurePlot()

subgridspecs = [(2, 1), (2, 1)]
subgs = []
for i in range(2):
    height_ratio = [1] + [2] * (subgridspecs[i][0] - 1)
    subgs.append(fig.gs[i + 2].subgridspec(*subgridspecs[i], wspace=0.1,
                                       hspace=0.05,
                                       height_ratios=height_ratio))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(len(subgs))]

data_dir = '../../simulation_data/background/'
fig_dir = '../../figures/background/'

model_list = ['ORd-CiPA', 'Lei']
# Load steady state IKr data for drug free, addition of a dofetilide-like drug
# and a verapamil-like drug conditions (AP clamp protocol)
# APclamp = myokit.DataLog.load_csv(
#     data_dir + 'APclamp.csv')
cipa_log = myokit.DataLog.load_csv(
    data_dir + model_list[0] + '_APclamp_current.csv')
lei_log = myokit.DataLog.load_csv(
    data_dir + model_list[1] + '_APclamp_current.csv')
log_all = [cipa_log, lei_log]

# Load AP model
APmodel = '../../math_model/AP_model/ohara-cipa-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
pulse_time = 1000

# Plot state occupancies of the hERG channel
color_seq = ['#7e7e7e', '#986464', '#989864', '#986496', '#988364',
             '#64986a', '#74a9cf', '#045a8d', '#2b8cbe']

panel_cipa = axs[0]
panel_lei = axs[1]

cipa_states = ['IC1', 'IC2', 'C1', 'C2', 'O', 'IO']
lei_states = ['y1', 'y2', 'y3', 'y4']

# Plot state occupancies of the hERG channel
panel_cipa[1][0].stackplot(cipa_log.time(),
                           *[cipa_log['ikr.' + s] for s in cipa_states],
                           labels=['ikr.' + s for s in cipa_states],
                           colors=color_seq[:6], zorder=-10)
panel_lei[1][0].stackplot(lei_log.time(),
                          *[lei_log['ikr.' + s] for s in lei_states],
                          labels=['ikr.' + s for s in lei_states],
                          colors=color_seq[:4], zorder=-10)

# for i in range(len(log_all)):
#     for j in range(len(log_all)):
#         if i == j:
#             plot.add_single(panel_cipa[0][i], log_all[j], 'membrane.V')
#             plot.add_single(panel4[1][i], log_all[j], 'ikr.IKr')
#         else:
#             plot.add_single(panel4[1][i], log_all[j], 'ikr.IKr',
#                             color='grey', alpha=0.5)
#     panel4[1][i].text(980, 0.9, drug_label[i], fontsize=8,
#                       ha='right', va='top')
plot.add_single(panel_cipa[0][0], log_all[0], 'ikr.IKr')
plot.add_single(panel_lei[0][0], log_all[1], 'ikr.IKr')


# # Adjust axes
# for col in range(3):
#     for tick in panel4[2][col].get_xticklabels():
#         tick.set_ha('right')
# fig.sharex(['Time (ms)'] * (len(drugs) + 1),
#            [(0, pulse_time)] * (len(drugs) + 1),
#            axs=panel4, subgridspec=subgridspecs[2])
# fig.sharey(['Voltage\n(mV)', 'Current\n(A/F)', 'State\noccupancy'],
#            axs=panel4, subgridspec=subgridspecs[2])

fig.savefig(fig_dir + "hERGmodels.pdf")
