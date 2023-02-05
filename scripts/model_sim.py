import myokit

import modelling

# Define protocol
pulse_time = 1000
protocol = myokit.pacing.blocktrain(pulse_time, 0.5, offset=10)

#
# Grandi (2010) model
#
# # Set up AP model
# APmodel = '../math_model/AP_model/Grd-2010.mmt'
# APmodel, _, x = myokit.load(APmodel)
# AP_model = modelling.Simulation(APmodel)
# AP_model.protocol = protocol

abs_tol = 1e-7
rel_tol = 1e-8
# log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)

# # Plot output
# fig = modelling.figures.FigureStructure(figsize=(5, 3), gridspec=(2, 1))
# plot = modelling.figures.FigurePlot()

# plot.add_single(fig.axs[0][0], log, 'interface.I_app')
# plot.add_single(fig.axs[1][0], log, 'membrane_potential.V_m')
# # plot.state_occupancy_plot(fig.axs[2][0], log, AP_model)

# fig.sharex(['Time (ms)'], [(0, pulse_time)])
# fig.sharey(['Current\n(A/F)', 'Voltage\n(mV)'])

# fig_dir = '../figures/basic_sim/'
# fig.savefig(fig_dir + "grd_sim.pdf")

# # Plot output
# fig = modelling.figures.FigureStructure(figsize=(8, 5), gridspec=(6, 2),
#                                         height_ratios=[1] * 6, hspace=0.2)

# current_list = ['membrane_potential.V_m', 'Ca_Concentrations.Ca_i',
#                 'Ca_Concentrations.Ca_j', 'Ca_Concentrations.Ca_sl',
#                 'I_Na.I_Na', 'I_Ca.I_Catot', 'I_Kr.I_kr', 'I_Ks.I_ks',
#                 'I_to.I_to', 'I_Ki.I_ki', 'I_NCX.I_ncx', 'I_NaK.I_nak']
# current_name = ['Vm', '[Ca]i', '[Ca]jct', '[Ca]SL',
#                 'I_Na', 'I_CaL', 'I_Kr', 'I_Ks',
#                 'I_to', 'I_K1', 'I_NCX', 'I_NaK']

# for i in range(len(current_list)):
#     row = i % 6
#     col = i // 6
#     plot.add_single(fig.axs[row][col], log, current_list[i])
#     fig.axs[row][col].spines[['right', 'top']].set_visible(False)
#     fig.axs[row][col].set_ylabel(current_name[i])
# # plot.state_occupancy_plot(fig.axs[2][0], log, AP_model)

# fig.sharex(['Time (ms)'] * 2, [(0, 500)] * 2)
# # fig.sharey(['Current\n(A/F)', 'Voltage\n(mV)'])

# fig_dir = '../figures/basic_sim/'
# fig.savefig(fig_dir + "grd_sim_currents.pdf")

#
# ten Tusscher (2006) model
#
# Set up AP model
APmodel = '../math_model/AP_model/TTP-2006.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel)
AP_model.protocol = protocol

log = AP_model.model_simulation(1000, abs_tol=abs_tol, rel_tol=rel_tol)

print(min(log['membrane.V'][500:]), max(log['membrane.V']))
# Plot TTP's action potential
fig = modelling.figures.FigureStructure(figsize=(5, 3), gridspec=(2, 1))
plot = modelling.figures.FigurePlot()

plot.add_single(fig.axs[0][0], log, 'membrane.i_Stim')
plot.add_single(fig.axs[1][0], log, 'membrane.V')
# plot.state_occupancy_plot(fig.axs[2][0], log, AP_model)

fig.sharex(['Time (ms)'], [(0, pulse_time)])
fig.sharey(['Current\n(A/F)', 'Voltage\n(mV)'])

fig_dir = '../figures/basic_sim/'
fig.savefig(fig_dir + "ttp_sim.pdf")