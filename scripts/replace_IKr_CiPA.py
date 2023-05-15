# Introduces the idea of trapping and justifies the use of the Milnes protocol
import matplotlib
import myokit
import myokit.lib.plots as mp

import modelling

fig_dir = '../figures/basic_sim/IKr_replacement/'

# Define protocol
cycle_length = 1000
protocol_offset = 50
protocol = myokit.pacing.blocktrain(cycle_length, 0.5, offset=protocol_offset)

# Define constants
repeats = 1000
abs_tol = 1e-7
rel_tol = 1e-8

#######################
#
# ORd model
#
#######################
# Load ORd model
APmodel_name = 'ORd'
APmodel = '../math_model/AP_model/ORd-2011.mmt'
# APmodel = '../math_model/AP_model/ohara-2011.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, base_constant=None)
AP_model.protocol = protocol

model_keys = modelling.ModelDetails().current_keys[APmodel_name]
time_key = model_keys['time']
Vm_key = model_keys['Vm']
current_key = model_keys['IKr']

# Simulate AP clamp protocol
APclamp = AP_model.model_simulation(repeats, abs_tol=abs_tol, rel_tol=rel_tol,
                                    timestep=None,
                                    log_var=[time_key, Vm_key])

times = APclamp[time_key]
clamp_protocol = APclamp[Vm_key]
t_max = times[-1]

hERGmodel = '../math_model/current_model/ohara-cipa-2017-IKr.mmt'
APsim = myokit.Simulation(APmodel)
APsim.set_fixed_form_protocol(times, clamp_protocol)
APsim.reset()
APsim.set_tolerance(abs_tol=abs_tol, rel_tol=rel_tol)

for pace in range(1000):
    APsim.run(t_max, log=myokit.LOG_NONE)
    APsim.set_time(0)
log = APsim.run(t_max)
log = log.npview()

plot = modelling.figures.FigurePlot()
fig = modelling.figures.FigureStructure(figsize=(6, 4), gridspec=(1, 2),
                                        wspace=0.1, plot_in_subgrid=True)
subgridspecs = [(2, 1)] * 2
subgs = []
for i in range(2):
    subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.1,
                                       hspace=0.1))
axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
    subgridspecs[k][1])] for i in range(subgridspecs[k][0])] for
    k in range(len(subgs))]

cycle_length_plot = 800
SD_panel = axs[0]
plot.add_single(SD_panel[0][0], log, Vm_key)
plot.add_single(SD_panel[1][0], log, current_key)
SD_panel[0][0].set_title("ORd (2011)")
AP_y_bottom1, AP_y_top1 = SD_panel[0][0].get_ylim()
IKr_y_bottom1, IKr_y_top1 = SD_panel[1][0].get_ylim()
fig.sharex(['Time (ms)'], [(0, cycle_length_plot)],
           axs=SD_panel, subgridspec=subgridspecs[0])

#######################
#
# AP-SD model
#
#######################
# Load AP-SD model
APmodel_name = 'AP-SD'
APmodel = '../math_model/AP_model/ohara-cipa-2017.mmt'
APmodel, _, x = myokit.load(APmodel)
AP_model = modelling.Simulation(APmodel, current_head_key='ikr')

model_keys = modelling.ModelDetails().current_keys[APmodel_name]
time_key = model_keys['time']
Vm_key = model_keys['Vm']
current_key = model_keys['IKr']

APsim = myokit.Simulation(APmodel)
APsim.set_fixed_form_protocol(times, clamp_protocol)
APsim.reset()
APsim.set_tolerance(abs_tol=abs_tol, rel_tol=rel_tol)

for pace in range(1000):
    APsim.run(t_max, log=myokit.LOG_NONE)
    APsim.set_time(0)
log = APsim.run(t_max)
log = log.npview()

SD_panel = axs[1]
plot.add_single(SD_panel[0][0], log, Vm_key)
plot.add_single(SD_panel[1][0], log, current_key)
SD_panel[0][0].set_title("OHara-CiPA (2017)")
AP_y_bottom2, AP_y_top2 = SD_panel[0][0].get_ylim()
IKr_y_bottom2, IKr_y_top2 = SD_panel[1][0].get_ylim()
fig.sharex(['Time (ms)'], [(0, cycle_length_plot)],
           axs=SD_panel, subgridspec=subgridspecs[1])

#######################
#
# Adjust figures
#
#######################
AP_y_min = min(AP_y_bottom1, AP_y_bottom2)
AP_y_max = max(AP_y_top1, AP_y_top2)
IKr_y_min = min(IKr_y_bottom1, IKr_y_bottom2)
IKr_y_max = max(IKr_y_top1, IKr_y_top2)

for i in range(6):
    axs[i][1][0].set_ylim(IKr_y_min, IKr_y_max)
    axs[i][0][0].set_ylim(AP_y_min, AP_y_max)
    if i == 0:
        axs[i][0][0].set_ylabel('AP')
        axs[i][1][0].set_ylabel(r"$I_\mathrm{Kr}$")
    else:
        axs[i][0][0].set_yticklabels([])
        axs[i][1][0].set_yticklabels([])

fig.savefig(fig_dir + 'hERG_compare.pdf')
