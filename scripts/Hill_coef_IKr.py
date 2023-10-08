## Get the Hill coefficient
## Save it in csv format - based on drug or model or per model per drug

# Compare the AP, APD90 and qNet between the ORd-Lei-SD models and the ORd-Lei-CS models
# for a given drug at a chosen IKr tuning method
# Output:
# AP
# 1. 2 pulses of action potential simulated from both models (steady state).
# 2. APD90 of both pulses for both models (steady state).
# APD
# 1. APD90 for both models at various drug concentration.
# qNet
# 1. qNet for both models at various drug concentration.

import argparse
import glob
import matplotlib
import myokit
import numpy as np
import os
import pandas as pd

import modelling

# Define AP model, drug and tuning method
parser = argparse.ArgumentParser(
    description="Comparison between ORd-Lei-SD model and the ORd-Lei-CS model"
    " for a given drug")
parser.add_argument("drug", help="Drug")
parser.add_argument('--plot', action='store_true',
                    help="Plot comparison of IKr")
parser.add_argument("--cache", action='store_true',
                    help='Use cache of Hill curve')
args = parser.parse_args()

drug = args.drug

# Define Hill curve directory
data_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                        'Lei', drug)
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
Hill_curve_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                              'Hill_curves', drug, 'Lei_Hill.csv')

# Set up IKr model, protocol period and drug concentration
IKr_sim = modelling.ModelSimController('Lei')
pulse_time = 25e3
drug_conc = modelling.SD_details.drug_concentrations[drug]['coarse']
Hill_model = modelling.HillModel()

if not args.cache:

    # Define drug and protocol
    Milnes_protocol = modelling.ProtocolLibrary().Milnes(pulse_time)

    # Load IKr model
    IKr_sim.sim.set_protocol(Milnes_protocol)
    IKr_sim._cycle_length = pulse_time

    # Simulate IKr of the SD model for a range of drug concentrations
    # Extract the peak of IKr
    peaks = []
    IKr_sim.set_SD_parameters(drug, ikr_model='Lei')
    for i in range(len(drug_conc)):
        IKr_sim.set_conc(drug_conc[i])
        log = IKr_sim.simulate()
        peak = IKr_sim.extract_peak(log)
        peaks.append(peak[-1])

        log.save_csv(os.path.join(data_dir,
                                  'SD_current_' + str(drug_conc[i]) + '.csv'))

    # Normalise drug response (peak current)
    peaks = (peaks - min(peaks)) / (max(peaks) - min(peaks))

    # Fit drug response to Hill curve
    optimiser = modelling.HillModelOpt(Hill_model)
    estimates, _ = optimiser.optimise(drug_conc, peaks)
    Hill_df = pd.DataFrame(estimates, index=['Hill coefficient', 'IC50'])
    Hill_df.to_csv(Hill_curve_dir)

    # Compare peak current
    IKr_sim.reset_parameters()
    for i in range(len(drug_conc)):
        reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
        IKr_sim.set_conc(0)
        IKr_sim.set_CS_parameter(reduction_scale)
        log = IKr_sim.simulate()
        IKr_sim.set_CS_parameter(1)

        log.save_csv(os.path.join(data_dir,
                                  'CS_current_' + str(drug_conc[i]) + '.csv'))
else:
    Hill_coef_df = pd.read_csv(Hill_curve_dir, index_col=[0],
                               skipinitialspace=True)
    estimates = Hill_coef_df.T.to_numpy()[0]

if args.plot:
    # Plot results
    fig = modelling.figures.FigureStructure(figsize=(8, 2), gridspec=(1, 2),
                                            width_ratios=[2.5, 1.2],
                                            plot_in_subgrid=True)
    plot = modelling.figures.FigurePlot()
    cmap = matplotlib.cm.get_cmap('viridis')

    subgridspecs = [(1, 2), (1, 1)]
    subgs = []
    for i in range(2):
        subgs.append(fig.gs[i].subgridspec(*subgridspecs[i], wspace=0.08,
                                           hspace=0.08))
    axs = [[[fig.fig.add_subplot(subgs[k][i, j]) for j in range(
            subgridspecs[k][1])] for i in range(subgridspecs[k][0])]
           for k in range(len(subgs))]

    # Read files name of IKr data
    SD_data_files = glob.glob(os.path.join(data_dir, 'SD_current_*.csv'))
    CS_data_files = glob.glob(os.path.join(data_dir, 'CS_current_*.csv'))

    SD_data_files = sorted(SD_data_files,
                           key=lambda x: float(os.path.basename(x)[11:-4]))
    CS_data_files = sorted(CS_data_files,
                           key=lambda x: float(os.path.basename(x)[11:-4]))

    conc_label = [os.path.basename(fpath)[11:-4] for fpath in SD_data_files]
    conc_label = [r"$10^{:d}$ nM".format(int(np.log10(float(i))))
                  if float(i) >= 1e3 else i + 'nM' for i in conc_label]

    # Load IKr data
    SD_hERG_log = []
    CS_hERG_log = []
    for i in range(len(SD_data_files)):
        SD_hERG_log.append(myokit.DataLog.load_csv(SD_data_files[i]))
        CS_hERG_log.append(myokit.DataLog.load_csv(CS_data_files[i]))

    # Plot Ikr
    panel1 = axs[0]
    plot.add_multiple(panel1[0][0], SD_hERG_log, IKr_sim.ikr_key,
                      labels=conc_label, color=cmap)
    plot.add_multiple(panel1[0][1], CS_hERG_log, IKr_sim.ikr_key,
                      labels=conc_label, color=cmap)

    # Adjust figure details
    panel1[0][1].legend(handlelength=0.9, ncol=2, columnspacing=0.9)
    panel1[0][0].set_title('SD model')
    panel1[0][1].set_title('CS model')
    fig.sharex(['Time (s)'] * 2, [(0, pulse_time)] * 2,
               axs=panel1, subgridspec=subgridspecs[0])
    fig.sharey(['Current (A/F)'],
               axs=panel1, subgridspec=subgridspecs[0])
    panel1[0][0].spines[['right', 'top']].set_visible(False)
    panel1[0][1].spines[['right', 'top']].set_visible(False)
    fig.adjust_ticks(panel1[0][0], pulse_time)
    fig.adjust_ticks(panel1[0][1], pulse_time)

    # Plot fitting result of Hill curve
    max_grid = np.ceil(np.log(drug_conc[-1]))
    conc_grid = np.arange(-3, max_grid + 1, 0.5)

    panel2 = axs[1]
    if not args.cache:
        panel2[0][0].plot(np.log(drug_conc[1:]), peaks[1:], 'o',
                          label='peak current')
    panel2[0][0].plot(conc_grid,
                      Hill_model.simulate(estimates[:2], np.exp(conc_grid)),
                      'k', label='fitted Hill eq')
    panel2[0][0].set_xlabel('Drug concentration (log)')
    panel2[0][0].set_ylabel('Normalised peak current')

    fig_dir = os.path.join(modelling.FIG_DIR, 'kinetics_comparison',
                           'Lei', 'SD2CS_Hill')
    if not os.path.isdir(fig_dir):
        os.makedirs(fig_dir)
    fig.savefig(os.path.join(fig_dir, drug + '.pdf'))

# ################################
# # Propagate to action potential
# ################################
# model_keys = modelling.model_naming.model_current_keys["ORd-Lei"]
# current_key = model_keys['IKr']
# APsim = modelling.ModelSimController('ORd-Lei')

# scale_df = pd.read_csv(os.path.join(modelling.RESULT_DIR,
#                        'ORd-Lei_conductance_scale.csv'),
#                        index_col=[0], skipinitialspace=True)
# conductance_scale = scale_df.loc[args.ikr_tuning].values[0]
# APsim.set_ikr_rescale(conductance_scale)

# if MODE in ["AP", "APD", "APD-qNet"]:
#     save_signal = 2
#     apd_CS = []
#     apd_SD = []

#     # Simulate AP of the AP-SD model and the AP-CS model
#     # Compute APD90
#     APsim.set_SD_parameters(drug, ikr_model='Lei')
#     for i in range(len(drug_conc)):
#         print('simulating concentration: ' + str(drug_conc[i]))

#         APsim.set_conc(drug_conc[i])
#         log = APsim.simulate(save_signal=save_signal)
#         if MODE == "AP":
#             log.save_csv(os.path.join(data_dir,
#                                       'SD_AP_' + str(drug_conc[i]) + '.csv'))

#         apd90 = APsim.APD90(log)
#         apd_SD.append(apd90)

#         reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
#         APsim.set_conc(0)
#         APsim.set_CS_parameter(reduction_scale)
#         log = APsim.simulate(save_signal=save_signal)
#         APsim.set_CS_parameter(1)

#         log.save_csv(os.path.join(data_dir, 'CS_AP_' + str(drug_conc[i]) + '.csv'))

#         apd90 = APsim.APD90(log)
#         apd_CS.append(apd90)

#         print('done concentration: ' + str(drug_conc[i]))

# # Save simulated APD90 of both the AP-SD model and the AP-CS model
# column_name = ['pulse ' + str(i) for i in range(save_signal)]
# APD_trapping_df = pd.DataFrame(apd_SD, columns=column_name)
# APD_trapping_df['drug concentration'] = drug_conc
# APD_trapping_df.to_csv(os.path.join(data_dir, 'SD_APD_pulses' +
#                        str(int(save_signal)) + '.csv'))
# APD_conductance_df = pd.DataFrame(apd_CS, columns=column_name)
# APD_conductance_df['drug concentration'] = drug_conc
# APD_conductance_df.to_csv(os.path.join(data_dir, 'CS_APD_pulses' +
#                           str(int(save_signal)) + '.csv'))

# # Define drug concentration range for steady state APD90 comparison between
# # models
# drug_conc = modelling.SD_details.drug_concentrations[drug]['fine']

# APD_conductance = []
# APD_trapping = []

# for i in range(len(drug_conc)):
#     print('simulating concentration: ' + str(drug_conc[i]))
#     APsim.set_conc(drug_conc[i])
#     log = APsim.simulate(prepace=999, save_signal=save_signal)

#     apd90 = APsim.APD90(log)
#     APD_trapping.append(apd90)

#     reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
#     APsim.set_conc(0)
#     APsim.set_CS_parameter(reduction_scale)
#     # save signal --> output pulse num
#     log = APsim.simulate(save_signal=save_signal)
#     APsim.set_CS_parameter(1)

#     apd90 = APsim.APD90(log)
#     APD_conductance.append(apd90)

#     print('done concentration: ' + str(drug_conc[i]))

# # Compute APD90 with AP behaviour in alternating cycles
# APD_trapping = [float('nan') if np.isnan(i).any() else max(i)
#                 for i in APD_trapping]
# APD_conductance = [float('nan') if np.isnan(i).any() else max(i)
#                    for i in APD_conductance]

# # Compute qNet
# pulse_time = 2000
# base_log_key = [APsim.time_key, APsim.Vm_key]
# qNet_current_list = [model_keys['ICaL'], model_keys['INaL'], current_key,
#                      model_keys['IKs'], model_keys['IK1'], model_keys['Ito']]
# qNet_current_list = [i for i in qNet_current_list if i is not None]

# APsim.sim.set_protocol(myokit.pacing.blocktrain(period=pulse_time,
#                                                 duration=0.5,
#                                                 offset=50))
# APsim._cycle_length = pulse_time
# APsim.reset_parameters()
# APsim.update_initial_state()

# qNet_SD_arr = []
# qNet_CS_arr = []
# print('Computing qNet...')
# APsim.set_SD_parameters(drug)
# for i in range(len(drug_conc)):
#     print('simulating concentration: ' + str(drug_conc[i]))

#     APsim.set_conc(drug_conc[i])
#     log = APsim.simulate(timestep=0.01,
#                          log_var=base_log_key + qNet_current_list)

#     qNet = APsim.qNet(log)
#     qNet_SD_arr.append(qNet)

#     # Run simulation for the ORd-CS model till steady state
#     reduction_scale = Hill_model.simulate(estimates[:2], drug_conc[i])
#     APsim.set_conc(0)
#     APsim.set_CS_parameter(reduction_scale)
#     log = APsim.simulate(timestep=0.01,
#                          log_var=base_log_key + qNet_current_list)
#     APsim.set_CS_parameter(1)

#     qNet = APsim.qNet(log)
#     qNet_CS_arr.append(qNet)

#     print('done concentration: ' + str(drug_conc[i]))

# # Save APD90 data
# APD_trapping_df = pd.DataFrame(np.array(APD_trapping), columns=['APD'])
# APD_trapping_df['drug concentration'] = drug_conc
# APD_trapping_df['qNet'] = qNet_SD_arr
# APD_trapping_df.to_csv(os.path.join(data_dir, 'SD_APD_fine.csv'))
# APD_conductance_df = pd.DataFrame(np.array(APD_conductance), columns=['APD'])
# APD_conductance_df['drug concentration'] = drug_conc
# APD_conductance_df['qNet'] = qNet_CS_arr
# APD_conductance_df.to_csv(os.path.join(data_dir, 'CS_APD_fine.csv'))
