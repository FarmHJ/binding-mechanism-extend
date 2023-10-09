# Get the Hill coefficient
# Save it in csv format - based on drug or model or per model per drug

import argparse
import glob
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
Hill_dir = os.path.join(modelling.RESULT_DIR, 'kinetics_comparison',
                        'Hill_curves', drug)
if not os.path.isdir(Hill_dir):
    os.makedirs(Hill_dir)
Hill_fname = 'Lei_Hill.csv'

# Set up IKr model, protocol period and drug concentration
IKr_sim = modelling.ModelSimController('Lei')
pulse_time = 25e3
drug_conc = modelling.SD_details.drug_concentrations[drug]['coarse']
Hill_model = modelling.HillModel()

if not args.cache:

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
    Hill_df.to_csv(os.path.join(Hill_dir, Hill_fname))

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
    Hill_coef_df = pd.read_csv(os.path.join(Hill_dir, Hill_fname),
                               index_col=[0], skipinitialspace=True)
    estimates = Hill_coef_df.T.to_numpy()[0]

if args.plot:
    # Plot results
    fig = modelling.figures.FigureStructure(figsize=(8, 2), gridspec=(1, 2),
                                            width_ratios=[2.5, 1.2],
                                            plot_in_subgrid=True)
    plot = modelling.figures.FigurePlot()

    subgridspecs = [(1, 2), (1, 1)]
    fig.subgrid(*subgridspecs[i], wspace=0.08, hspace=0.08)

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
    panel1 = fig.axs[0]
    plot.add_multiple(panel1[0][0], SD_hERG_log, IKr_sim.ikr_key,
                      labels=conc_label)
    plot.add_multiple(panel1[0][1], CS_hERG_log, IKr_sim.ikr_key,
                      labels=conc_label)

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

    panel2 = fig.axs[1]
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
