import numpy as np
import pandas as pd

import modelling


class ModelComparison(object):
    """
    To create a class to run the model comparison.
    """

    def __init__(self, APsimulator, IKr_simulator=None):
        super(ModelComparison, self).__init__()

        # check if the simulator is from the ModelSimController class
        if not isinstance(APsimulator, modelling.ModelSimController):
            raise TypeError(
                'The model has to be an instance of a '
                'modelling.ModelSimController.')
        self.APsim = APsimulator

        if IKr_simulator is not None and not \
                isinstance(IKr_simulator, modelling.ModelSimController):
            raise TypeError(
                'The model has to be an instance of a '
                'modelling.ModelSimController.')
        self.IKr_sim = IKr_simulator

        # Set normalising constant for drug input
        self.drug_norm_constant = 1
        self.Hill_model = modelling.HillModel()

    def prepare_inputs(self, input, skip_ikr=False):
        self.id_index = input['id']
        param_values = input['param_values']
        self._set_parameters(param_values)

        if skip_ikr:
            drug_conc = input['drug_conc_Hill'].values[0]
            self.drug_conc_ikr = drug_conc[~np.isnan(drug_conc)]
            self.Hill_coef = input['Hill_curve'].values[0]

    def _set_parameters(self, param_df):
        if 'Cmax' in param_df.columns:
            param_df = param_df.drop(columns=['Cmax'])
        self.param_df = param_df

    def normalise_drug_conc(self):
        '''
        Set normalising constant for drug input
        '''
        if not hasattr(self, 'param_df'):
            raise ValueError('Please set parameter values with '
                             'method `set_parameters()`')

        self.orig_halfmax = self.param_df['halfmax'][0]
        self.param_df['halfmax'] = 1

        # Calculate the normalising constant
        Hill_n = self.param_df['n'][0]
        self.drug_norm_constant = np.power(self.orig_halfmax, 1 / Hill_n)

    def _get_peak(self, drug_conc):

        self.IKr_sim.set_conc(drug_conc)
        log = self.IKr_sim.simulate(
            log_var=[self.IKr_sim.time_key, self.IKr_sim.ikr_key])
        peak = self.IKr_sim.extract_peak(log)

        return peak[-1]

    def get_drug_effect(self, drug_conc=None,
                        upper_thres=0.9, lower_thres=0.05,
                        max_counter=20, parallel=True):

        self.IKr_sim.reset_parameters()
        self.IKr_sim.set_parameters(self.param_df)

        if drug_conc is None:
            drug_conc = list(np.append(0, 10**np.linspace(-1, 1, 5)))

        drug_conc = [i / self.drug_norm_constant for i in drug_conc]

        peaks = []
        for i in range(len(drug_conc)):
            peak = self._get_peak(drug_conc[i])
            peaks.append(peak)

        peaks_norm = (peaks - min(peaks)) / (max(peaks) - min(peaks))

        # Make sure there are enough data points for the head of Hill curve
        data_pt_checker = [True if i > upper_thres else False
                           for i in peaks_norm]
        counter = 0
        while sum(data_pt_checker) < 3 and counter < max_counter:
            drug_conc.insert(1, drug_conc[1] / np.sqrt(10))
            peak = self._get_peak(drug_conc[1])
            peaks.insert(1, peak)

            peaks_norm = (peaks - min(peaks)) / (max(peaks) - min(peaks))
            data_pt_checker = [True if i > upper_thres else False
                               for i in peaks_norm]
            counter += 1

        if counter == max_counter:
            return 'Hill curve did not form.', drug_conc, peaks_norm

        # Make sure there are enough data points for the tail of Hill curve
        data_pt_checker = [True if i < lower_thres else False
                           for i in peaks_norm]
        counter = 0
        while sum(data_pt_checker) < 3 and counter < max_counter:
            drug_conc = drug_conc + [max(drug_conc) * np.sqrt(10)]
            peak = self._get_peak(drug_conc[-1])
            peaks.append(peak)

            peaks_norm = (peaks - min(peaks)) / (max(peaks) - min(peaks))
            data_pt_checker = [True if i < lower_thres else False
                               for i in peaks_norm]
            counter += 1

        if counter == max_counter:
            return 'Hill curve did not form.', drug_conc, peaks_norm

        # Fit Hill curve
        optimiser = modelling.HillModelOpt(self.Hill_model)
        Hill_curve, _ = optimiser.optimise(drug_conc, peaks_norm,
                                           parallel=parallel)

        self.drug_conc_ikr = drug_conc
        self.peaks_norm = peaks_norm
        self.Hill_coef = Hill_curve[:2]

        del drug_conc
        del peaks_norm

    def _get_apd90(self, drug_conc, save_signal=2):

        self.APsim.set_conc(drug_conc)
        log = self.APsim.simulate(
            save_signal=save_signal,
            log_var=[self.APsim.time_key, self.APsim.Vm_key])

        apd90_SD = self.APsim.APD90(log)

        # Run simulation for conductance model
        reduction_scale = self.Hill_model.simulate(
            self.Hill_coef, drug_conc)
        self.APsim.set_conc(0)
        self.APsim.set_CS_parameter(reduction_scale)
        log = self.APsim.simulate(
            save_signal=save_signal,
            log_var=[self.APsim.time_key, self.APsim.Vm_key])
        self.APsim.set_CS_parameter(1)

        apd90_CS = self.APsim.APD90(log)

        del log

        return apd90_SD, apd90_CS

    def get_APD(self, drug_conc=None, save_signal=2, data_points=20, EAD=True):

        self.APsim.reset_parameters()
        self.APsim.set_parameters(self.param_df)

        APD_trapping = []
        APD_conductance = []
        if drug_conc is None:
            if hasattr(self, 'drug_conc_ikr'):
                drug_conc = 10**np.linspace(np.log10(self.drug_conc_ikr[1]),
                                            np.log10(max(self.drug_conc_ikr)),
                                            data_points)
            else:
                drug_conc = 10**np.linspace(-1, 5, data_points)
        drug_conc = list(drug_conc)

        for i in range(len(drug_conc)):
            apd90_SD, apd90_CS = self._get_apd90(drug_conc[i],
                                                 save_signal=save_signal)

            APD_trapping.append(apd90_SD)
            APD_conductance.append(apd90_CS)

        APD_trapping = [float('nan') if np.isnan(i).any() else max(i)
                        for i in APD_trapping]
        APD_conductance = [float('nan') if np.isnan(i).any() else max(i)
                           for i in APD_conductance]
        if EAD:
            checker_trapping = [np.isnan(i) for i in APD_trapping]
            checker_conductance = [np.isnan(i) for i in APD_conductance]
            checker_count = min(sum(checker_trapping),
                                sum(checker_conductance))
            counter = 0
            while checker_count < 2 and counter < 10:
                drug_conc = drug_conc + [max(drug_conc) * np.sqrt(10)]
                apd90_SD, apd90_CS = self._get_apd90(drug_conc[-1],
                                                     save_signal=save_signal)

                apd90_max = float('nan') if np.isnan(apd90_SD).any() \
                    else max(apd90_SD)
                APD_trapping.append(apd90_max)

                apd90_max = float('nan') if np.isnan(apd90_CS).any() \
                    else max(apd90_CS)
                APD_conductance.append(apd90_max)

                checker_trapping = [np.isnan(i) for i in APD_trapping]
                checker_conductance = [np.isnan(i) for i in APD_conductance]
                checker_count = min(sum(checker_trapping),
                                    sum(checker_conductance))

                counter += 1

        self.drug_conc_AP = drug_conc
        self.APD_trapping = APD_trapping
        self.APD_conductance = APD_conductance

        del APD_trapping
        del APD_conductance

    def RMSE(self):

        square_sum = 0
        count = 0
        for i in range(len(self.APD_trapping)):
            if not np.isnan(self.APD_trapping[i]) and not \
                    np.isnan(self.APD_conductance[i]):
                square_sum += (self.APD_trapping[i] -
                               self.APD_conductance[i])**2
                count += 1

        if count != 0:
            self.RMSError = square_sum / count
            self.RMSError = np.sqrt(self.RMSError)
        else:
            self.RMSError = float('nan')

    def ME(self):

        diff_sum = 0
        count = 0
        for i in range(len(self.APD_trapping)):
            if not np.isnan(self.APD_trapping[i]) and not \
                    np.isnan(self.APD_conductance[i]):
                diff_sum += (self.APD_trapping[i] - self.APD_conductance[i])
                count += 1

        if count != 0:
            self.MAError = diff_sum / count
        else:
            self.MAError = float('nan')

    def process_data(self, save_orig_halfmax=False):

        conc_ikr_ind = ['conc_' + str(i) for i, _ in
                        enumerate(self.drug_conc_ikr)]
        conc_AP_ind = ['conc_' + str(i) for i, _ in
                       enumerate(self.drug_conc_AP)]

        index_dict = {'drug_conc_ikr': conc_ikr_ind,
                      'Hill_curve': ['Hill_coef', 'IC50'],
                      'param_values': self.param_df.keys(),
                      'drug_conc_AP': conc_AP_ind,
                      'APD_trapping': conc_AP_ind,
                      'APD_conductance': conc_AP_ind,
                      'RMSE': ['RMSE'], 'ME': ['ME']}
        if hasattr(self, 'peaks_norm'):
            index_dict.update({'peak_current': conc_ikr_ind})

        all_index = [(i, j) for i in index_dict.keys() for j in index_dict[i]]
        df_index = pd.MultiIndex.from_tuples(all_index)

        if save_orig_halfmax:
            self.param_df['halfmax'] = self.orig_halfmax

        data_list = list(self.drug_conc_ikr) + list(self.Hill_coef) + \
            list(self.param_df.values[0]) + list(self.drug_conc_AP) + \
            self.APD_trapping + self.APD_conductance + \
            [self.RMSError] + [self.MAError]
        if hasattr(self, 'peaks_norm'):
            data_list += list(self.peaks_norm)
        df = pd.DataFrame(data_list, index=df_index, columns=[self.id_index])

        return df.T
