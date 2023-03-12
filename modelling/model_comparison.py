import numpy as np


class ModelComparison(object):
    """
    To create a class to run the model comparison.
    """

    def __init__(self, drug_param_values):
        super(ModelComparison, self).__init__()

        self.drug_param_values = drug_param_values
        param_names = ['Vhalf', 'Kmax', 'Ku', 'N', 'EC50']

        # Check if the parameters are correct
        if sorted(list(self.drug_param_values.columns)) != sorted(param_names):
            ValueError("Parameters are incorrect. The parameters are " +
                       "Vhalf, Kmax, Ku, N and EC50")

    def APD_sim(self, AP_model, Hill_curve_coefs, drug_conc=None,
                steady_state_pulse=1000, save_signal=2, offset=50,
                data_points=20, EAD=False, norm_constant=1,
                abs_tol=1e-7, rel_tol=1e-8):

        base_conductance = AP_model.original_constants['gKr']
        APD_trapping = []
        APD_conductance = []
        if drug_conc is None:
            drug_conc = 10**np.linspace(-1, 5, data_points)
        drug_conc = list(drug_conc)

        for i in range(len(drug_conc)):
            # Run simulation for trapping model
            log = AP_model.custom_simulation(
                self.drug_param_values, drug_conc[i], steady_state_pulse,
                timestep=0.1, save_signal=save_signal, abs_tol=abs_tol,
                rel_tol=rel_tol, log_var=['engine.time', 'membrane.V'])

            # Compute APD90
            APD_trapping_pulse = []
            for pulse in range(save_signal):
                apd90 = AP_model.APD90(log['membrane.V', pulse], offset, 0.1)
                APD_trapping_pulse.append(apd90)
            APD_trapping.append(APD_trapping_pulse)

            # Run simulation for conductance model
            reduction_scale = self.Hill_model.simulate(
                Hill_curve_coefs, drug_conc[i] * norm_constant)
            d2 = AP_model.conductance_simulation(
                base_conductance * reduction_scale, steady_state_pulse,
                timestep=0.1, save_signal=save_signal, abs_tol=abs_tol,
                rel_tol=rel_tol, log_var=['engine.time', 'membrane.V'])

            # Compute APD90
            APD_conductance_pulse = []
            for pulse in range(save_signal):
                apd90 = AP_model.APD90(d2['membrane.V', pulse], offset, 0.1)
                APD_conductance_pulse.append(apd90)
            APD_conductance.append(APD_conductance_pulse)

        APD_trapping = [max(i) for i in APD_trapping]
        APD_conductance = [max(i) for i in APD_conductance]
        if EAD:
            checker_trapping = [True if i >= 1000 else False
                                for i in APD_trapping]
            checker_conductance = [True if i >= 1000 else False
                                   for i in APD_conductance]
            checker_count = min(sum(checker_trapping),
                                sum(checker_conductance))
            counter = 0
            while checker_count < 2 and counter < 10:
                drug_conc = drug_conc + [max(drug_conc) * np.sqrt(10)]
                log = AP_model.custom_simulation(
                    self.drug_param_values, drug_conc[-1], steady_state_pulse,
                    timestep=0.1, save_signal=save_signal, abs_tol=abs_tol,
                    rel_tol=rel_tol, log_var=['engine.time', 'membrane.V'])

                # Compute APD90
                APD_trapping_pulse = []
                for pulse in range(save_signal):
                    apd90 = AP_model.APD90(log['membrane.V', pulse],
                                           offset, 0.1)
                    APD_trapping_pulse.append(apd90)
                APD_trapping.append(max(APD_trapping_pulse))

                # Run simulation for conductance model
                reduction_scale = self.Hill_model.simulate(
                    Hill_curve_coefs, drug_conc[-1])
                d2 = AP_model.conductance_simulation(
                    base_conductance * reduction_scale, steady_state_pulse,
                    timestep=0.1, save_signal=save_signal, abs_tol=abs_tol,
                    rel_tol=rel_tol, log_var=['engine.time', 'membrane.V'])

                # Compute APD90
                APD_conductance_pulse = []
                for pulse in range(save_signal):
                    apd90 = AP_model.APD90(d2['membrane.V', pulse],
                                           offset, 0.1)
                    APD_conductance_pulse.append(apd90)
                APD_conductance.append(max(APD_conductance_pulse))

                checker_trapping = [True if i >= 1000 else False
                                    for i in APD_trapping]
                checker_conductance = [True if i >= 1000 else False
                                       for i in APD_conductance]
                checker_count = min(sum(checker_trapping),
                                    sum(checker_conductance))

                counter += 1

        return APD_trapping, APD_conductance, drug_conc

    def compute_RMSE(self, APD_trapping, APD_conductance):

        square_sum = 0
        count = 0
        for i in range(len(APD_trapping)):
            if APD_trapping[i] < 1000 and APD_conductance[i] < 1000:
                square_sum += (APD_trapping[i] - APD_conductance[i])**2
                count += 1

        if count != 0:
            RMSError = square_sum / count
            RMSError = np.sqrt(RMSError)
        else:
            RMSError = float('nan')

        return RMSError

    def compute_ME(self, APD_trapping, APD_conductance):

        diff_sum = 0
        count = 0
        for i in range(len(APD_trapping)):
            if APD_trapping[i] < 1000 and APD_conductance[i] < 1000:
                diff_sum += (APD_trapping[i] - APD_conductance[i])
                count += 1

        if count != 0:
            MAError = diff_sum / count
        else:
            MAError = float('nan')

        return MAError
