import myokit
import numpy as np
import os

import modelling

APmodel_dir = os.path.join(modelling.MAIN_DIR, 'math_model', 'AP_model')

AP_file_names = {
    'Grandi': {
        'AP': 'Grd-2010.mmt',
        'AP_IKr': 'Grd-2010-Li-SD.mmt',
        'label': 'Grandi (2010)'
    },
    'ORd-Li': {
        'AP': 'ORd-2011.mmt',
        'AP_IKr': 'ohara-cipa-2017.mmt',
        'label': "CiPAORdv1.0 (2017)"
    },
    'TTP': {
        'AP': 'TTP-2006.mmt',
        'AP_IKr': 'TTP-2006-Li-SD.mmt',
        'label': 'ten Tusscher (2006)'
    },
    'Tomek': {
        'AP': 'Tomek-2019.mmt',
        'AP_IKr': 'Tomek-2019-Li-SD.mmt',
        'label': 'Tomek (2019)'
    },
    'Tomek-Cl': {
        'AP': 'Tomek-Cl-2020.mmt',
        'AP_IKr': 'Tomek-Cl-2020-Li-SD.mmt',
        'label': 'Tomek (2019)'
    },
    'ORd-Lei': {
        'AP': 'ORd-2011.mmt',
        'AP_IKr': 'ORd-2011-Lei-SD.mmt',
        'label': 'Lei (2019)'
    },}


def APmodel_mmt(model, ikr_modified=True):
    """
    Take model name and returns mmt file directory
    """
    # If original model is requested
    if ikr_modified:
        file_key = 'AP_IKr'
    else:
        file_key = 'AP'
    return os.path.join(APmodel_dir, AP_file_names[model][file_key])


class ModelSimController(object):
    """
    To control the simulations of models
    """

    def __init__(self, APmodel_name, ikr_modified=True, period=1000):
        super(ModelSimController, self).__init__()

        # Load model
        self.model = myokit.load_model(APmodel_mmt(APmodel_name,
                                                   ikr_modified=ikr_modified))
        # parameters
        self._parameters = {}

        self.sim = myokit.Simulation(self.model)
        self.sim.set_tolerance(abs_tol=modelling.ABS_TOL,
                               rel_tol=modelling.REL_TOL)
        self.period = period
        protocol = myokit.pacing.blocktrain(period=self.period, duration=0.5,
                                            offset=50, level=1)
        self.sim.set_protocol(protocol)
        del(protocol)

        self.initial_state = self.sim.state()

        model_current_keys = modelling.model_naming.model_current_keys
        current_keys = model_current_keys[APmodel_name]
        self.ikr_key = current_keys['IKr']
        self.Vm_key = current_keys['Vm']
        self.ikr_key_head = self.ikr_key[:self.ikr_key.index('.')]
        self.ikr_component = self.model.get(self.ikr_key_head)

    def set_SD_parameters(self, drug, ikr_model='Li'):
        SD_params = modelling.BindingParameters().load_SD_parameters(
            drug, ikr_model=ikr_model)
        param_dict = SD_params.to_dict(orient='index')[drug]

        param_dict['Kt'] = 3.5e-5
        del param_dict['Cmax']

        dict_contents = list(param_dict.items())
        param_dict.clear()
        for k, content in dict_contents:
            param_dict[self.ikr_key_head + '.' + k] = content

        self._parameters.update(param_dict)

    def set_conductance(self, ikr_conductance):
        self._conductance = float(ikr_conductance)
        param_dict = {self.ikr_key_head + '.gKr': self._conductance}
        self._parameters.update(param_dict)

    def set_ikr_rescale(self, ikr_rescale):
        self._ikr_tune = float(ikr_rescale)
        param_dict = {'tune.ikr_rescale': self._ikr_tune}
        self._parameters.update(param_dict)

    def reset_parameters(self):
        param_dict = {
            self.ikr_key_head + "Vhalf":
                self.model.get(self.ikr_component.var('Vhalf')).eval(),
            self.ikr_key_head + "Kmax":
                self.model.get(self.ikr_component.var('Kmax')).eval(),
            self.ikr_key_head + "Ku":
                self.model.get(self.ikr_component.var('Ku')).eval(),
            self.ikr_key_head + "n":
                self.model.get(self.ikr_component.var('n')).eval(),
            self.ikr_key_head + "halfmax":
                self.model.get(self.ikr_component.var('halfmax')).eval(),
            self.ikr_key_head + "Kt":
                self.model.get(self.ikr_component.var('Kt')).eval(),
            self.ikr_key_head + "gKr":
                self.model.get(self.ikr_component.var('gKr')).eval(),
            'tune.ikr_rescale': 1}
        self._parameters.update(param_dict)
        del(param_dict)

    def get_parameters(self):
        return self._parameters

    def set_conc(self, concentration):
        if concentration < 0:
            ValueError('Drug concentration is lower than 0.')
        self._conc = float(concentration)
        param_dict = {self.ikr_key_head + '.D': self._conc}
        self._parameters.update(param_dict)

    def _set_parameters(self, parameters):
        """
        Set the parameter values
        """
        for p in parameters.keys():
            self.sim.set_constant(p, parameters[p])

    def set_protocol(self, times, voltages):
        self.sim.set_fixed_form_protocol(times, voltages)

    def update_initial_state(self, paces=1000):
        """
        Mainly to simulate till steady state for drug-free conditions
        """

        # Update fixed parameters
        self._set_parameters(self._parameters)

        # Ensure initial condition
        self.sim.set_time(0)
        self.sim.reset()
        self.sim.set_state(self.initial_state)

        self.sim.pre(paces * self.period)
        self.initial_state = self.sim.state()

    def simulate(self, prepace=1000, save_signal=1, timestep=0.1,
                 log_var=None, reset=True):

        self._set_parameters(self._parameters)

        self.sim.set_time(0)
        if reset:
            self.sim.reset()
            self.sim.set_state(self.initial_state)

        self.sim.pre(self.period * prepace)
        log = self.sim.run(self.period * save_signal, log=log_var,
                           log_interval=timestep)
        log = log.npview()
        if save_signal > 1:
            log = log.fold(self.period)

        return log

    def APD90(self, log, offset, protocol_duration):
        """
        Compute APD90
        """
        time_log = log.time()
        Vm_log = log[self.Vm_key]
        timestep = (max(time_log) - min(time_log)) / len(time_log)
        pulse_num = round(len(Vm_log) / len(time_log))

        APD90s = []
        for i in range(pulse_num):
            AP_Vm = Vm_log[
                round(i * protocol_duration / timestep):
                round((i + 1) * (protocol_duration + offset) / timestep)]
            APamp = max(AP_Vm) - min(AP_Vm)
            APD90_v = min(AP_Vm) + 0.1 * APamp

            time_log = np.concatenate((time_log,
                                       time_log + max(time_log) + timestep))

            min_APD = int(5 / timestep)
            offset_index = int(offset / timestep)

            start_time = time_log[offset_index]
            end_time = None
            for ind, v in enumerate(AP_Vm[offset_index + min_APD:]):
                if v < APD90_v:
                    t_prev = time_log[offset_index + min_APD:][ind - 1]
                    v_prev = AP_Vm[offset_index + min_APD:][ind - 1]
                    t_current = time_log[offset_index + min_APD:][ind]
                    end_time = t_prev + (APD90_v - v_prev) / (v - v_prev) * \
                        (t_current - t_prev)
                    break

            if end_time is not None:
                APD90_value = end_time - start_time
            else:
                APD90_value = float('Nan')

            APD90s.append(APD90_value)

        return APD90s

    def extract_peak(self, signal_log, current_name=None):
        peaks = []
        if current_name is None:
            current_name = self.ikr_key
        pulses = len(signal_log.keys_like(current_name))

        if pulses == 0:
            peaks.append(np.max(signal_log[current_name]))

        for i in range(pulses):
            peaks.append(np.max(signal_log[current_name, i]))

        if peaks[0] == 0:
            peak_reduction = 0
        else:
            peak_reduction = (peaks[0] - peaks[-1]) / peaks[0]

        return peaks, peak_reduction
