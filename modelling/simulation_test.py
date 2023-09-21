import myokit
import numpy as np
import os

import modelling

APmodel_dir = os.path.join(modelling.MAIN_DIR, 'math_model', 'AP_model')


def APmodel_mmt(model, ikr_modified=True):
    """
    Take model name and returns mmt file directory
    """
    AP_filenames = modelling.model_naming.AP_file_names
    # If original model is requested
    if ikr_modified:
        file_key = 'AP_IKr'
    else:
        file_key = 'AP'
    return os.path.join(APmodel_dir, AP_filenames[model][file_key])


class ModelSimController(object):
    """
    To control the simulations of models
    """

    def __init__(self, APmodel_name, ikr_modified=True, cycle_length=1000,
                 protocol_offset=50):
        super(ModelSimController, self).__init__()

        # Load model
        self.model = myokit.load_model(APmodel_mmt(APmodel_name,
                                                   ikr_modified=ikr_modified))
        # parameters
        self._parameters = {}

        self.sim = myokit.Simulation(self.model)
        self.sim.set_tolerance(abs_tol=modelling.ABS_TOL,
                               rel_tol=modelling.REL_TOL)
        self._cycle_length = cycle_length
        self.protocol_offset = protocol_offset
        protocol = myokit.pacing.blocktrain(period=self._cycle_length,
                                            duration=0.5,
                                            offset=self.protocol_offset,
                                            level=1)
        self.sim.set_protocol(protocol)
        del(protocol)

        self.initial_state = self.sim.state()

        model_current_keys = modelling.model_naming.model_current_keys
        self.current_keys = model_current_keys[APmodel_name]
        self.ikr_key = self.current_keys['IKr']
        self.Vm_key = self.current_keys['Vm']
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

    def set_CS_parameter(self, ikr_conductance_scale):
        self._conductance_scale = float(ikr_conductance_scale)
        param_dict = {self.ikr_key_head + '.gKr':
                      self.model.get(self.ikr_component.var('gKr')).eval()
                      * self._conductance_scale}
        self._parameters.update(param_dict)

    def set_ikr_rescale(self, ikr_rescale):
        self._ikr_tune = float(ikr_rescale)
        param_dict = {'tune.ikr_rescale': self._ikr_tune}
        self._parameters.update(param_dict)

    def reset_parameters(self):
        param_dict = {}
        for k in modelling.SD_details.SD_param_names + ['Kt', 'gKr']:
            param_dict[self.ikr_key_head + '.' + k] = \
                self.model.get(self.ikr_component.var(k)).eval()

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

        self.sim.pre(paces * self._cycle_length)
        self.initial_state = self.sim.state()

    def simulate(self, prepace=1000, save_signal=1, timestep=0.1,
                 log_var=None, reset=True):

        self._set_parameters(self._parameters)
        self.prepace = prepace

        self.sim.set_time(0)
        if reset:
            self.sim.reset()
            self.sim.set_state(self.initial_state)

        self.sim.pre(self._cycle_length * self.prepace)
        log = self.sim.run(self._cycle_length * save_signal, log=log_var,
                           log_interval=timestep).npview()
        if save_signal > 1:
            log = log.fold(self._cycle_length)

        return log

    def APD90(self, log):
        """
        Compute APD90
        """
        time_log = log.time()
        timestep = (max(time_log) - min(time_log)) / len(time_log)

        Vm_key_list = [x for x in log.keys() if x.endswith(str(self.Vm_key))]
        Vm_key_list.sort()
        pulse_num = len(Vm_key_list)

        Vm_signal = []
        for i in range(pulse_num):
            Vm_signal += list(log[Vm_key_list[i]])

        APD90s = []
        for i in range(pulse_num):
            AP_Vm = Vm_signal[
                round(i * self._cycle_length / timestep):
                round((i + 1) * (self._cycle_length + self.protocol_offset) /
                      timestep)]
            APamp = max(AP_Vm) - min(AP_Vm)
            APD90_v = min(AP_Vm) + 0.1 * APamp

            time_log = np.concatenate((time_log,
                                       time_log + max(time_log) + timestep))

            min_APD = int(5 / timestep)
            offset_index = int(self.protocol_offset / timestep)
            # offset index can be estimated from log

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

        if len(APD90s) == 1:
            return APD90s[0]
        else:
            return APD90s

    def qNet(self, log):
        """
        Compute qNet
        """
        time_log = log.time()
        timestep = (max(time_log) - min(time_log)) / len(time_log)

        # Make sure requirements to compute qNet are satisfied
        if np.abs(timestep - 0.01) > 1e-8:
            print('Warning: Time step should be 0.01ms instead of ' \
                  f'{timestep}ms.')
        if np.abs(time_log[-1] - time_log[0] - 1999.99) > 1e-8:
            print('Warning: qNet should be calculated with time length ' \
                  f'2000ms instead of {time_log[-1] - time_log[0]}ms.')
        if np.abs(self._cycle_length - 2000) > 1e-8:
            print(f'The pacing cycle length is set to {self._cycle_length}ms.')
            print('Warning: qNet should be calculated at 0.5Hz (cl=2000).')
        if self.prepace != 1000:
            print(f'The number of prepaces is set to {self.prepace}.')
            print('Warning: Dutta et al. 2017 used 1000 prepace.')

        qNet_current = ['ICaL', 'INaL', 'IKr', 'IKs', 'IK1', 'Ito']
        missing_current_keys = [i for i in qNet_current
                                if self.current_keys[i] is None]
        if missing_current_keys:
            print('There are missing currents')
            # Rephrase
        qNet_current_keys = [self.current_keys[i] for i in qNet_current
                             if self.current_keys[i] is not None]
        inet = 0
        for c in qNet_current_keys:
            inet += log[c]
            # Assuming the data log has only one period

        return np.trapz(inet, x=log.time()) * 1e-3  # pA/pF*ms -> pA/pF*s

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
