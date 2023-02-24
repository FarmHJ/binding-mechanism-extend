import myokit
import numpy as np

import modelling


class Simulation(object):
    """
    To create a class to control the simulations of models.
    """

    def __init__(self, model, protocol=None, current_head_key=None):
        super(Simulation, self).__init__()

        self.model = model
        self.protocol = protocol
        self.sim = myokit.Simulation(self.model, self.protocol)
        self.initial_state = self.sim.state()
        self.current_head_key = current_head_key

        if current_head_key is None:
            self.current_head = next(iter(self.model.states())).parent()
        else:
            self.current_head = self.model.get(current_head_key)
        # Save model's original constants
        self.original_constants = {
            "Vhalf": self.model.get(self.current_head.var('Vhalf')).eval(),
            "Kmax": self.model.get(self.current_head.var('Kmax')).eval(),
            "Ku": self.model.get(self.current_head.var('Ku')).eval(),
            "n": self.model.get(self.current_head.var('n')).eval(),
            "EC50": self.model.get(self.current_head.var('halfmax')).eval(),
            "Kt": self.model.get(self.current_head.var('Kt')).eval(),
            "gKr": self.model.get(self.current_head.var('gKr')).eval(), }

    def model_simulation(self, repeats, conductance_name=None,
                         conductance_value=None, t_max=None,
                         timestep=0.1, save_signal=1, log_var=None,
                         abs_tol=1e-6, rel_tol=1e-4):

        if t_max is None:
            t_max = self.protocol.characteristic_time()

        self.sim = myokit.Simulation(self.model, self.protocol)
        self.sim.reset()
        self.sim.set_tolerance(abs_tol=abs_tol, rel_tol=rel_tol)

        if conductance_name is not None:
            self.sim.set_constant(conductance_name, conductance_value)

        self.sim.pre(t_max * (repeats - save_signal))
        log = self.sim.run(t_max * save_signal, log=log_var,
                           log_interval=timestep)
        d2 = log.npview()
        if save_signal > 1:
            d2 = d2.fold(t_max)

        return d2

    def drug_simulation(self, drug, drug_conc, repeats,
                        timestep=0.1, save_signal=1, log_var=None,
                        set_state=None, abs_tol=1e-6, rel_tol=1e-4):

        param_lib = modelling.BindingParameters()

        Vhalf = param_lib.binding_parameters[drug]['Vhalf']
        Kmax = param_lib.binding_parameters[drug]['Kmax']
        Ku = param_lib.binding_parameters[drug]['Ku']
        N = param_lib.binding_parameters[drug]['N']
        EC50 = param_lib.binding_parameters[drug]['EC50']

        t_max = self.protocol.characteristic_time()

        concentration = self.model.get(self.current_head_key + '.D')
        concentration.set_state_value(drug_conc)

        self.sim = myokit.Simulation(self.model, self.protocol)
        self.sim.reset()
        self.sim.set_tolerance(abs_tol=abs_tol, rel_tol=rel_tol)
        if set_state:
            set_state[self.current_head_key + '.D'][-1] = drug_conc
            self.sim.set_state(set_state)

        self.sim.set_constant(self.current_head.var('Vhalf'), Vhalf)
        self.sim.set_constant(self.current_head.var('Kmax'), Kmax)
        self.sim.set_constant(self.current_head.var('Ku'), Ku)
        self.sim.set_constant(self.current_head.var('n'), N)
        self.sim.set_constant(self.current_head.var('halfmax'), EC50)
        self.sim.set_constant(self.current_head.var('Kt'), 3.5e-5)
        self.sim.set_constant(self.current_head.var('gKr'),
                              self.original_constants["gKr"])
        print('conductance rescale: ',
              self.model.get('drug.ikr_rescale').eval())

        self.sim.pre(t_max * (repeats - save_signal))
        log = self.sim.run(t_max * save_signal, log=log_var,
                           log_interval=timestep)
        d2 = log.npview()
        if save_signal > 1:
            d2 = d2.fold(t_max)

        return d2

    def conductance_simulation(self, conductance, repeats,
                               timestep=0.1, save_signal=1, log_var=None,
                               abs_tol=1e-6, rel_tol=1e-4, set_state=None):
        self.sim = myokit.Simulation(self.model, self.protocol)
        self.sim.reset()
        if set_state:
            self.sim.set_state(set_state)
        self.sim.set_tolerance(abs_tol=abs_tol, rel_tol=rel_tol)

        self.sim.set_constant(self.current_head.var('Vhalf'),
                              self.original_constants["Vhalf"])
        self.sim.set_constant(self.current_head.var('Kmax'),
                              self.original_constants["Kmax"])
        self.sim.set_constant(self.current_head.var('Ku'),
                              self.original_constants["Ku"])
        self.sim.set_constant(self.current_head.var('n'),
                              self.original_constants["n"])
        self.sim.set_constant(self.current_head.var('halfmax'),
                              self.original_constants["EC50"])
        self.sim.set_constant(self.current_head.var('Kt'),
                              self.original_constants["Kt"])
        self.sim.set_constant(self.current_head.var('gKr'), conductance)
        t_max = self.protocol.characteristic_time()

        self.sim.pre(t_max * (repeats - save_signal))
        log = self.sim.run(t_max * save_signal, log=log_var,
                           log_interval=timestep)
        d2 = log.npview()
        if save_signal > 1:
            d2 = d2.fold(t_max)

        return d2

    def extract_peak(self, signal_log, current_name):
        peaks = []
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

    def APD90(self, signal, offset, timestep):
        APA = max(signal) - min(signal)
        APD90_v = min(signal) + 0.1 * APA
        min_APD = int(50 / timestep)
        offset_index = int(offset / timestep)
        index = np.abs(np.array(signal[offset_index + min_APD:]) -
                       APD90_v).argmin()
        APD90 = index * timestep + min_APD * timestep  # - offset
        # index = np.abs(np.array(signal[offset:]) - APD90_v).argmin()
        # APD90 = index * timestep - offset
        if APD90 < 1:
            print('check')
            APD90 = len(signal) * timestep

        return APD90, APD90_v
