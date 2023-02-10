import myokit
import numpy as np

import modelling


class Simulation(object):
    """
    To create a class to control the simulations of models.
    """

    def __init__(self, model, protocol=None):
        super(Simulation, self).__init__()

        self.model = model
        self.protocol = protocol
        self.sim = myokit.Simulation(self.model, self.protocol)
        self.initial_state = self.sim.state()

    def model_simulation(self, repeats, t_max=None,
                         timestep=0.1, save_signal=1, log_var=None,
                         abs_tol=1e-6, rel_tol=1e-4):

        if t_max is None:
            t_max = self.protocol.characteristic_time()

        self.sim = myokit.Simulation(self.model, self.protocol)
        self.sim.reset()
        self.sim.set_tolerance(abs_tol=abs_tol, rel_tol=rel_tol)

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
        index = np.abs(np.array(signal[offset + 5000:]) - APD90_v).argmin()
        print(index)
        APD90 = index * timestep - offset
        if APD90 < 1:
            APD90 = len(signal) * timestep

        return APD90, APD90_v
