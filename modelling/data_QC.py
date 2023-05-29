import numpy as np


class QualityControl(object):
    """
    Methods to remove unwanted experimental traces:
    General
    1. Check seal resistance, membrane capacitance and series resistance
    For Pharm protocol
    1. Traces look like the protocol - 2nd segment higher than 3rd segment
    (control)
    """

    def __init__(self):
        super(QualityControl, self).__init__()

        # Define thresholds
        # General
        self.Rseal_thres = [1e8, 1e12]
        self.Cm_thres = [1e-12, 1e-10]
        # self.Rseries_thres = [1e6, 2.5e7]
        self.Rseries_thres = [5e6, 2e7]

        # Stability of traces
        self.rmsd0c = 0.2

    def qc_general(self, Rseal, Cm, Rseries):
        if any(np.array(Rseal) < self.Rseal_thres[0]) or \
                any(np.array(Rseal) > self.Rseal_thres[1]):
            # print('Rseal: ', Rseal)
            qc_Rseal = False
        else:
            qc_Rseal = True

        if any(np.array(Cm) < self.Cm_thres[0]) or \
                any(np.array(Cm) > self.Cm_thres[1]):
            # print('Cm: ', Cm)
            qc_Cm = False
        else:
            qc_Cm = True

        if any(np.array(Rseries) < self.Rseries_thres[0]) or \
                any(np.array(Rseries) > self.Rseries_thres[1]):
            # print('Rseries: ', Rseries)
            qc_Rseries = False
        else:
            qc_Rseries = True

        return [qc_Rseal, qc_Cm, qc_Rseries]

    def qc_stable(self, trace1, trace2):
        rmsd0_1 = np.sqrt(np.mean((trace1) ** 2))
        rmsd0_2 = np.sqrt(np.mean((trace2) ** 2))
        rmsdc = np.mean([rmsd0_1, rmsd0_2]) * self.rmsd0c

        rmsd_trace = np.sqrt(np.mean((trace1 - trace2) ** 2))
        if rmsd_trace > rmsdc or not \
                (np.isfinite(rmsd_trace) and np.isfinite(rmsdc)):
            print('rmsd: ', rmsd_trace)
            print('rmsd to zero: ', rmsdc)
            return False
        return True
