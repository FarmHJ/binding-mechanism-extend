class QualityControl(object):
    """
    Methods to remove unwanted experimental traces:
    General
    1. Check seal resistance, membrane capacitance and series resistance
    For Pharm protocol
    1. Traces look like the protocol - 2nd segment higher than 3rd segment (control)
    """

    def __init__(self):
        super(QualityControl, self).__init__()

        # Define thresholds
        # General
        self.Rseal_thres = [1e8, 1e12]
        self.Cm_thres = [1e-12, 1e-10]
        self.Rseries_thres = [1e6, 2.5e7]

    def qc_general(self, Rseal, Cm, Rseries):
        if Rseal < self.Rseal_thres[0] or Rseal > self.Rseal_thres[1]:
            print('Rseal: ', Rseal)
            qc_Rseal = False
        else:
            qc_Rseal = True

        if Cm < self.Cm_thres[0] or Cm > self.Cm_thres[1]:
            print('Cm: ', Cm)
            qc_Cm = False
        else:
            qc_Cm = True

        if Rseries < self.Rseries_thres[0] or Rseries > self.Rseries_thres[1]:
            print('Rseries: ', Rseries)
            qc_Rseries = False
        else:
            qc_Rseries = True

        return [qc_Rseal, qc_Cm, qc_Rseries]
