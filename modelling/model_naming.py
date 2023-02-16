class ModelDetails(object):
    """
    To create a library of all the current namings for different AP models.
    """

    def __init__(self):
        super(ModelDetails, self).__init__()

        self.APmodels = ['Grandi', 'TTP', 'AP-SD', 'Grandi-SD']
        self.current_list = ['IKr', 'IKs', 'Ito', 'IKb', 'IK1', 'INaK', 'INa',
                             'INaL', 'ICaL', 'INaCa', 'INaB', 'ICaB', 'IClCa',
                             'IClB', 'ICaP', 'IKACh', 'IKATP']
        self.current_keys = {
            'Grandi': {
                'IKr': 'I_Kr.I_kr',
                'IKs': 'I_Ks.I_ks',
                'Ito': 'I_to.I_to',
                'IKb': None,  # same as plateau potassium current?
                'IK1': 'I_Ki.I_ki',
                'INaK': 'I_NaK.I_nak',
                'INa': 'I_Na.I_Na',
                'INaL': None,
                'ICaL': 'I_Ca.I_Catot',
                'INaCa': 'I_NCX.I_ncx',
                'INaB': 'I_NaBK.I_nabk',
                'ICaB': 'I_CaBK.I_cabk',
                'IClCa': 'I_ClCa.I_ClCa',
                'IClB': 'I_ClCa.I_Clbk',
                'ICaP': 'I_PCa.I_pca',
                'IKACh': None,
                'IKATP': None
            },
            'AP-SD': {
                'IKr': 'ikr.IKr',
                'IKs': 'iks.IKs',
                'Ito': 'ito.Ito',
                'IKb': 'ikb.IKb',
                'IK1': 'ik1.IK1',
                'INaK': 'inak.INaK',
                'INa': 'ina.INa',
                'INaL': 'inal.INaL',
                'ICaL': 'ical.ICaL_total',
                'INaCa': 'inaca.INaCa_total',
                'INaB': 'inab.INab',
                'ICaB': 'icab.ICab',
                'IClCa': None,
                'IClB': None,
                'ICaP': 'ipca.IpCa',
                'IKACh': None,
                'IKATP': None
            },
        }
