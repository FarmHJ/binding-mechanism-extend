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
            # 'Grandi': {
            #     'IKr': 'I_Kr.I_kr',
            #     'IKs': 'I_Ks.I_ks',
            #     'Ito': 'I_to.I_to',
            #     'IKb': 'I_Kp.I_kp',  # same as plateau potassium current?
            #     'IK1': 'I_Ki.I_ki',
            #     'INaK': 'I_NaK.I_nak',
            #     'INa': 'I_Na.I_Na',
            #     'INaL': None,
            #     'ICaL': 'I_Ca.I_Catot',
            #     'INaCa': 'I_NCX.I_ncx',
            #     'INaB': 'I_NaBK.I_nabk',
            #     'ICaB': 'I_CaBK.I_cabk',
            #     'IClCa': 'I_ClCa.I_ClCa',
            #     'IClB': 'I_ClCa.I_Clbk',
            #     'ICaP': 'I_PCa.I_pca',
            #     'IKACh': None,
            #     'IKATP': None
            # },
            'Grandi': {
                'IClB': 'I_ClCa.I_Clbk',
                'IClCa': 'I_ClCa.I_ClCa',
                'Ito': 'I_to.I_to',
                'IKb': 'I_Kp.I_kp',  # same as plateau potassium current?
                'IKs': 'I_Ks.I_ks',
                'IKr': 'I_Kr.I_kr',
                'ICaP': 'I_PCa.I_pca',
                'IK1': 'I_Ki.I_ki',
                'INaK': 'I_NaK.I_nak',
                'ICaL': 'I_Ca.I_Catot',
                'INaCa': 'I_NCX.I_ncx',
                'ICaB': 'I_CaBK.I_cabk',
                'INaB': 'I_NaBK.I_nabk',
                'INa': 'I_Na.I_Na',
                'INaL': None,
                'IKACh': None,
                'IKATP': None,
                'time': 'environment.time',
                'Vm': 'membrane_potential.V_m',
            },
            'AP-SD': {
                'Ito': 'ito.Ito',
                'IKb': 'ikb.IKb',
                'IKs': 'iks.IKs',
                'IKr': 'ikr.IKr',
                'ICaP': 'ipca.IpCa',
                'IK1': 'ik1.IK1',
                'INaK': 'inak.INaK',
                'ICaL': 'ical.ICaL_total',
                'INaL': 'inal.INaL',
                'INaCa': 'inaca.INaCa_total',
                'ICaB': 'icab.ICab',
                'INaB': 'inab.INab',
                'INa': 'ina.INa',
                'IClCa': None,
                'IClB': None,
                'IKACh': None,
                'IKATP': None
            },
            'TTP': {
                'Ito': 'transient_outward_current.i_to',
                'IKb': 'potassium_pump_current.i_p_K',
                'IKs': 'slow_time_dependent_potassium_current.i_Ks',
                'IKr': 'rapid_time_dependent_potassium_current.i_Kr',
                'IK1': 'inward_rectifier_potassium_current.i_K1',
                'INaCa': 'sodium_calcium_exchanger_current.i_NaCa',
                'INaK': 'sodium_potassium_pump_current.i_NaK',
                'ICaP': 'calcium_pump_current.i_p_Ca',
                'ICaL': 'L_type_Ca_current.i_CaL',
                'ICaB': 'calcium_background_current.i_b_Ca',
                'INaB': 'sodium_background_current.i_b_Na',
                'INa': 'fast_sodium_current.i_Na',
                'INaL': None,
                'IClCa': None,
                'IClB': None,
                'IKACh': None,
                'IKATP': None,
                'time': 'environment.time',
                'Vm': 'membrane.V'
            },
        }

        self.current_colours = dict({
            'IKr': 0,
            'IKs': 1,
            'Ito': 2,
            'IKb': 3,
            'IK1': 4,
            'INaK': 5,
            'INa': 16,
            'INaL': 17,
            'ICaL': 10,
            'INaCa': 12,
            'INaB': 14,
            'ICaB': 15,
            'IClCa': 6,
            'IClB': 7,
            'ICaP': 13,
            'IKACh': 18,
            'IKATP': 19,
        })

        self.current_names = {
            'IKr': 'IKr',
            'IKs': 'IKs',
            'Ito': 'Ito (+Isus)',
            'IKb': 'IKb (IbK)',
            'IK1': 'IK1',
            'INaK': 'INaK',
            'INa': 'INa',
            'INaL': 'INaL',
            'ICaL': 'ICaL',
            'INaCa': 'INaCa (INCX)',
            'INaB': 'INa,B',
            'ICaB': 'ICa,B',
            'IClCa': 'IClCa',
            'IClB': 'ICl,B',
            'ICaP': 'ICa,P (IpCa)',
            'IKACh': 'IK,ACh',
            'IKATP': 'IK,ATP',
        }