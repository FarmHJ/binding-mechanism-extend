class ModelDetails(object):
    """
    To create a library of all the current namings for different AP models.
    """

    def __init__(self):
        super(ModelDetails, self).__init__()

        self.APmodels = ['Grandi', 'TTP', 'ORd-CiPA', 'Grandi-SD']
        self.current_list = ['IKr', 'IKs', 'Ito', 'IKb', 'IK1', 'INaK', 'INa',
                             'INaL', 'ICaL', 'INaCa', 'INaB', 'ICaB', 'IClCa',
                             'IClB', 'ICaP', 'IKACh', 'IKATP']
        self.current_keys = {
            'ORd': {
                'Ito': 'Ito.Ito',
                'IKb': 'IKb.IKb',
                'IKs': 'IKs.IKs',
                'IKr': 'IKr.IKr',
                'ICaP': 'IpCa.IpCa',
                'IK1': 'IK1.IK1',
                'INaK': 'INaK.INaK',
                'ICaL': 'ICaL.ICaL_total',
                'INaL': 'INaL.INaL',
                'INaCa': 'INaCa_i.INaCa_total',
                'ICaB': 'ICab.ICab',
                'INaB': 'INab.INab',
                'INa': 'INa.INa',
                'IClCa': None,
                'IClB': None,
                'IKACh': None,
                'IKATP': None,
                'time': 'environment.time',
                'Vm': 'membrane.v'
            },
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
            'ORd-CiPA': {
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
                'IKATP': None,
                'time': 'engine.time',
                'Vm': 'membrane.V'
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
            'Tomek': {
                'IClB': 'ICl.IClb',
                'IClCa': 'ICl.IClCa',
                'Ito': 'Ito.Ito',
                'IKb': 'IKb.IKb',
                'IKs': 'IKs.IKs',
                'IKr': 'IKr.IKr',
                'ICaP': 'IpCa.IpCa',
                'IK1': 'IK1.IK1',
                'INaK': 'INaK.INaK',
                'ICaL': 'ICaL.ICaL',
                'INaL': 'INaL.INaL',
                'INaCa': 'INaCa.INaCa_total',
                'ICaB': 'ICab.ICab',
                'INaB': 'INab.INab',
                'INa': 'INa.INa',
                'IKACh': None,
                'IKATP': None,
                'time': 'environment.time',
                'Vm': 'membrane.v'
            },
            'Tomek-Cl': {
                'IClB': 'ICl.IClb',
                'IClCa': 'ICl.IClCa',
                'Ito': 'Ito.Ito',
                'IKb': 'IKb.IKb',
                'IKs': 'IKs.IKs',
                'IKr': 'IKr.IKr',
                'ICaP': 'IpCa.IpCa',
                'IK1': 'IK1.IK1',
                'INaK': 'INaK.INaK',
                'ICaL': 'ICaL.ICaL',
                'INaL': 'INaL.INaL',
                'INaCa': 'INaCa.INaCa_total',
                'ICaB': 'ICab.ICab',
                'INaB': 'INab.INab',
                'INa': 'INa.INa',
                'IKACh': None,
                'IKATP': None,
                'time': 'environment.time',
                'Vm': 'membrane.v'
            },
            'Lei': {
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
                'IKATP': None,
                'time': 'engine.time',
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

        self.file_names = {
            'ORd': {
                'AP_path': 'math_model/AP_model/ORd-2011.mmt',
                'AP_SD_path': 'math_model/AP_model/ohara-cipa-2017.mmt',
                'label': 'ORd (2011)'
            },
            'Grandi': {
                'AP_path': 'math_model/AP_model/Grd-2010.mmt',
                'AP_SD_path': 'math_model/AP_model/Grd-2010-IKr-SD.mmt',
                'label': 'Grandi (2010)'
            },
            'ORd-CiPA': {
                'AP_path': 'math_model/AP_model/ohara-cipa-2017.mmt',
                'label': "O'Hara-CiPA (2017)"
            },
            'TTP': {
                'AP_path': 'math_model/AP_model/TTP-2006.mmt',
                'AP_SD_path': 'math_model/AP_model/TTP-2006-IKr-SD.mmt',
                'label': 'ten Tusscher (2006)'
            },
            'Tomek': {
                'AP_path': 'math_model/AP_model/Tomek-2019.mmt',
                'AP_SD_path': 'math_model/AP_model/Tomek-2019-IKr-SD.mmt',
                'label': 'Tomek (2019)'
            },
            'Tomek-Cl': {
                'AP_path': 'math_model/AP_model/Tomek-Cl-2020.mmt',
                'AP_SD_path': 'math_model/AP_model/Tomek-Cl-2020-IKr-SD.mmt',
                'label': 'Tomek (2020)'
            },
            'Lei': {
                'IKr_path': 'math_model/current_model/lei2019.mmt',
                'AP_path': 'math_model/AP_model/ohara-cipa-2017.mmt',
                'AP_SD_path': 'math_model/AP_model/ORd-CiPA-Lei-SD.mmt',
                'label': 'Lei (2019)'
            },
        }


class SDModelDetails(object):
    """
    To create a library on constant and variable namings of the SD model
    """

    def __init__(self):
        super(SDModelDetails, self).__init__()

        self.param_names = ['Vhalf', 'Kmax', 'Ku', 'N', 'EC50']

        self.current_keys = {
            'Lei': {
                'IKr': 'ikr.IKr',
                'time': 'engine.time',
                'Vm': 'membrane.V',
            }}
