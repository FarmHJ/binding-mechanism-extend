import os
import pandas as pd

import modelling


def load_ikr_scaling(APmodel, tuning):
    scale_df = pd.read_csv(
        os.path.join(modelling.RESULT_DIR, 'ikr_calibration.csv'),
        index_col=[0], skipinitialspace=True)

    return scale_df.loc[APmodel, tuning]


APmodel_list = ['ORd-Li', 'Grandi', 'TTP', 'Tomek', 'ORd-Lei']

AP_file_names = {
    'ORd': {
        'AP': 'ORd-2011.mmt',
        'title': "ORd (2011)",
        'label': 'ORd'
    },
    'Grandi': {
        'AP': 'Grd-2010.mmt',
        'AP_IKr': 'Grd-2010-Li-SD.mmt',
        'title': 'Grandi (2010)',
        'label': 'Grandi-Li'
    },
    'ORd-Li': {
        'AP': 'ORd-2011.mmt',
        'AP_IKr': 'ohara-cipa-2017-opt.mmt',
        'title': "ORd-CiPAv1.0 (2017)",
        'label': 'ORd-Li'
    },
    'TTP': {
        'AP': 'TTP-2006.mmt',
        'AP_IKr': 'TTP-2006-Li-SD.mmt',
        'title': 'ten Tusscher (2006)',
        'label': 'ten Tusscher-Li'
    },
    'Tomek-old': {
        'AP': 'Tomek-2019.mmt',
        'AP_IKr': 'Tomek-2019-Li-SD.mmt',
        'title': 'Tomek (2019)',
        'label': 'Tomek-Li'
    },
    'Tomek': {
        'AP': 'Tomek-Cl-2020.mmt',
        'AP_IKr': 'Tomek-Cl-2020-Li-SD.mmt',
        'title': 'Tomek (2020)',
        'label': 'Tomek-Li'
    },
    'ORd-Lei': {
        'AP': 'ohara-cipa-2017-opt.mmt',
        # 'AP': 'ORd-2011.mmt',
        'AP_IKr': 'ORd-2011-Lei-SD.mmt',
        'title': 'Lei (2019)',
        'label': 'ORd-Lei'
    }}

IKr_file_names = {
    'Li': {
        'IKr-SD': 'li2016-SD.mmt',
        'label': "Li-SD"
    },
    'Lei': {
        'IKr-SD': 'lei2019-SD.mmt',
        'label': 'Lei-SD'
    }}

current_list = ['IKr', 'IKs', 'Ito', 'IKb', 'IK1', 'INaK', 'INa',
                'INaL', 'ICaL', 'INaCa', 'INaB', 'ICaB', 'IClCa',
                'IClB', 'ICaP', 'IKACh', 'IKATP']
qNet_current_list = ['ICaL', 'INaL', 'IKr', 'IKs', 'IK1', 'Ito']
model_current_keys = {
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
    'ORd-Li': {
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
    'Tomek-old': {
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
    'ORd-Lei': {
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

contribution_current_colours = dict({
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

qNet_current_colours = dict({
    'IKr': 0,
    'IKs': 0,
    'Ito': 0,
    'IKb': 19,
    'IK1': 0,
    'INaK': 19,
    'INa': 19,
    'INaL': 0,
    'ICaL': 0,
    'INaCa': 19,
    'INaB': 19,
    'ICaB': 19,
    'IClCa': 19,
    'IClB': 19,
    'ICaP': 19,
    'IKACh': 19,
    'IKATP': 19,
})

IKr_current_colours = dict({
    'IKr': 6,
    'IKs': 15,
    'Ito': 15,
    'IKb': 15,
    'IK1': 15,
    'INaK': 15,
    'INa': 15,
    'INaL': 15,
    'ICaL': 15,
    'INaCa': 15,
    'INaB': 15,
    'ICaB': 15,
    'IClCa': 15,
    'IClB': 15,
    'ICaP': 15,
    'IKACh': 15,
    'IKATP': 15,
})

current_labels = {
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
