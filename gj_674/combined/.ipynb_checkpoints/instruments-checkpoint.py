"""
Instrument bitwise codes for Mega-MUSCLES, based on sections from Parke Loyd's rc.py



"""


#listed in normalization order
instruments = ['hst_cos_g130m','hst_cos_g160m','hst_cos_g230l','hst_sts_g140m','hst_sts_e140m','hst_sts_e230m',
               'hst_sts_e230h','hst_sts_g230l','hst_sts_g430l','hst_sts_g430m','mod_gap_fill-',
               'xmm_epc_multi','xmm_epc_pn---', 'cxo_acs_-----', 'mod_euv_young', 'mod_apc_-----',
               'mod_lya_young','mod_phx_-----', 'oth_---_other', 'hst_sts_g230lb', 'hst_sts_g750l', 'hst_fos_g570h',
               'hst_fos_g780h', 'hst_cos_g140l', 'hst_sts_g140l']
instvals = [2**i for i in range(len(instruments))]
default_order = ['hst_cos_g130m','hst_cos_g160m','hst_cos_g230l','hst_sts_g140m','hst_sts_e140m','hst_sts_e230m',
                 'hst_sts_e230h','hst_sts_g230l', 'hst_sts_g230lb', 'xmm_epc_multi','xmm_epc_pn---', 'cxo_acs_-----',
                 'mod_euv_young', 'mod_apc_-----', 'mod_lya_young', 'mod_phx_-----', 'hst_sts_g750l', 'hst_sts_g430l',
                 'hst_sts_g430m', 'mod_gap_fill-', 'oth_---_other']

# for use in making plots
#instranges = {'xobs': [5., 60.], 'xmm': [5., 60.], 'cxo': [1.0, 2.0], 'euv': [100., 1170.], 'hst': [1170., 5700.],
              #'apec': [60., 100.], 'phx': [5700., 55000.], 'lya': lyacut, 'c130m': [1170., 1460.],
       #       'c160m': [1390., 1785.], 'c230l': [1670., 3190.], 's140m': [1170., 1710.], 's230m': [1605., 3070.],
        #      's230h': [2380., 2890.], 's230l': [1570., 3190.], 's430m': [3795., 4080.], 's430l': [2895., 5710.]}
instnames = {'xobs': 'XMM or Chandra', 'xmm': 'XMM', 'cxo': 'Chandra', 'euv': 'Empirical EUV Estimate',
             'hst': 'HST', 'apec': 'APEC Model Corona', 'phx': 'PHOENIX Model Photosphere',
             'lya': 'Ly$\\alpha$ Reconstruction', 'c130m': 'HST COS G130M', 'c160m': 'HST COS G160M',
             'c230l': 'HST COS G230L', 's140m': 'HST STIS E140M', 's230m': 'HST STIS E230M',
             's230h': 'HST STIS E230H', 's230l': 'HST STIS G230L', 's430m': 'HST STIS G430M',
             's430l': 'HST STIS G430L', 'acs':'ACIS'}
instabbrevs = {'xobs':'XMM or Chandra', 'apec':'APEC', 'euv':'Empirical EUV Estimate', 'hst':'HST', 'phx':'PHOENIX'}

# for use in making FITS headers
HLSPtelescopes = {'hst':'HST', 'cxo':'CXO', 'xmm':'XMM', 'mod':'MODEL', 'oth':'OTHER'}
HLSPinstruments = {'cos':'COS', 'sts':'STIS', 'euv':'EUV-SCALING', 'lya':'LYA-RECONSTRUCTION', 'phx':'PHX',
                   'epc':'EPIC', 'gap':'POLYNOMIAL-FIT', 'apc':'APEC', '---':'NA', 'acs':'ACIS', 'fos':'FOS'}
HLSPgratings = {'g130m':'G130M', 'g160m':'G160M', 'g430l':'G430L', 'g430m':'G430M', 'g140m':'G140M', 'e230m':'E230M',
                'e230h':'E230H', 'g230l':'G230L', 'e140m':'E140M', 'fill-':'NA', '-----':'NA', 'young':'NA',
                'pn---':'PN', 'multi':'MULTI', 'other':'NA', 'g750l':'G750L', 'g230lb':'G230LB', 'g570h':'g570H',
                'g780h':'G780H', 'g140l':'G140L'}


def getinststr(inst_val):
    """Return the string version of an instrument value."""
    return instruments[instvals.index(inst_val)]

def getinsti(instrument):
    """Return the identifying number for instrument, where instrument is
    of the form, e.g., 'hst_cos_g130m'."""
    try:
        return instvals[instruments.index(instrument)]
    except ValueError:
        return -99
