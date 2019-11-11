from sewingmachine import equivalentwidths, linelist
import numpy as np
from astropy.io import fits
import apogee.tools.read as apread
import apogee.tools.path as appath
import apogee.spec.plot as splot
from tqdm import tqdm
import os

_DEFAULT_DR = appath._default_dr()


def measure_apogee(allStar, linelist_obj, output_fits=False, *args, **kwargs):
    if isinstance(linelist_obj, str):
        linelist_obj = linelist.Linelist(linelist_obj)
    loc_ids, apogee_ids = make_speclist(allStar)
    lams = splot.apStarWavegrid()
    ews = np.empty([np.shape(allStar)[0], np.shape(linelist_obj.labels)[0]])
    errs = np.empty([np.shape(allStar)[0], np.shape(linelist_obj.labels)[0]])
    try:
        dr = int(_DEFAULT_DR)
    except ValueError:
        dr = 16
    if dr <= 13:
        lockey = 'LOCATION_ID'
    else:
        lockey = 'FIELD'
    for i in tqdm(range(len(apogee_ids))):
        try:
            if isinstance(allStar[lockey][i], np.bytes_):
                specs, hdr = apread.aspcapStar(allStar[lockey][i].decode(), allStar['APOGEE_ID'][i].decode(), ext=1)
                errspec, hdr = apread.aspcapStar(allStar[lockey][i].decode(), allStar['APOGEE_ID'][i].decode(), ext=2)
            else:
                specs, hdr = apread.aspcapStar(allStar[lockey][i], allStar['APOGEE_ID'][i], ext=1)
                errspec, hdr = apread.aspcapStar(allStar[lockey][i], allStar['APOGEE_ID'][i], ext=2)
            spec = np.dstack([lams, specs, errspec])[0]
            out = equivalentwidths.measurelinelist(spec, linelist_obj, error=True, *args, **kwargs)
            ews[i], errs[i] = out[0], out[1]
        except IOError:
            print('Spectrum missing from SAS?')
            ews[i], errs[i] = np.ones(np.shape(linelist_obj.labels)[0])*np.nan, np.ones(np.shape(linelist_obj.labels)[0])*np.nan
    return ews, errs

def make_speclist(allStar):
    loc_ids = allStar['LOCATION_ID']
    apogee_ids = allStar['APOGEE_ID']
    return loc_ids, apogee_ids
