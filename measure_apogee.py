from sewingmachine import equivalentwidths, linelist
import numpy as np
from astropy.io import fits
import apogee.tools.read as apread
import apogee.spec.plot as splot

def measure_apogee(allStar, linelist_obj, output_fits=False):
    if isinstance(linelist_obj, str):
        linelist_obj = linelist.Linelist(linelist_obj)
    loc_ids, apogee_ids = make_speclist(allStar)
    lams = splot.apStarWavegrid()
    ews = np.empty([np.shape(allStar)[0], np.shape(linelist_obj.labels)[0]])
    errs = np.empty([np.shape(allStar)[0], np.shape(linelist_obj.labels)[0]])
    for i in tqdm(range(len(apogee_ids))):
        specs, hdr = apread.aspcapStar(allStar['LOCATION_ID'][i], allStar['APOGEE_ID'][i], ext=1)
        errspec, hdr = apread.aspcapStar(allStar['LOCATION_ID'][i], allStar['APOGEE_ID'][i], ext=2)
        spec = np.dstack([lams, specs, errspec])[0]
        ews[i], errs[i] = equivalentwidths.measurelinelist(spec, linelist_obj, error=True)
    return ews, errs        

def make_speclist(allStar):
    loc_ids = allStar['LOCATION_ID']
    apogee_ids = allStar['APOGEE_ID']
    return loc_ids, apogee_ids
