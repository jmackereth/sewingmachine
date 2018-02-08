import numpy as np
import matplotlib.pyplot as plt
import linelist
from scipy.interpolate import interp1d

def measurelinelist(spec, line_obj,
                    sigmaclip=True,
                    sigma=2,
                    exclude_bad = True,
                    error = False,
                    plot=False,
                    verbose=False,
                    return_flags=False):
    '''
    measurelinelist
    ----------------
    Measure all EWs of lines in a linelist for a given spectrum. Plots optional.
    ----------------
    INPUT:
    spec - [N_lambda,2] (or [N_lambda,3]) shape array where [:,0] contains wavelengths,
    [:,1] is flux. For errors a 3rd column with the error spectrum should be included.
    line_obj - either a linelist.Linelist object, or a path to the desired linelist
    sigmaclip - if True, clip continuum pixels outside of sigma range
    sigma - N_sigma outside which continuum pixels are clipped
    exclude_bad - if True, exclude bad pixels (with flux = 0) from the continuum fit
    errors - if True, integrate an error spectrum, provided in the spec variable
    plot - plot all the measured EWs and continuum fits into successive axes
    verbose - if True print messages about the measurements
    return_flags - if True, flags are returned containing issues with measurement
    ----------------
    OUTPUT:
    EWs - array of measured EWs, in the order of input linelist
    flags (optional) - list of warning flags for each measurement
    ----------------
    HISTORY:
    2018/29/01 - Written - Mackereth (ARI, LJMU)
    '''
    if type(line_obj) is str:
        line_obj = linelist.Linelist(line_obj)
    nlines = len(line_obj.labels)
    EWs = np.empty(len(line_obj.labels))
    if return_flags:
        flags = []
    if error:
        errs = np.empty(len(line_obj.labels))
    if plot:
        fig =  plt.figure()
        fig.set_size_inches(7,3*nlines)
        width = 0.9
        height = 0.85/float(nlines)
        gap = 0.1/float(nlines)
    axes = []
    for i in range(nlines):
        if not plot:
            if return_flags:
                if error:
                    EWs[i], errs[i], flagi = trapz_ew(spec, line_obj.integration[i], line_obj.windows[i],
                                                      sigmaclip = sigmaclip, sigma = sigma, verbose=verbose,
                                                      exclude_bad=exclude_bad, error=True, return_flags=return_flags)
                else:
                    EWs[i], flagi = trapz_ew(spec, line_obj.integration[i], line_obj.windows[i],
                                            sigmaclip = sigmaclip, sigma = sigma, verbose=verbose,
                                            exclude_bad=exclude_bad, return_flags=return_flags)
                flags.append(flagi)
            else:
                if error:
                    EWs[i], errs[i] = trapz_ew(spec, line_obj.integration[i], line_obj.windows[i],
                              sigmaclip = sigmaclip, sigma = sigma, verbose=verbose,
                              exclude_bad=exclude_bad, error=True)
                else:
                    EWs[i] = trapz_ew(spec, line_obj.integration[i], line_obj.windows[i],
                              sigmaclip = sigmaclip, sigma = sigma, verbose=verbose,
                              exclude_bad=exclude_bad)
        if plot:
            ax = fig.add_axes([0.1,(nlines-i)*height+(nlines-i)*gap, width, height])
            plt.axes(ax)
            axes.append(ax)
            if return_flags:
                if error:
                    EWs[i], errs[i], flagi = trapz_ew(spec, line_obj.integration[i], line_obj.windows[i],
                                  sigmaclip = sigmaclip, sigma = sigma, verbose=verbose,plot=True,
                                  exclude_bad=exclude_bad, error=True, return_flags=return_flags)
                else:
                    EWs[i], flagi = trapz_ew(spec, line_obj.integration[i], line_obj.windows[i],
                                  sigmaclip = sigmaclip, sigma = sigma, verbose=verbose,plot=True,
                                  exclude_bad=exclude_bad, return_flags=return_flags)
                flags.append(flagi)
            else:
                if error:
                    EWs[i], errs[i] = trapz_ew(spec, line_obj.integration[i], line_obj.windows[i],
                                sigmaclip = sigmaclip, sigma = sigma, plot=True, verbose=verbose,
                                exclude_bad = exclude_bad, error=True)
                else:
                    EWs[i] = trapz_ew(spec, line_obj.integration[i], line_obj.windows[i],
                                sigmaclip = sigmaclip, sigma = sigma, plot=True, verbose=verbose,
                                exclude_bad = exclude_bad)
            ax.text(0.06,0.1,line_obj.labels[i]+r' EW $= $'+str(round(EWs[i],3))+r'$\mathrm{\ \AA}$', fontsize=10, transform=ax.transAxes,
                    bbox={'facecolor':'White', 'alpha':0.8, 'pad':5})
    if plot:
        for i in range(len(axes)):
            if i != 2:
                axes[i].set_xlabel('')
    if return_flags:
        if error:
            return EWs, errs, flags
        return EWs, flags
    if error:
        return EWs, errs
    return EWs

def trapz_ew(spec, integration, windows,
                  sigmaclip=True,
                  sigma=2,
                  exclude_bad = True,
                  error = False,
                  plot=False,
                  verbose=False,
                  return_flags = False):
    '''
    trapz_ew
    ----------------
    Measure an EW using a simple trapezium integration, including continuum fit
    and sigma clipping (optional)
    ----------------
    INPUT:
    spec - [N_lambda,2] shape array where [:,0] contains wavelengths and [:,1] is flux
    (optional third column for error array)
    integration - Length 2 tuple with integration region
    windows - list of 2-tuples defining continuum windows (any length)
    sigmaclip - if True, clip continuum pixels outside of sigma range
    sigma - N_sigma outside which continuum pixels are clipped
    exclude_bad - if True, exclude bad pixels (with flux = 0) from the continuum fit
    error - if True, integrate an error array (in 3rd column of spec)
    plot - plot the measured EW and continuum onto the current axes
    verbose - if True print messages about the measurement
    return_flags - if True, flags are returned containing issues with measurement
    ----------------
    OUTPUT:
    EW - the measured EW
    flags (optional) - list of warning flags
    ----------------
    HISTORY:
    2018/29/01 - Written - Mackereth (ARI, LJMU)
    '''
    all_clipped = False
    flags = []
    spec_x = spec[:,0] #create arrays for spectra
    spec_y = spec[:,1]
    norm_pix = []
    norm_lambda = []
    # mask out continuum pixels
    windowmask = np.zeros(len(spec_x), dtype=bool)
    all_clipped = False
    for i in windows:
        windowmask[(spec_x <= i[1]) & (spec_x >= i[0])] = 1
        if exclude_bad and len(spec_y[windowmask][spec_y[windowmask] < 0.0001]) > 0:
            bad_excluded = True
            if len(spec_y[windowmask][spec_y[windowmask] > 0.0001]) < 2:
                all_clipped = True
                flags.append('CONTINUUM_ALL_BAD')
                if verbose:
                    print('Continuum all at 0 flux - returning NaN EW')
            else:
                windowmask[spec_y < 0.0001] = 0
                flags.append('CONTINUUM_BAD_PIXEL')
    # do initial continuum fit
    cont_fit = np.polyfit(spec_x[windowmask],spec_y[windowmask], 1)
    cont_poly = np.poly1d(cont_fit)
    if sigmaclip:
        #clip pixels outside sigma stds of the continuum fit - then re-fit
        std = sigma*np.nanstd(spec_y[windowmask]- cont_poly(spec_x[windowmask]))
        residual = np.fabs(spec_y - cont_poly(spec_x))
        in_clip = residual <= std
        clipped_inds = np.where(windowmask & ~in_clip)
        windowmask = windowmask & in_clip
        if len(spec_x[windowmask]) < 2:
            if verbose:
                print('Bad continuum, returning NaN EW - consider re-defining windows?')
            if return_flags:
                flags.append('CONTINUUM_ALL_CLIPPED')
            all_clipped = True
        if not all_clipped:
            if verbose:
                print 'Re-fitting continuum after sigma clipping...'
            cont_fit = np.polyfit(spec_x[windowmask],spec_y[windowmask], 1)
            cont_poly = np.poly1d(cont_fit)
    # mask out integration region
    intmask = np.zeros(len(spec_x), dtype=bool)
    intmask[(spec_x <= integration[1]) & (spec_x >= integration[0])] = 1
    #linear interpolation of spectrum - for fractional pixels
    interpol = interp1d(spec_x, spec_y, kind='linear')
    line_x = np.zeros(len(spec_x[intmask])+2)
    line_y = np.zeros(len(spec_x[intmask])+2)
    # line inside fractional pixels
    line_x[1:-1] = spec_x[intmask]
    line_y[1:-1] = spec_y[intmask]
    #fractional pixels on blue and red side of line
    line_x[0] = integration[0]
    line_x[-1]= integration[1]
    line_y[0] = interpol(integration[0])
    line_y[-1] = interpol(integration[1])
    if return_flags and (line_y < 0.0001).any():
        flags.append('INTEGRATION_BAD_PIXEL')
    #continuum over the integration region
    cont_y = cont_poly(line_x)
    #normalise the flux over this range
    nline_y = line_y/cont_y
    contarea = max(line_x)-min(line_x)
    linearea = np.trapz((1-nline_y), x=line_x)
    if error:
        errspec = np.dstack([spec[:,0], spec[:,2]])[0]
        err = trapz_error(errspec, integration)
    if plot:
        # make plot of line and fitted continuum
        flat_windows = np.array(windows).ravel()
        minmax_windows = [np.min(flat_windows), np.max(flat_windows)]
        window_range = minmax_windows[1]-minmax_windows[0]
        plotrange = [minmax_windows[0]-window_range*0.05, minmax_windows[1]+window_range*0.05]
        plotmask = np.zeros(len(spec_x), dtype=bool)
        plotmask[(spec_x < plotrange[1]) & (spec_x > plotrange[0])] = 1
        plt.plot(spec_x[plotmask], spec_y[plotmask], zorder=1)
        for ii, window in enumerate(windows):
            plt.axvline(window[0], color='Black', linestyle='dotted')
            plt.axvline(window[1], color='Black', linestyle='dotted')
            plt.axvspan(window[0], window[1], color='Gray', alpha=0.1)
        plt.fill_between(line_x,line_y,cont_y, color='Black', alpha=0.2)
        plt.scatter(spec_x[windowmask], spec_y[windowmask], color='Black', s=20, lw=0.,zorder=2)
        if sigmaclip:
            plt.scatter(spec_x[clipped_inds], spec_y[clipped_inds], color='Red', s=30, lw=1., marker='x',zorder=2)
        plt.plot(spec_x[plotmask], cont_poly(spec_x[plotmask]), color='Black', linestyle='dashed')
        plt.xlabel(r'$\lambda [\mathrm{\AA}]$')
        plt.ylabel(r'$f/f_c(\lambda)$')
    if all_clipped:
        # if bad-continuum return NaN
        if return_flags:
            if error:
                return np.nan, np.nan, flags
            return np.nan, flags
        if error:
            return np.nan, np.nan
        return np.nan
    if return_flags:
        if error:
            return linearea, err, flags
        return linearea, flags
    if error:
        return linearea, err
    return linearea


def trapz_error(errspec, integration):
    intmask = np.zeros(len(errspec[:,0]), dtype=bool)
    intmask[(errspec[:,0] <= integration[1]) & (errspec[:,0] >= integration[0])] = 1
    interpol = interp1d(errspec[:,0], errspec[:,1], kind='linear')
    line_x = np.zeros(len(errspec[:,0][intmask])+2)
    line_y = np.zeros(len(errspec[:,1][intmask])+2)
    # line inside fractional pixels
    line_x[1:-1] = errspec[:,0][intmask]
    line_y[1:-1] = errspec[:,1][intmask]
    #fractional pixels on blue and red side of line
    line_x[0] = integration[0]
    line_x[-1]= integration[1]
    line_y[0] = interpol(integration[0])
    line_y[-1] = interpol(integration[1])
    err_line = np.sqrt(np.sum(line_y**2))
    return err_line
