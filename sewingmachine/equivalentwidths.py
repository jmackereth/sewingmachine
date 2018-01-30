import numpy as np
import matplotlib.pyplot as plt
import linelist
from scipy.interpolate import interp1d

def measurelinelist(spec, line_obj,
                    sigmaclip=True,
                    sigma=2,
                    plot=False,
                    verbose=False):
    '''
    measurelinelist
    ----------------
    Measure all EWs of lines in a linelist for a given spectrum. Plots optional.
    ----------------
    INPUT:
    spec - [N_lambda,2] shape array where [:,0] contains wavelengths and [:,1] is flux
    line_obj - either a linelist.Linelist object, or a path to the desired linelist
    sigmaclip - if True, clip continuum pixels outside of sigma range
    sigma - N_sigma outside which continuum pixels are clipped
    plot - plot all the measured EWs and continuum fits into successive axes
    verbose - if True print messages about the measurements
    ----------------
    OUTPUT:
    EWs - array of measured EWs, in the order of input linelist
    ----------------
    HISTORY:
    2018/29/01 - Written - Mackereth (ARI, LJMU)
    '''
    if type(line_obj) is str:
        line_obj = linelist.Linelist(line_obj)
    nlines = len(line_obj.labels)
    EWs = np.empty(len(line_obj.labels))
    if plot:
        fig =  plt.figure()
        fig.set_size_inches(7,3*nlines)
        width = 0.9
        height = 0.85/float(nlines)
        gap = 0.1/float(nlines)
    axes = []
    for i in range(nlines):
        if not plot:
            EWs[i] = trapz_ew(spec, line_obj.integration[i], line_obj.windows[i],
                          sigmaclip = sigmaclip, sigma = sigma, verbose=verbose)
        if plot:
            ax = fig.add_axes([0.1,(nlines-i)*height+(nlines-i)*gap, width, height])
            plt.axes(ax)
            axes.append(ax)
            EWs[i] = trapz_ew(spec, line_obj.integration[i], line_obj.windows[i],
                          sigmaclip = sigmaclip, sigma = sigma, plot=True, verbose=verbose)
            ax.text(0.06,0.1,line_obj.labels[i]+r' EW $= $'+str(round(EWs[i],3))+r'$\mathrm{\ \AA}$', fontsize=10, transform=ax.transAxes,
                    bbox={'facecolor':'White', 'alpha':0.8, 'pad':5})
    if plot:
        # hide the spines between ax and ax2
        for i in range(len(axes)):
            if i != 2:
                axes[i].set_xlabel('')
    return EWs

def trapz_ew(spec, integration, windows,
                  sigmaclip=True,
                  sigma=2,
                  plot=False,
                  verbose=False):
    '''
    trapz_ew
    ----------------
    Measure an EW using a simple trapezium integration, including continuum fit
    and sigma clipping (optional)
    ----------------
    INPUT:
    spec - [N_lambda,2] shape array where [:,0] contains wavelengths and [:,1] is flux
    integration - Length 2 tuple with integration region
    windows - list of 2-tuples defining continuum windows (any length)
    sigmaclip - if True, clip continuum pixels outside of sigma range
    sigma - N_sigma outside which continuum pixels are clipped
    plot - plot the measured EW and continuum onto the current axes
    verbose - if True print messages about the measurement
    ----------------
    OUTPUT:
    EW - the measured EW
    ----------------
    HISTORY:
    2018/29/01 - Written - Mackereth (ARI, LJMU)
    '''
    spec_x = spec[:,0] #create arrays for spectra
    spec_y = spec[:,1]
    norm_pix = []
    norm_lambda = []
    # mask out continuum pixels
    windowmask = np.zeros(len(spec_x), dtype=bool)
    for i in windows:
        windowmask[(spec_x <= i[1]) & (spec_x >= i[0])] = 1
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
        all_clipped = False
        if len(spec_x[windowmask]) < 2:
            if verbose:
                print('Bad continuum, returning NaN EW - consider re-defining windows?')
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
    #continuum over the integration region
    cont_y = cont_poly(line_x)
    #normalise the flux over this range
    nline_y = line_y/cont_y
    contarea = max(line_x)-min(line_x)
    linearea = np.trapz((1-nline_y), x=line_x)
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
        return np.nan
    return linearea
