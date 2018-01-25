
def specewmeasure(spec, integration, windows, synthetic=False, sigmaclip= True, sigma=2):
    spec_x = spec[:,0] #create arrays for spectra
    spec_y = spec[:,1]
    norm_pix = []
    norm_lambda = []
    for i in windows:
        region = (spec_x > i[0]) & (spec_x < i[1])
        pix = spec_y[region]
        lambd = spec_x[region]
        norm_pix.extend(pix)
        norm_lambda.extend(lambd)
    norm_pix = np.array(norm_pix)
    norm_lambda = np.array(norm_lambda)
    cont_fit = np.polyfit(norm_lambda,norm_pix, 1)
    cont_poly = np.poly1d(cont_fit)
    if sigmaclip == True:
        std = sigma*np.std(norm_pix)
        residual = np.abs(norm_pix - cont_poly(norm_lambda))
        clip = residual < std
        norm_lambda = norm_lambda[clip]
        norm_pix = norm_pix[clip]
        cont_fit = np.polyfit(norm_lambda,norm_pix, 1)
        cont_poly = np.poly1d(cont_fit)
    int_reg = [(spec_x >= integration[0]) & (spec_x <= integration[1])]
    interpol = interp1d(spec_x, spec_y, kind='linear')
    line_x = spec_x[int_reg]
    line_y = spec_y[int_reg]
    bfrac_x=np.array([integration[0], line_x[0]])
    rfrac_x=np.array([line_x[-1], integration[1]])
    cont_y = cont_poly(line_x)
    norm_y = line_y/cont_y
    bfrac_y=np.array([interpol(integration[0])/cont_poly(integration[0]), norm_y[0]])
    rfrac_y=np.array([norm_y[-1], interpol(integration[1])/cont_poly(integration[1])])
    contarea = max(line_x)-min(line_x)
    linearea = simps((1-norm_y), line_x)
    bfrac_area = simps((1-bfrac_y), bfrac_x)
    rfrac_area = simps((1-rfrac_y), rfrac_x)
    if synthetic == False:
        err_interpol = interp1d(errspec_x, errspec_y, kind='linear')
        err_x = errspec_x[int_reg]
        err_y = errspec_y[int_reg]
        errbfrac_x = np.array([integration[0], err_x[0]])
        errrfrac_x = np.array([err_x[-1], integration[1]])
        errbfrac_y = np.array([err_interpol(integration[0]), err_y[0]])
        errrfrac_y = np.array([err_y[-1], err_interpol(integration[1])])
        err_line = simps(err_y, err_x)+ simps(errbfrac_y, errbfrac_x)+ simps(errrfrac_y, errrfrac_x)	
        EW = linearea+bfrac_area+rfrac_area
        return [EW, err_line]
    else:
        EW = linearea+bfrac_area+rfrac_area
    return [EW, 0]