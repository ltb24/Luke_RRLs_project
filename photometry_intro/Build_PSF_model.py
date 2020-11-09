def Build_PSF_model(data, channel, j, fwhm, sigma_val, model_threshold, round_val, sharp_val, plot):
    
    ## PARAMETERS ##
    
    fwhm = fwhm
    sigma_val = sigma_val
    model_threshold = model_threshold
    roundlo = -round_val
    roundhi = round_val
    sharphi = sharp_val
    
    ## COUNTER ##
    j = j
    
    ## CHANNEL ##
    
    # aperture corrections for 337 (6,6,14) apertures in channels 1 & 2, given in IRAC handbook §4.10
    # zmags from IRAC handbook §4.8
    if channel == '3p6um':
        aper_corr = 1.125
        zmag = 18.8
    elif channel == '4p5um':
        aper_corr = 1.120
        zmag = 18.32
    else:
        return print('Input either channel 1 or channel 2 only')
    
    mean, median, std = sigma_clipped_stats(data, sigma = sigma_val)

    starfind_init = DAOStarFinder(fwhm = fwhm, threshold = model_threshold * std, roundlo = roundlo, roundhi = roundhi, sharphi = sharphi)
    epsf_sources = starfind_init(data)
    print('Number of model stars = {}'.format(len(epsf_sources)))
    
    if plot == True:

        # plot detected stars for ePSF model to verify good stars
        positions = np.transpose((epsf_sources['xcentroid'], epsf_sources['ycentroid']))
        apertures = CircularAperture(positions, r = 6.)

        plt.imshow(data, cmap = 'viridis', origin = 'lower', norm = LogNorm(), interpolation = 'nearest')
        apertures.plot(color = 'black', lw = 1.)
        plt.colorbar(fraction = 0.05)
        plt.title('Selected ePSF model stars for epoch{}, channel {}: threshold {} * std, fwhm = {}, roundlo = {}, roundhi = {}'
                  .format(j, channel, model_threshold, fwhm, roundlo, roundhi))
        plt.grid(b = True, which = 'major', lw = .5, color = 'black')
        plt.grid(b = True, which = 'minor', lw = .5, color = 'black')
        plt.gcf().set_size_inches(15, 8)
        plt.show()
        plt.close()
    
    elif plot == False:
        pass
    
    ## APERTURE PHOTOMETRY ON MODEL STARS ##
    
    positions = np.transpose((epsf_sources['xcentroid'], epsf_sources['ycentroid']))
    circular_apertures = CircularAperture(positions, r = 6.)
    annuli_apertures = CircularAnnulus(positions, r_in = 6., r_out = 14.)
    apertures = [circular_apertures, annuli_apertures]

    # initial aperture photometry table
    phot_init = aperture_photometry(data, apertures)

    # background subtraction using sigma-clipped median and annuli
    annulus_masks = annuli_apertures.to_mask(method = 'center')

    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0] # extract 1D array of data values
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d) # utilise sigma clipping on the annulus masks
        bkg_median.append(median_sigclip)

    bkg_median = np.array(bkg_median)
    # now append bkg_median, aperture background and aperture sum background values to photometry data
    phot_init['annulus_median'] = bkg_median
    phot_init['aper_bkg'] = bkg_median * circular_apertures.area
    phot_init['aper_sum_bkgsub'] = phot_init['aperture_sum_0'] - phot_init['aper_bkg']
    print(phot_init)
    
    ## STAR CUTOUTS FOR ePSF ##

    cutout_size = 200
    hsize = (cutout_size - 1) / 2
    x = epsf_sources['xcentroid']
    y = epsf_sources['ycentroid']
    mask = ((x > hsize) & (x < (data.shape[1] - 1 - hsize)) &
           (y > hsize) & (y < (data.shape[0] - 1 - hsize)))

    # table of star positions
    star_tbl = Table()
    star_tbl['x'] = x[mask]
    star_tbl['y'] = y[mask]
    print('Number of refined model stars = {}'.format(len(star_tbl)))

    if plot == True:
        # visualise stars to verify - not strictly necessary as visual checks should be done before using this script
        cutout_pos = np.transpose((star_tbl['x'], star_tbl['y']))
        cutout_apers = CircularAperture(cutout_pos, r = 6.)

        plt.imshow(data, cmap = 'viridis', origin = 'lower', norm = LogNorm(), interpolation = 'nearest')
        cutout_apers.plot(color = 'black', lw = 1.)
        plt.colorbar(fraction = 0.05)
        plt.title('Selected cutout ePSF model stars with params: threshold {} * std, fwhm = {}, roundlo = {}, roundhi = {}'
                  .format(model_threshold, fwhm, roundlo, roundhi))
        plt.grid(b = True, which = 'major', lw = .5, color = 'black')
        plt.grid(b = True, which = 'minor', lw = .5, color = 'black')
        plt.gcf().set_size_inches(15, 8)
        plt.show()
        plt.close()
        
    elif plot == False:
        pass
    
    ## EXTRACT STARS ##

    mean_val, median_val, std_val = sigma_clipped_stats(data, sigma = sigma_val)
    sub_data = data - median_val

    nddata = NDData(data = sub_data)
    stars = extract_stars(nddata, star_tbl, size = 25)
    
    if plot == True:
        # visualise 36 extracted stars
        nrows = 3
        ncols = 3
        fig, ax = plt.subplots(nrows = nrows, ncols = ncols, figsize = (20,20), squeeze = True)
        ax = ax.ravel()
        for i in range(nrows * ncols):
            norm = simple_norm(stars[i], 'log', percent = 99.)
            ax[i].imshow(stars[i], cmap = 'viridis', norm = norm, origin = 'lower')
        plt.show()
        plt.close()
    
    elif plot == False:
        pass

    ## BUILD ePSF ##
    global epsf, fitter

    epsf_builder = EPSFBuilder(oversampling = 2, maxiters = 10, progress_bar = True)
    epsf, fitter = epsf_builder(stars)

    norm = simple_norm(epsf.data, 'log', percent = 99.)
    plt.imshow(epsf.data, cmap = 'viridis', norm = norm, origin = 'lower')
    plt.colorbar()
    plt.show()
    plt.close()
    
    return epsf, fitter