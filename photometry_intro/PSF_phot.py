def PSF_phot(data, channel, epoch, fwhm, sigma_val, round_val, sharp_val, niters, plot):
    
    ## PARAMETERS ##
    channel = channel
    fwhm = fwhm
    sigma_val = sigma_val
    roundlo = -round_val
    roundhi = round_val
    sharohi = sharp_val
    ## COUNTER ##
    j = epoch

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
    
    ## SOURCE DETECTION ON ORIGINAL IMAGE ##

    mean_val, median_val, std_val = sigma_clipped_stats(data, sigma = sigma_val)
    
    psf_daofind = DAOStarFinder(fwhm = fwhm, threshold = sigma_val * std_val, roundlo = roundlo, roundhi = roundhi)
    psf_sources = psf_daofind(data)

    psf_positions = np.transpose((psf_sources['xcentroid'], psf_sources['ycentroid']))
    psf_apertures = CircularAperture(psf_positions, r = 6.)

    if plot == True:
        plt.imshow(data, cmap = 'viridis', origin = 'lower', norm = LogNorm(), interpolation = 'nearest')
        psf_apertures.plot(color = 'black', lw = 1.)
        plt.colorbar(fraction = 0.05)
        plt.title('Detected PSF stars for epoch{} in channel {}: threshold {} * std, fwhm = {}, roundlo = {}, roundhi = {}'
                  .format(j, channel, sigma_val, fwhm, roundlo, roundhi))
        plt.grid(b = True, which = 'major', lw = .5, color = 'black')
        plt.grid(b = True, which = 'minor', lw = .5, color = 'black')
        plt.gcf().set_size_inches(15, 8)
        plt.show()
    
    elif plot == False:
        pass
    
    print('Number of stars detected = {}'.format(len(psf_sources)))

    ## GROUP ##

    psf_sources['xcentroid'].name = 'x_0'
    psf_sources['ycentroid'].name = 'y_0'

    daogroup = DAOGroup(crit_separation = sigma_val * fwhm)
    bkg_estimator = MMMBackground()
    fitter = LevMarLSQFitter()

    data_psf = np.nan_to_num(data, nan = 1**-7)

    ## PHOTOMETRY ##

    PSF_photometry = IterativelySubtractedPSFPhotometry(finder = psf_daofind,
                                                        group_maker = daogroup,
                                                        bkg_estimator = bkg_estimator,
                                                        psf_model = epsf,
                                                        fitter = fitter,
                                                        niters = niters,
                                                        aperture_radius = 6.,
                                                        fitshape = (11, 11))

    result_phot = PSF_photometry(image = data_psf)
    residual_image = PSF_photometry.get_residual_image()
    print(len(result_phot))
    
    if plot == True:
        # visualise data
        plt.subplot(1, 2, 1)
        plt.imshow(data_psf, cmap = 'viridis', norm = LogNorm(), interpolation = 'nearest', origin = 'lower', vmin = 0.000001, vmax = 10**6)
        plt.title('input data')
        plt.colorbar(orientation = 'horizontal')

        plt.subplot(1, 2, 2)
        plt.imshow(residual_image, cmap = 'viridis', norm = LogNorm(), interpolation = 'nearest', origin = 'lower', vmin = 0.000001, vmax = 10**6)
        plt.title('residual image')
        plt.colorbar(orientation = 'horizontal')
        plt.gcf().set_size_inches(20, 14)
        plt.show()
        plt.close()
        
    elif plot == False:
        pass
    
    ## APPARENT MAGNITUDES ##

    phot = result_phot                  # redefine photometry table for ease
    phot['bkgsub_flux'] = float('NaN')  # populate new column to convert into flux
    phot['apparent_mag'] = float('NaN') # populate a new table (very quirky here, nans?)

    for i in range(0, len(phot)):
        phot['bkgsub_flux'][i] = phot['flux_0'][i] * fluxconv / exptime
        for i in range(0, len(phot)):
            if phot['bkgsub_flux'][i] >= 0:
                phot['apparent_mag'][i] = zmag - 2.5 * math.log10(phot['bkgsub_flux'][i] * aper_corr)

    # export into csv file
    phot['id', 'x_0', 'y_0', 'apparent_mag'].write(r'C:\Users\lukeb\Documents\MPhys_RRLs\output_files\psfphot_epoch{}_{}.txt'.format(j, channel), format = 'csv', overwrite = True)

    # format columns
    for col in phot.colnames:
        phot[col].info.format = '%.8g'
    print(phot['id', 'x_0', 'y_0', 'flux_0', 'apparent_mag'])