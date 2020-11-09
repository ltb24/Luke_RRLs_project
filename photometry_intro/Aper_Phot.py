def Aper_Phot(data, channel, epoch, fwhm, sigma_val, round_val, plot):
    ## APERTURE PHOTOMETRY ##
    ## PERFORMS APERTURE PHOTOMETRY ON SPITZER MOSAIC IMAGES ##
    # data: input data to be analysed in FITS file data format
    # channel: Spitzer channel of desired image(s) - 1 (or 3.6um) or 2 (or 4.5um) allowed
    # circular_aperture: the aperture to be used in the photometry analysis
    # annuli_aperture: the annuli to be used in the photometry analysis
    
    global phot
    
    ## PARAMETERS ##
    channel = channel
    j = epoch
    fwhm = fwhm
    sigma_val = sigma_val
    roundlo = -round_val
    roundhi = round_val
    
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
    
    ## SOURCE DETECTION ##
    mean_val, median_val, std_val = sigma_clipped_stats(data, sigma = sigma_val)
    
    daofind = DAOStarFinder(fwhm = fwhm, threshold = sigma_val * std_val, roundlo = roundlo, roundhi = roundhi)
    sources = daofind(data) #- median_val) # necessary here?
    print('Number of stars detected: {}'.format(len(sources)))
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    apertures = CircularAperture(positions, r = 6.)

    if plot == True:
        plt.imshow(data, cmap = 'viridis', origin = 'lower', norm = LogNorm(), interpolation = 'nearest')
        plt.colorbar(fraction = 0.05)
        apertures.plot(color = 'black', lw = 1., alpha = .75)
        plt.title('Aperture photometry on epoch {} in channel {}: sigma = {}, fwhm = {}, roundlo = {}, roundhi = {}'
                  .format(j, channel, sigma_val, fwhm, roundlo, roundhi))
        plt.grid(b = True, which = 'major', lw = .5, alpha = .4, color = 'black')
        plt.gcf().set_size_inches(15, 8)
        plt.savefig(r'aperture_photometry_output_images/aperture_phot_epoch{}_channel{}.png'.format(j, channel))
        plt.show()
        plt.close()
    
    elif plot == False:
        pass
    
    ## APERTURE PHOTOMETRY ##
    
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    circular_apertures = CircularAperture(positions, r = 6.)
    annuli_apertures = CircularAnnulus(positions, r_in = 6., r_out = 14.)    
    apertures = [circular_apertures, annuli_apertures]

    # aperture photometry
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

    ## APPARENT MAGNITUDE ##

    phot = phot_init                    # redefine photometry table for ease
    phot['bkgsub_flux'] = float('NaN')  # populate new column to convert into flux
    phot['apparent_mag'] = float('NaN') # populate a new table (very quirky here, nans?)

    # convert data back into flux then calculate apparent magnitude
    for i in range(0, len(phot)):
        phot['bkgsub_flux'][i] = phot['aper_sum_bkgsub'][i] * fluxconv / exptime
        for i in range(0, len(phot)):
            if phot['bkgsub_flux'][i] >= 0:
                phot['apparent_mag'][i] = zmag - 2.5 * math.log10(phot['bkgsub_flux'][i] * aper_corr)

    # export into csv file
    phot['id', 'xcenter', 'ycenter', 'apparent_mag'].write(r'C:\Users\lukeb\Documents\MPhys_RRLs\output_files\aperphot_fn_epoch{}_{}.txt'.format(epoch, channel), format = 'csv', overwrite = True)
    return phot