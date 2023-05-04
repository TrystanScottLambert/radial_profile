import sys

sys.path.append('../../py-genfunctions/')

import tool_images as tim
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.centroids import centroid_com, centroid_quadratic
from photutils.centroids import centroid_1dg, centroid_2dg
import matplotlib.pylab as plt
from matplotlib.patches import Ellipse
from photutils.aperture import EllipticalAperture
from photutils.aperture import EllipticalAnnulus
from photutils import aperture_photometry

import warnings
# warnings.simplefilter('ignore', UserWarning)
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

def photometry_JWST(data, wcs, pixel_scale, radius, position, type_aper, theta_galaxy, verbose=False):
    
    # Convert position from world to pixel coordinates
    position_pix = wcs.wcs_world2pix(position[0], position[1], 1)
    
    
    # Convert radius from arcseconds to pixels
    radius = np.array(radius) 
    radius_pix = radius / pixel_scale

    
    # Create aperture
    if type_aper == 'ellipse':
        aper = EllipticalAperture((position_pix[0], position_pix[1]), radius_pix[0], radius_pix[1], theta=theta_galaxy)
    elif type_aper == 'annulus':
        aper = EllipticalAnnulus((position_pix[0], position_pix[1]), radius_pix[0], radius_pix[1], radius_pix[2], radius_pix[3], theta=theta_galaxy)
    
        
    # Create mask and apply it to data
    masks = aper.to_mask(method='center')
    data_weighted = masks.multiply(data)

    # Plot image if verbose is True
    if verbose:
        plt.figure()
        plt.imshow(data_weighted, origin='lower', vmin=-3*sigma_clipped_stats(data, mask=np.isnan(data))[2], vmax=15*sigma_clipped_stats(data, mask=np.isnan(data))[2])

    # Perform aperture photometry
    phot_table = aperture_photometry(data, aper, wcs=wcs)
    
    # Return the total flux in the aperture
    return float(phot_table['aperture_sum'])

def error_sigma_JWST(n_pix, rms, flux_eletron_s, A):

    error_aperture = np.sqrt(n_pix * rms + np.abs(flux_eletron_s))
    return error_aperture / A


def sigma_JWST(data, header, wcs, pixel_scale, radius, position, type_aper, theta_galaxy, verbose=False):

    counts = photometry_JWST(data, wcs, pixel_scale, radius, position, type_aper, theta_galaxy)

    flux_eletron_s = counts

    if type_aper == 'ellipse':
        theta = radius
        A = (np.pi * radius[0] * radius[1])
        sigma = flux_eletron_s / A

    if type_aper == 'annulus':
        A = np.pi * (radius[1] * radius[2] - radius[0] * radius[3])
        sigma = flux_eletron_s / A

    n_pix = A / (pixel_scale ** 2)
    rms_calculate_aperture = sigma_clipped_stats(data, mask=np.isnan(data))[2]

    error_sigma = error_sigma_UV(n_pix, rms_calculate_aperture *38184.516, flux_eletron_s *38184.516, A)

    #conversion_mjy = (header['PHOTFLAM'] * ((header['PHOTPLAM']) ** 2) * 3.33564095E+04) * 1e3
    
    error_sigma = error_sigma/error_sigma
    
    sigma_mjy = sigma #* conversion_mjy
    error_sigma_mjy = error_sigma #* conversion_mjy

    if verbose:
        plt.figure()
        plt.imshow(data_weighted, origin='lower', vmin=-3 * sigma_clipped_stats(data, mask=np.isnan(data))[2],
                   vmax=15 * sigma_clipped_stats(data, mask=np.isnan(data))[2])

    return sigma_mjy, error_sigma_mjy






def run_JWST(path_JWST,path_ALMA,b_array,eccen_galaxy,bpa_galaxy, show_diagnosis = True,use_convolved_image= True):#path_ALMA,path_HST, b_array, position, bmaj_galaxy, bmin_galaxy, bpa_galaxy, verbose = False):
    
    
    
    # Importing ALMA data:
    data_ALMA, header_ALMA, wcs_ALMA, pix_scale_ALMA, bmaj_ALMA, bmin_ALMA, bpa_ALMA, rms_ALMA, imagesize_ALMA = tim.import_fits_ALMA(path_ALMA[0])

    for file in range(len(path_JWST)):#len(path_HST[3:])):

        # Taking HST data
        name_JWST =  'F770W' 

        
        # Importing HST data:
        data_JWST, header_JWST, wcs_JWST, pixscale_JWST,  rms_JWST, imagesize_JWST, fhwm_JWST, fhwm_angle_JWST = tim.import_fits_JWST(path_JWST[file])


        
        # Convolve both images:
        convolved_ALMA, convolved_JWST = tim.convolution(data_HST_toconvolve     =  data_JWST,
                                                         data_ALMA_toconvolve     =  data_ALMA, 
                                                         beam_ALMA_toconvolve     =  (bmaj_ALMA, bmin_ALMA, bpa_ALMA), 
                                                         fhwm_HST_toconvolve     =  (fhwm_JWST,fhwm_JWST, fhwm_angle_JWST),
                                                         pixscale_HST_toconvolve  =  pixscale_JWST, 
                                                         pixscale_ALMA_toconvolve =  pix_scale_ALMA)
        
        rms_JWST_convolved = sigma_clipped_stats(convolved_JWST, mask = np.isnan(convolved_JWST))[2]
        
        # Calculating centroids of original image:
        ra_center_orig_pix, dec_center_orig_pix = tim.centroid(data            = data_JWST, 
                                                               x_init          = int(imagesize_JWST/2 + (0.1/pixscale_JWST)), 
                                                               y_init          = int(imagesize_JWST/2 + (0./pixscale_JWST)),
                                                               box_size        = int(1.5/pixscale_JWST), 
                                                               centroid_func   = centroid_quadratic, 
                                                               rms             = rms_JWST , 
                                                               name_cent_func  = ' ', 
                                                               verbose_output  = False)
        
        # Converting centroid in pixels to world coordinates
        ra_center_orig, dec_center_orig = wcs_JWST.all_pix2world(ra_center_orig_pix[0], dec_center_orig_pix[0],1)

        # Calculating centroids of convolved image:
        ra_center_conv_pix, dec_center_conv_pix = tim.centroid(data            = convolved_JWST, 
                                                               x_init          = int(imagesize_JWST/2 + (0.1/pixscale_JWST)), 
                                                               y_init          = int(imagesize_JWST/2 + (0.1/pixscale_JWST)),
                                                               box_size        = int(1.5/pixscale_JWST), 
                                                               centroid_func   = centroid_quadratic, 
                                                               rms             = rms_JWST_convolved, 
                                                               name_cent_func  = ' ', 
                                                               verbose_output  = False)
        
        
        #print('ra_center_conv_pix: ',ra_center_conv_pix, 'dec_center_conv_pix: ', dec_center_conv_pix)
        ra_center_conv, dec_center_conv = wcs_JWST.all_pix2world(ra_center_conv_pix, dec_center_conv_pix,1) 
        #print('ra_center_conv: ',ra_center_conv, 'dec_center_conv: ', dec_center_conv)
        
        #Diagnosis
        if show_diagnosis == True:
            plot_diagnosis(name_HST,data_HST,wcs_HST,rms_HST,ra_center_orig_pix, dec_center_orig_pix, imagesize_HST,pixscale_HST,eccen_galaxy,bpa_galaxy, convolved_HST,  rms_HST_convolved, ra_center_conv_pix, dec_center_conv_pix, b_array)
        
        if use_convolved_image == True:
            data = convolved_JWST
            position = (ra_center_conv[0], dec_center_conv[0])
            m = 'p'
            l = name_JWST + ' - convolved'
        else:
            data = data_JWST
            position = (ra_center_orig, dec_center_orig)
            m = 'X'
            l = name_JWST
        
        sigma_array = np.array([name_JWST])
        sigma_array_error = np.array([name_JWST])
        
        for rad in range(len(b_array)):
            if rad == 0:


                major_a = b_array[rad]
                minor_b = b_array[rad]*np.sqrt(1 - eccen_galaxy **2)

                array_radius = [major_a,minor_b]

                t_aper = 'ellipse'
                #print(minor_b)
            else:

                major_a_in = b_array[rad - 1]
                minor_b_in = b_array[rad - 1]*np.sqrt(1 - eccen_galaxy **2)
                major_a_out = b_array[rad]
                minor_b_out = b_array[rad]*np.sqrt(1 - eccen_galaxy **2)

                array_radius = [major_a_in ,major_a_out, minor_b_out,minor_b_in]

                t_aper = 'annulus'
            
            
            
            sig, err = sigma_JWST(data           = data,
                                header           = header_JWST,
                                wcs              = wcs_JWST,
                                pixel_scale      = pixscale_JWST,
                                radius           = array_radius,
                                position         = position,
                                type_aper        = t_aper,
                                theta_galaxy     = np.pi*(bpa_galaxy)/180) # <<<<< AQUI O AJUSTE INCLINACAO
            

            sigma_array = np.append(sigma_array ,sig)
            sigma_array_error = np.append(sigma_array_error ,err)
        
        if file == 0:
            sigma_output_array = sigma_array 
            sigma_output_array_error = sigma_array_error 
        else:
            sigma_output_array = np.vstack((sigma_output_array, sigma_array) ) #
            sigma_output_array_error =  np.vstack((sigma_output_array_error, sigma_array_error) ) # np.array([sigma_output_array_error, sigma_array_error] ,dtype=object)#
             
    return sigma_output_array, sigma_output_array_error






def photometry_HST(data, wcs, pixel_scale, radius, position, type_aper, theta_galaxy, verbose=False):
    """
    Performs aperture photometry on HST data.

    Parameters
    ----------
    data : numpy array
        HST data to be analyzed.
    wcs : WCS object
        WCS information for the data.
    pixel_scale : float
        Pixel scale in arcseconds/pixel.
    radius : array-like
        Radii of the aperture. For an elliptical aperture, radius should be a 2-element array [a,b] 
        with a and b being the semi-major and semi-minor axes, respectively. For an annular aperture, 
        radius should be a 4-element array [a_in, a_out, b_in, b_out] with a_in and b_in being the inner
        semi-major and semi-minor axes and a_out and b_out being the outer semi-major and semi-minor axes.
    position : tuple
        Position of the aperture center in world coordinates (RA, Dec).
    type_aper : str
        Type of aperture. Can be 'ellipse' or 'annulus'.
    theta_galaxy : float
        Position angle of the galaxy, in degrees. Used for elliptical apertures.
    verbose : bool, optional
        If True, plots the image after data weighting. Default is False.

    Returns
    -------
    float
        Total flux within the aperture.
    """
    
    # Convert position from world to pixel coordinates
    position_pix = wcs.wcs_world2pix(position[0], position[1], 1)
    
    
    # Convert radius from arcseconds to pixels
    radius = np.array(radius) 
    radius_pix = radius / pixel_scale

    
    # Create aperture
    if type_aper == 'ellipse':
        aper = EllipticalAperture((position_pix[0], position_pix[1]), radius_pix[0], radius_pix[1], theta=theta_galaxy)
    elif type_aper == 'annulus':
        aper = EllipticalAnnulus((position_pix[0], position_pix[1]), radius_pix[0], radius_pix[1], radius_pix[2], radius_pix[3], theta=theta_galaxy)

        
    # Create mask and apply it to data
    masks = aper.to_mask(method='center')
    data_weighted = masks.multiply(data)

    # Plot image if verbose is True
    if verbose:
        plt.figure()
        plt.imshow(data_weighted, origin='lower', vmin=-3*sigma_clipped_stats(data, mask=np.isnan(data))[2], vmax=15*sigma_clipped_stats(data, mask=np.isnan(data))[2])

    # Perform aperture photometry
    phot_table = aperture_photometry(data, aper, wcs=wcs)

    # Return the total flux in the aperture
    return float(phot_table['aperture_sum'])


def error_sigma_UV(n_pix, rms, flux_eletron_s, A):
    """
    Calculates the error in flux measurements in UV data.

    Parameters
    ----------
    n_pix : float
        Number of pixels in the aperture.
    rms : float
        RMS noise in the data.
    flux_electron_s : float
        Flux of electrons per second.
    A : float
        Aperture area.

    Returns
    -------
    float
        Error in the flux measurement.
    """
    error_aperture = np.sqrt(n_pix * rms + np.abs(flux_eletron_s))
                             
    return error_aperture / A



def sigma_UV(data, header, wcs, pixel_scale, radius, position, type_aper, theta_galaxy, verbose=False):
    """
    Calculates the flux and error in UV data.

    Parameters
    ----------
    data : ndarray
        Array containing the UV data.
    header : fits header
        FITS header containing relevant information for flux conversion.
    wcs : WCS object
        World Coordinate System object for the data.
    pixel_scale : float
        Pixel scale of the data.
    radius : tuple
        Tuple containing the radii of the aperture or annulus.
    position : tuple
        Tuple containing the position of the aperture or annulus.
    type_aper : str
        Type of aperture ('ellipse' or 'annulus').
    theta_galaxy : float
        Angle of rotation of the aperture in degrees.
    verbose : bool, optional
        Whether or not to display a plot of the weighted data, by default False.

    Returns
    -------
    tuple
        Tuple containing the flux and error in MJy/sr units.
    """
    counts = photometry_HST(data, wcs, pixel_scale, radius, position, type_aper, theta_galaxy)
    flux_eletron_s = counts

    if type_aper == 'ellipse':
        theta = radius
        A = (np.pi * radius[0] * radius[1])
        sigma = flux_eletron_s / A

    if type_aper == 'annulus':
        A = np.pi * (radius[1] * radius[2] - radius[0] * radius[3])
        sigma = flux_eletron_s / A

    n_pix = A / (pixel_scale ** 2)
    rms_calculate_aperture = 0.0062#sigma_clipped_stats(data, mask=np.isnan(data))[2]
    
    conversion_mjy = (header['PHOTFLAM'] * ((header['PHOTPLAM']) ** 2) * 3.33564095E+04) * 1e3
    #print(rms_calculate_aperture,rms_calculate_aperture*conversion_mjy)
    error_sigma = error_sigma_UV(n_pix, rms_calculate_aperture*header['EXPTIME'], flux_eletron_s*header['EXPTIME'], A)

    error_sigma = error_sigma/header['EXPTIME']

    sigma_mjy = sigma * conversion_mjy
    error_sigma_mjy = error_sigma * conversion_mjy

    if verbose:
        plt.figure()
        plt.imshow(data_weighted, origin='lower', vmin=-3 * sigma_clipped_stats(data, mask=np.isnan(data))[2],
                   vmax=15 * sigma_clipped_stats(data, mask=np.isnan(data))[2])

    return sigma_mjy, error_sigma_mjy



def photometry_ALMA(data, 
                    wcs, pixel_scale, 
                    radius, position, 
                    type_aper, theta_galaxy, 
                    verbose=False):
    """
    Calculate aperture photometry for ALMA data.
    
    Args:
    - data (2D array): ALMA data
    - wcs (astropy.wcs.WCS): WCS solution
    - pixel_scale (float): pixel scale in arcseconds
    - radius (list or tuple): aperture radius or annulus radii in arcseconds
    - position (list or tuple): position of the center of the aperture in degrees
    - type_aper (str): type of aperture ('ellipse' or 'annulus')
    - theta_galaxy (float): angle in degrees of the major axis of the ellipse
    - verbose (bool): if True, plot the weighted data in the aperture

    Returns:
    - float: aperture photometry value
    """
    # Convert center position from degrees to pixels
    position_pix = wcs.wcs_world2pix(position[0], position[1], 1)

    # Convert aperture radius or annulus radii from arcseconds to pixels
    radius = np.array(radius) 
    radius_pix = radius / pixel_scale

    # Define aperture or annulus
    if type_aper == 'ellipse':
        aper = EllipticalAperture((position_pix[0], position_pix[1]), radius_pix[0], radius_pix[1], theta=theta_galaxy)
    elif type_aper == 'annulus':
        aper = EllipticalAnnulus((position_pix[0], position_pix[1]), radius_pix[0], radius_pix[1], radius_pix[2], radius_pix[3], theta=theta_galaxy)

    # Create mask and apply it to the data
    masks = aper.to_mask(method='exact')
    data_weighted = masks.multiply(data)
    
    # Plot if verbose is True
    if verbose:
        plt.figure()
        plt.imshow(data_weighted, origin='lower')
    
    # Count number of pixels used in aperture
    number_total = np.count_nonzero(~np.isnan(masks.multiply(data)))

    # Perform aperture photometry
    phot_table = aperture_photometry(data, aper, wcs=wcs)
    return float(phot_table['aperture_sum'])


def sigma_CII(data, wcs, pixel_scale, radius, position, type_aper, bmaj, bmin, theta_galaxy, verbose=False):
    
    # 1. Get counts using ALMA photometry
    counts = photometry_ALMA(data, wcs, pixel_scale, radius, position, type_aper, theta_galaxy, verbose)  
    
    # 2. Convert counts to integrated flux in Jy
    # 2.1 Determine number of pixels per beam
    area_beam_arcsec = ((np.pi *(bmaj) * (bmin)) / (4 * np.log(2)))
    num_pix_beam = area_beam_arcsec / (pixel_scale)**2
    conversion_Jy_beam_to_Jy = 1 / num_pix_beam
    integrated_flux = counts * conversion_Jy_beam_to_Jy
    
    # 3. Calculate sigma in mJy
    # 3.1 Circular aperture
    if type_aper == 'ellipse':
        A = (np.pi * radius[0] * radius[1])
        sigma = integrated_flux / A
        
    # 3.2 Annular aperture
    if type_aper == 'annulus':
        A = np.pi * (radius[1] * radius[2] - radius[0] * radius[3])
        sigma = integrated_flux / A
    
    # 4. Convert sigma to mJy and calculate error in mJy
    sigma_mJy = sigma * 1e3
    error_sigma_mJy = ((np.nanstd(data) * np.sqrt(A / area_beam_arcsec))) * 1e3
    
    # 5. Return sigma and error in mJy
    return sigma_mJy, error_sigma_mJy






def plot_diagnosis(name_HST,data_HST,wcs_HST,rms_HST,ra_center_orig_pix, dec_center_orig_pix, imagesize_HST,pixscale_HST,eccen_galaxy,bpa_galaxy, convolved_HST,  rms_HST_convolved, ra_center_conv_pix, dec_center_conv_pix, b_array):
    fig_diag = plt.figure(figsize = (28,5))
    
    ax1 = fig_diag.add_subplot(141, projection = wcs_HST)
    box_ax1 = 1

    ax1.imshow(data_HST, vmin = -3*rms_HST, vmax = 15*rms_HST, cmap = 'viridis')
    ax1.contour(data_HST, levels = np.arange(3,30,8)*rms_HST, colors = 'black', linewidths = 1.)
    ax1.scatter(ra_center_orig_pix, dec_center_orig_pix, marker = '+', s = 100, color = 'brown' )

    ax1.set_xlim( imagesize_HST/2 - box_ax1/pixscale_HST, imagesize_HST/2 + box_ax1/pixscale_HST  )
    ax1.set_ylim( imagesize_HST/2 - box_ax1/pixscale_HST, imagesize_HST/2 + box_ax1/pixscale_HST  )
    
    ax1.set_xlabel('RA', fontsize = 15)
    ax1.set_ylabel('DEC', fontsize = 15)
    ax1.tick_params(fontsize = 12)


    ax2 = fig_diag.add_subplot(142, projection = wcs_HST)

    box_ax2 = 1

    ax2.imshow(data_HST, vmin = -3*rms_HST, vmax = 15*rms_HST, cmap = 'viridis')
    ax2.contour(data_HST, levels = np.arange(3,30,8)*rms_HST, colors = 'black', linewidths = 1., alpha = 0.5)
    for b in b_array/pixscale_HST:

        b = b*2

        e = Ellipse((ra_center_orig_pix, dec_center_orig_pix), b, b*np.sqrt(1 - eccen_galaxy**2), 
                    angle=bpa_galaxy , edgecolor = 'red', facecolor = 'None')
        ax2.add_patch(e)


    ax2.set_xlim( imagesize_HST/2 - box_ax2/pixscale_HST, imagesize_HST/2 + box_ax2/pixscale_HST  )
    ax2.set_ylim( imagesize_HST/2 - box_ax2/pixscale_HST, imagesize_HST/2 + box_ax2/pixscale_HST  )

    ax2.set_xlabel('RA', fontsize = 15)
    ax2.set_ylabel('DEC', fontsize = 15)
    ax2.tick_params(fontsize = 12)

    ax3 = fig_diag.add_subplot(143, projection = wcs_HST)

    box_ax3 = 1

    ax3.imshow(convolved_HST, vmin = -3*rms_HST_convolved, vmax = 15*rms_HST_convolved, cmap = 'viridis')
    ax3.contour(convolved_HST, levels = np.arange(3,30,8)*rms_HST_convolved, colors = 'black', linewidths = 1.)
    ax3.scatter(ra_center_conv_pix, dec_center_conv_pix, marker = '+', s = 100, color = 'brown' )

    ax3.set_xlim( imagesize_HST/2 - box_ax3/pixscale_HST, imagesize_HST/2 + box_ax3/pixscale_HST  )
    ax3.set_ylim( imagesize_HST/2 - box_ax3/pixscale_HST, imagesize_HST/2 + box_ax3/pixscale_HST  )

    ax3.set_xlabel('RA', fontsize = 15)
    ax3.set_ylabel('DEC', fontsize = 15)
    ax3.tick_params(fontsize = 12)

    ax4 = fig_diag.add_subplot(144, projection = wcs_HST)

    box_ax4 = 1
    ax4.imshow(convolved_HST, vmin = -3*rms_HST_convolved, vmax = 15*rms_HST_convolved, cmap = 'viridis')
    ax4.contour(convolved_HST, levels = np.arange(3,30,8)*rms_HST_convolved, colors = 'black', linewidths = 1., alpha = 0.5)
    ax4.scatter(ra_center_conv_pix, dec_center_conv_pix, marker = '+', s = 100, color = 'brown' )

    for b in b_array/pixscale_HST:
        b = b*2
        e = Ellipse((ra_center_conv_pix, dec_center_conv_pix), b, b*np.sqrt(1 - eccen_galaxy**2),
                    angle=bpa_galaxy , edgecolor = 'red', facecolor = 'None')
        ax4.add_patch(e)

    ax4.set_xlim( imagesize_HST/2 - box_ax4/pixscale_HST, imagesize_HST/2 + box_ax4/pixscale_HST  )
    ax4.set_ylim( imagesize_HST/2 - box_ax4/pixscale_HST, imagesize_HST/2 + box_ax4/pixscale_HST  )

    ax4.set_xlabel('RA', fontsize = 15)
    ax4.set_ylabel('DEC', fontsize = 15)
    ax4.tick_params(fontsize = 12)

    fig_diag.savefig('diagnosis_'+name_HST+'.png', dpi = 300)

    plt.close(fig_diag)





def run_ALMA(path_HST,path_ALMA,b_array,eccen_galaxy,bpa_galaxy, show_diagnosis = True,use_convolved_image= True):
    
    
    # Importing the HST data:
    data_HST, header_HST, wcs_HST, pixscale_HST, rms_HST,imagesize_HST,fhwm_HST,fhwm_angle_HST  = tim.import_fits_HST(path_HST[0])

    # Defining the HST sigma array
    sigma_output_array =  np.array([])
    sigma_output_array_error = np.array([])
    
    for file in range(len(path_ALMA)):
        
        # Importing ALMA data:
        data_ALMA, header_ALMA, wcs_ALMA, pixscale_ALMA, bmaj_ALMA, bmin_ALMA, bpa_ALMA, rms_ALMA, imagesize_ALMA = tim.import_fits_ALMA(path_ALMA[file])
        
        # Taking the name of the ALMA dataset
        name_ALMA = path_ALMA[file].split('ALL_')[-1].split('.fits')[0]#[-3] +'-'+path_ALMA[file].split('_')[-2]

        
        convolved_ALMA, convolved_HST = tim.convolution(data_HST_toconvolve = data_HST,
                                                    data_ALMA_toconvolve = data_ALMA, 
                                                    beam_ALMA_toconvolve = (bmaj_ALMA, bmin_ALMA, bpa_ALMA), 
                                                    fhwm_HST_toconvolve = (fhwm_HST,fhwm_HST, fhwm_angle_HST),
                                                    pixscale_HST_toconvolve = pixscale_HST, 
                                                    pixscale_ALMA_toconvolve = pixscale_ALMA)
        rms_ALMA_convolved = sigma_clipped_stats(convolved_ALMA, mask = np.isnan(convolved_ALMA))[2]




        
        ra_center_orig_pix, dec_center_orig_pix = tim.centroid(data = data_ALMA, 
                                                       x_init = int(imagesize_ALMA/2 + (0.1/pixscale_ALMA)), 
                                                       y_init= int(imagesize_ALMA/2 + (0./pixscale_ALMA)),
                                                       box_size = int(1.5/pixscale_ALMA), 
                                                       centroid_func = centroid_2dg, 
                                                       rms = rms_ALMA , 
                                                       name_cent_func = ' ', 
                                                       verbose_output = False)
        ra_center_orig, dec_center_orig = wcs_ALMA.all_pix2world(ra_center_orig_pix[0], dec_center_orig_pix[0],1)


        ra_center_conv_pix, dec_center_conv_pix = tim.centroid(data = convolved_ALMA, 
                                                       x_init = int(imagesize_ALMA/2 + (0.1/pixscale_ALMA)), 
                                                       y_init= int(imagesize_ALMA/2 + (-0.1/pixscale_ALMA)),
                                                       box_size = int(1.5/pixscale_ALMA), 
                                                       centroid_func = centroid_1dg, 
                                                       rms = rms_ALMA_convolved, 
                                                       name_cent_func = ' ', 
                                                       verbose_output = False)
        ra_center_conv, dec_center_conv = wcs_ALMA.all_pix2world(ra_center_conv_pix, dec_center_conv_pix,1) 

         
    
        #Diagnosis
        if show_diagnosis == True:
            plot_diagnosis(name_ALMA,data_ALMA,wcs_ALMA,rms_ALMA,ra_center_orig_pix, dec_center_orig_pix, imagesize_ALMA,pixscale_ALMA,eccen_galaxy,bpa_galaxy, convolved_ALMA,  rms_ALMA_convolved, ra_center_conv_pix, dec_center_conv_pix, b_array)
        
        

   
        if use_convolved_image == True:
            data = convolved_ALMA
            position = (ra_center_conv[0], dec_center_conv[0])
            bmaj_to_sigma =  (bmaj_ALMA**2 + fhwm_HST**2)**0.5 #0.3647 # 
            bmin_to_sigma =   (bmin_ALMA**2 + fhwm_HST**2)**0.5 # 0.3270 #
            
            

                
        else:
            data = data_ALMA
            position = (ra_center_orig, dec_center_orig)
            bmaj_to_sigma = bmaj_ALMA
            bmin_to_sigma = bmin_ALMA

            
        
        
        sigma_array = np.array([name_ALMA])
        sigma_array_error = np.array([name_ALMA])
        
        for rad in range(len(b_array)):
            if rad == 0:


                major_a = b_array[rad]
                minor_b = b_array[rad]*np.sqrt(1 - eccen_galaxy **2)

                array_radius = [major_a,minor_b]
                
                t_aper = 'ellipse'


            else:

                major_a_in = b_array[rad - 1]
                minor_b_in = b_array[rad - 1]*np.sqrt(1 - eccen_galaxy **2)
                major_a_out = b_array[rad]
                minor_b_out = b_array[rad]*np.sqrt(1 - eccen_galaxy **2)

                array_radius = [major_a_in ,major_a_out, minor_b_out,minor_b_in]

                t_aper = 'annulus'
            
            sig,  err = sigma_CII(data = data,
                               wcs = wcs_ALMA,
                               pixel_scale = pixscale_ALMA,
                               radius = array_radius,
                               position = position,
                               type_aper = t_aper,
                               bmaj = 0.36,#bmaj_ALMA_conv,
                               bmin = 0.32,#bmin_ALMA_conv,
                                  theta_galaxy = np.pi*(bpa_galaxy)/180)



            sigma_array = np.append(sigma_array ,sig)
            sigma_array_error = np.append(sigma_array_error ,err)
        
        if file == 0:
            sigma_output_array = sigma_array 
            sigma_output_array_error = sigma_array_error 
        else:
            sigma_output_array = np.vstack((sigma_output_array, sigma_array) ) #
            sigma_output_array_error =  np.vstack((sigma_output_array_error, sigma_array_error) ) # np.array([sigma_output_array_error, sigma_array_error] ,dtype=object)#
             
    return sigma_output_array, sigma_output_array_error

def run_ALMA_fromJWST(path_JWST,path_ALMA,b_array,eccen_galaxy,bpa_galaxy, show_diagnosis = True,use_convolved_image= True):
    
    
    # Importing the HST data:
    data_JWST, header_JWST, wcs_JWST, pixscale_JWST, rms_JWST,imagesize_JWST,fhwm_JWST,fhwm_angle_JWST  = tim.import_fits_JWST(path_JWST[0]) 

    # Defining the HST sigma array
    sigma_output_array =  np.array([])
    sigma_output_array_error = np.array([])
    
    for file in range(len(path_ALMA)):
        
        # Importing ALMA data:
        data_ALMA, header_ALMA, wcs_ALMA, pixscale_ALMA, bmaj_ALMA, bmin_ALMA, bpa_ALMA, rms_ALMA, imagesize_ALMA = tim.import_fits_ALMA(path_ALMA[file])
        
        # Taking the name of the ALMA dataset
        name_ALMA = path_ALMA[file].split('ALL_')[-1].split('.fits')[0]#[-3] +'-'+path_ALMA[file].split('_')[-2]

        
        convolved_ALMA, convolved_JWST = tim.convolution(data_HST_toconvolve = data_JWST,
                                                    data_ALMA_toconvolve = data_ALMA, 
                                                    beam_ALMA_toconvolve = (bmaj_ALMA, bmin_ALMA, bpa_ALMA), 
                                                    fhwm_HST_toconvolve = (fhwm_JWST,fhwm_JWST, fhwm_angle_JWST),
                                                    pixscale_HST_toconvolve = pixscale_JWST, 
                                                    pixscale_ALMA_toconvolve = pixscale_ALMA)
        rms_ALMA_convolved = sigma_clipped_stats(convolved_ALMA, mask = np.isnan(convolved_ALMA))[2]




        
        ra_center_orig_pix, dec_center_orig_pix = tim.centroid(data = data_ALMA, 
                                                       x_init = int(imagesize_ALMA/2 + (0.1/pixscale_ALMA)), 
                                                       y_init= int(imagesize_ALMA/2 + (0./pixscale_ALMA)),
                                                       box_size = int(1.5/pixscale_ALMA), 
                                                       centroid_func = centroid_2dg, 
                                                       rms = rms_ALMA , 
                                                       name_cent_func = ' ', 
                                                       verbose_output = False)
        ra_center_orig, dec_center_orig = wcs_ALMA.all_pix2world(ra_center_orig_pix[0], dec_center_orig_pix[0],1)


        ra_center_conv_pix, dec_center_conv_pix = tim.centroid(data = convolved_ALMA, 
                                                       x_init = int(imagesize_ALMA/2 + (0.1/pixscale_ALMA)), 
                                                       y_init= int(imagesize_ALMA/2 + (-0.1/pixscale_ALMA)),
                                                       box_size = int(1.5/pixscale_ALMA), 
                                                       centroid_func = centroid_1dg, 
                                                       rms = rms_ALMA_convolved, 
                                                       name_cent_func = ' ', 
                                                       verbose_output = False)
        ra_center_conv, dec_center_conv = wcs_ALMA.all_pix2world(ra_center_conv_pix, dec_center_conv_pix,1) 

         
    
        #Diagnosis
        if show_diagnosis == True:
            plot_diagnosis(name_ALMA,data_ALMA,wcs_ALMA,rms_ALMA,ra_center_orig_pix, dec_center_orig_pix, imagesize_ALMA,pixscale_ALMA,eccen_galaxy,bpa_galaxy, convolved_ALMA,  rms_ALMA_convolved, ra_center_conv_pix, dec_center_conv_pix, b_array)
        
        

   
        if use_convolved_image == True:
            data = convolved_ALMA
            position = (ra_center_conv[0], dec_center_conv[0])
            bmaj_to_sigma =  (bmaj_ALMA**2 + fhwm_JWST**2)**0.5 #0.3647 # 
            bmin_to_sigma =   (bmin_ALMA**2 + fhwm_JWST**2)**0.5 # 0.3270 #
            
            

                
        else:
            data = data_ALMA
            position = (ra_center_orig, dec_center_orig)
            bmaj_to_sigma = bmaj_ALMA
            bmin_to_sigma = bmin_ALMA

            
        
        
        sigma_array = np.array([name_ALMA])
        sigma_array_error = np.array([name_ALMA])
        
        for rad in range(len(b_array)):
            if rad == 0:


                major_a = b_array[rad]
                minor_b = b_array[rad]*np.sqrt(1 - eccen_galaxy **2)

                array_radius = [major_a,minor_b]
                
                t_aper = 'ellipse'


            else:

                major_a_in = b_array[rad - 1]
                minor_b_in = b_array[rad - 1]*np.sqrt(1 - eccen_galaxy **2)
                major_a_out = b_array[rad]
                minor_b_out = b_array[rad]*np.sqrt(1 - eccen_galaxy **2)

                array_radius = [major_a_in ,major_a_out, minor_b_out,minor_b_in]

                t_aper = 'annulus'
            
            sig,  err = sigma_CII(data = data,
                               wcs = wcs_ALMA,
                               pixel_scale = pixscale_ALMA,
                               radius = array_radius,
                               position = position,
                               type_aper = t_aper,
                               bmaj = 0.36,#bmaj_ALMA_conv,
                               bmin = 0.32,#bmin_ALMA_conv,
                                  theta_galaxy = np.pi*(bpa_galaxy)/180)



            sigma_array = np.append(sigma_array ,sig)
            sigma_array_error = np.append(sigma_array_error ,err)
        
        if file == 0:
            sigma_output_array = sigma_array 
            sigma_output_array_error = sigma_array_error 
        else:
            sigma_output_array = np.vstack((sigma_output_array, sigma_array) ) #
            sigma_output_array_error =  np.vstack((sigma_output_array_error, sigma_array_error) ) # np.array([sigma_output_array_error, sigma_array_error] ,dtype=object)#
             
    return sigma_output_array, sigma_output_array_error




def run_HST(path_HST,path_ALMA,b_array,eccen_galaxy,bpa_galaxy, show_diagnosis = True,use_convolved_image= True):#path_ALMA,path_HST, b_array, position, bmaj_galaxy, bmin_galaxy, bpa_galaxy, verbose = False):
    
    
    
    # Importing ALMA data:
    data_ALMA, header_ALMA, wcs_ALMA, pix_scale_ALMA, bmaj_ALMA, bmin_ALMA, bpa_ALMA, rms_ALMA, imagesize_ALMA = tim.import_fits_ALMA(path_ALMA[0])

    for file in range(len(path_HST)):#len(path_HST[3:])):

        # Taking HST data
        name_HST =  path_HST[file].split('_')[-3] 

        
        # Importing HST data:
        data_HST, header_HST, wcs_HST, pixscale_HST,  rms_HST, imagesize_HST, fhwm_HST, fhwm_angle_HST = tim.import_fits_HST(path_HST[file])


        
        # Convolve both images:
        convolved_ALMA, convolved_HST = tim.convolution(data_HST_toconvolve      =  data_HST,
                                                        data_ALMA_toconvolve     =  data_ALMA, 
                                                        beam_ALMA_toconvolve     =  (bmaj_ALMA, bmin_ALMA, bpa_ALMA), 
                                                        fhwm_HST_toconvolve      =  (fhwm_HST,fhwm_HST, fhwm_angle_HST),
                                                        pixscale_HST_toconvolve  =  pixscale_HST, 
                                                        pixscale_ALMA_toconvolve =  pix_scale_ALMA)
        
        rms_HST_convolved = sigma_clipped_stats(convolved_HST, mask = np.isnan(convolved_HST))[2]
        
        # Calculating centroids of original image:
        ra_center_orig_pix, dec_center_orig_pix = tim.centroid(data            = data_HST, 
                                                               x_init          = int(imagesize_HST/2 + (0.1/pixscale_HST)), 
                                                               y_init          = int(imagesize_HST/2 + (0./pixscale_HST)),
                                                               box_size        = int(1.5/pixscale_HST), 
                                                               centroid_func   = centroid_quadratic, 
                                                               rms             = rms_HST , 
                                                               name_cent_func  = ' ', 
                                                               verbose_output  = False)
        
        # Converting centroid in pixels to world coordinates
        ra_center_orig, dec_center_orig = wcs_HST.all_pix2world(ra_center_orig_pix[0], dec_center_orig_pix[0],1)

        # Calculating centroids of convolved image:
        ra_center_conv_pix, dec_center_conv_pix = tim.centroid(data            = convolved_HST, 
                                                               x_init          = int(imagesize_HST/2 + (0.1/pixscale_HST)), 
                                                               y_init          = int(imagesize_HST/2 + (0.1/pixscale_HST)),
                                                               box_size        = int(1.5/pixscale_HST), 
                                                               centroid_func   = centroid_quadratic, 
                                                               rms             = rms_HST_convolved, 
                                                               name_cent_func  = ' ', 
                                                               verbose_output  = False)
        #print('ra_center_conv_pix: ',ra_center_conv_pix, 'dec_center_conv_pix: ', dec_center_conv_pix)
        ra_center_conv, dec_center_conv = wcs_HST.all_pix2world(ra_center_conv_pix, dec_center_conv_pix,1) 
        #print('ra_center_conv: ',ra_center_conv, 'dec_center_conv: ', dec_center_conv)
        
        #Diagnosis
        if show_diagnosis == True:
            plot_diagnosis(name_HST,data_HST,wcs_HST,rms_HST,ra_center_orig_pix, dec_center_orig_pix, imagesize_HST,pixscale_HST,eccen_galaxy,bpa_galaxy, convolved_HST,  rms_HST_convolved, ra_center_conv_pix, dec_center_conv_pix, b_array)
        
        if use_convolved_image == True:
            data = convolved_HST
            position = (ra_center_conv[0], dec_center_conv[0])
            m = 'p'
            l = name_HST + ' - convolved'
        else:
            data = data_HST
            position = (ra_center_orig, dec_center_orig)
            m = 'X'
            l = name_HST
        
        sigma_array = np.array([name_HST])
        sigma_array_error = np.array([name_HST])
        
        for rad in range(len(b_array)):
            if rad == 0:


                major_a = b_array[rad]
                minor_b = b_array[rad]*np.sqrt(1 - eccen_galaxy **2)

                array_radius = [major_a,minor_b]

                t_aper = 'ellipse'

            else:

                major_a_in = b_array[rad - 1]
                minor_b_in = b_array[rad - 1]*np.sqrt(1 - eccen_galaxy **2)
                major_a_out = b_array[rad]
                minor_b_out = b_array[rad]*np.sqrt(1 - eccen_galaxy **2)

                array_radius = [major_a_in ,major_a_out, minor_b_out,minor_b_in]

                t_aper = 'annulus'

            sig, err = sigma_UV(data = data,
                                header = header_HST,
                              wcs = wcs_HST,
                              pixel_scale =  pixscale_HST,
                              radius = array_radius,
                              position =position,
                              type_aper = t_aper,
                              theta_galaxy =  np.pi*(bpa_galaxy)/180)
            
            
            sigma_array = np.append(sigma_array ,sig)
            sigma_array_error = np.append(sigma_array_error ,err)
        
        if file == 0:
            sigma_output_array = sigma_array 
            sigma_output_array_error = sigma_array_error 
        else:
            sigma_output_array = np.vstack((sigma_output_array, sigma_array) ) #
            sigma_output_array_error =  np.vstack((sigma_output_array_error, sigma_array_error) ) # np.array([sigma_output_array_error, sigma_array_error] ,dtype=object)#
             
    return sigma_output_array, sigma_output_array_error



    
def halo(path_HST,path_ALMA,b_array,eccen_galaxy,bpa_galaxy, show_diagnosis_HST = True, show_diagnosis_ALMA = True):
	
	
	sig_CII,sig_CII_error = run_ALMA(path_HST,path_ALMA,b_array,eccen_galaxy,bpa_galaxy, show_diagnosis = show_diagnosis_ALMA,use_convolved_image= True)
	
	#sig_dust,sig_dust_error = run_ALMA(path_HST,path_ALMA,b_array,eccen_galaxy,bpa_galaxy, show_diagnosis = True,use_convolved_image= True)
	
	sig_UV,sig_UV_error = run_HST(path_HST,path_ALMA,b_array,eccen_galaxy,bpa_galaxy, show_diagnosis = show_diagnosis_HST ,use_convolved_image= True)
	
	return sig_CII,sig_CII_error, sig_UV,sig_UV_error 
	
	
    
    
def calculate_background_threshold(Area, path_HST, path_ALMA):
    
    array_area = []
    rms_alt = []
    
    for a in Area:
    
        # Calculate radius in arcseconds
        radius_arcsec = (a/np.pi)**0.5


        # Import FITS data 
        data_HST, header_HST, wcs_HST, pixscale_HST, rms_HST, imagesize_HST, fhwm_HST, fhwm_angle_HST  = tim.import_fits_HST(path_HST) 
        
        data_ALMA, header_ALMA, wcs_ALMA, pixscale_ALMA, bmaj_ALMA, bmin_ALMA, bpa_ALMA, rms_ALMA, imagesize_ALMA = tim.import_fits_ALMA(path_ALMA)
        
        convolved_ALMA, convolved_HST = tim.convolution(data_HST_toconvolve = data_HST,
                                                    data_ALMA_toconvolve = data_ALMA,
                                                    beam_ALMA_toconvolve = (bmaj_ALMA, bmin_ALMA, bpa_ALMA), 
                                                    fhwm_HST_toconvolve = (fhwm_HST,fhwm_HST, fhwm_angle_HST),
                                                    pixscale_HST_toconvolve = pixscale_HST, 
                                                    pixscale_ALMA_toconvolve = pixscale_ALMA)
        
        
        data_HST = convolved_HST     


        # Set x and y limits for HST data
        rms_cut = sigma_clipped_stats(convolved_HST[300:700,300:700])[2]
        
        
        data_HST[data_HST > 6*rms_cut ] = np.nan


        # Define central coordinates in degrees
        ra_center_deg, dec_center_deg = 150.0393101, 2.3371604


        # Define arrays of RA and Dec coordinates for aperture photometry
        array_ra = np.arange(
                            ra_center_deg - 200*0.06/3600, 
                            ra_center_deg + 200*0.06/3600, 
                            2*radius_arcsec/3600
                            )

        array_dec = np.arange(
                              dec_center_deg - 200*0.06/3600, 
                              dec_center_deg + 200*0.06/3600, 
                              2*radius_arcsec/3600
                              )

        # Perform aperture photometry on each coordinate in the arrays
        array_phot = []

        for ra in array_ra:
            for dec in array_dec:

                # Convert RA and Dec to pixel coordinates
                ra_pix, dec_pix =  wcs_HST.wcs_world2pix(ra,dec,1)

                # Perform aperture photometry on the pixel coordinates
                phot = photometry_HST(data = data_HST, 
                                      wcs =  wcs_HST, 
                                      pixel_scale = pixscale_HST,
                                      radius = (radius_arcsec,radius_arcsec), 
                                      position = (ra,dec), 
                                      type_aper = 'ellipse', 
                                      theta_galaxy = 0, 
                                      verbose = False)
                if phot > -10:
                    # Append the result of aperture photometry to the array_phot
                    array_phot.append(phot)

        # Calculate the background threshold using the sigma-clipped stats of array_phot
        background = sigma_clipped_stats(array_phot)[2]

        # Calculate Area
        A = (np.pi*radius_arcsec**2)

        # Calculate conversion
        conversion_mjy = (header_HST['PHOTFLAM'] * ((header_HST['PHOTPLAM'])**2) *3.33564095E+04 ) *1e3
        
        
        
        
        array_area.append(background*1*conversion_mjy/A)
        
        
        num_pixel_area = A/pixscale_HST**2
        rms_alt.append(rms_cut/(num_pixel_area )**0.5)
        
    return array_area, rms_alt