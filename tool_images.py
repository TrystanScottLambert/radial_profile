import warnings
from astropy.io import fits
warnings.filterwarnings('ignore', category=UserWarning, append=True)
warnings.simplefilter('ignore', UserWarning)
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astropy.
 import Gaussian2DKernel
from astropy.convolution import convolve, convolve_fft
import matplotlib.pyplot as plt
import numpy as np
from photutils.centroids import centroid_sources
from matplotlib.patches import Ellipse
#import pyregion

def import_fits_ALMA(path_file, cube = False):
    
    '''
    Import main parameters of ALMA cubes/images:
    
    Input: directory path of a single fits file. 
    
    Output parameters:
    
    1) data
    2) header 
    3) wcs 
    4) pixelscale 
    5) bmaj - arcsec
    6) bmin - arcsec
    7) bpa - degree
    8) rms 
    9) imagesize  
    
    
    '''
    f = fits.open(path_file)
    
    data = f[0].data[0][0]
    header = f[0].header
    wcs = WCS(header, naxis =2)

    pixelscale = 3600 * header['CDELT2']
    imagesize = f[0].data[0][0].shape[0]

    if cube == False:
        data = f[0].data[0][0]

        bmaj = 3600 * header['BMAJ']
        bmin = 3600 * header['BMIN']
        bpa =  header['BPA']
        rms = sigma_clipped_stats(data, mask = np.isnan(data))[2]

    else:

        bmaj = np.array(f[1].data.field('BMAJ'))
        bmin = np.array(f[1].data.field('BMIN'))
        bpa = np.array(f[1].data.field('BPA'))

        rms = []
        data = f[0].data[0]
        for i in range(len(data)):
            rms.append(np.nanstd(data[i])) 


    return data, header, wcs, pixelscale,bmaj,bmin,bpa,rms, imagesize

# def plot_channels(path_file, channel_range,vmin_rms,vmax_rms,cmap,rms_levels,box_size, moment0_path, cont_moment0):
    
#     data_ALMA, header_ALMA, wcs_ALMA, pixscale_ALMA,bmaj_ALMA,bmin_ALMA,bpa_ALMA,rms_ALMA, imagesize_ALMA = import_fits_ALMA(path_file, cube = True)
    
    
#     nrow = 0
#     ncol = 4
#     n_chan = len(channel_range)
    
#     vel_array = (1 + channel_range - header_ALMA['CRPIX3'])*header_ALMA['CDELT3'] + header_ALMA['CRVAL3']
    
#     vel_array = 2.99792458e5 * (header_ALMA['RESTFRQ'] - vel_array )/header_ALMA['RESTFRQ']
 
    
#     if n_chan < 4:
#         nrow = 1
#         ncol = n_chan
#         plt.figure(figsize= (4*n_chan,5))
#     elif n_chan%4 == 0:
#         nrow = int(n_chan/4)
#         plt.figure(figsize= (26,5*nrow))
#     else:
#         nrow = int(n_chan/4) + 1
#         plt.figure(figsize= (26,5*nrow))
    

#     for ind in range(n_chan):
        
        
#         set_plot_param()
#         ax = plt.subplot(nrow,4,ind + 1, projection = wcs_ALMA)
    
#         ax.imshow(data_ALMA[channel_range[ind]], vmin = vmin_rms*rms_ALMA[ind], vmax = vmax_rms*rms_ALMA[ind], cmap = cmap)
                   
#         cs = ax.contour(data_ALMA[channel_range[ind]], levels = rms_levels*rms_ALMA[ind], colors = 'black' )

#         ax.text( imagesize_ALMA/2 - (box_size*0.75)/pixscale_ALMA, imagesize_ALMA/2 + (box_size*0.8)/pixscale_ALMA,    str(round(vel_array[ind])) + ' km/s - '+str(channel_range[ind]), color = 'white', fontsize = 30,va='center' )
        
#         if moment0_path != False:
#             data_mom0, header_mom0, wcs_mom0, pixscale_mom0,bmaj_mom0,bmin_mom0,bpa_mom0,rms_mom0, imagesize_mom0 = import_fits_ALMA(moment0_path)

#             ax.contour(data_mom0, levels = cont_moment0*rms_mom0, colors = 'red' ,transform= ax.get_transform(wcs_mom0))
        
        
#         try:

  
#             region = 'Natural_10_Channel_'+str(int(channel_range[ind]))+'.reg'
#             r = pyregion.open(region).as_imagecoord(header_ALMA)
#             p,n = r.get_mpl_patcimport pyregionhes_texts()

#             ax.add_patch(p[0])


#             print(m,n)
#         except:
#             pass

#         ax.set_xlabel('Right Ascension', fontsize = 25)
#         ax.set_ylabel('Declination', fontsize = 25)
#         ax.tick_params(axis='both',  labelsize=25)
#         ax.set_xlim(imagesize_ALMA/2 - box_size/pixscale_ALMA, imagesize_ALMA/2 + box_size/pixscale_ALMA)
#         ax.set_ylim(imagesize_ALMA/2 - box_size/pixscale_ALMA, imagesize_ALMA/2 + box_size/pixscale_ALMA)

#     plt.show()
    

# #         beam = Ellipse((data_ALMA.shape[0]/2 - (box_size - bmaj_ALMA)/pixscale_ALMA, data_ALMA.shape[0]/2 - (box_size-bmaj_ALMA )/pixscale_ALMA), bmaj_ALMA/pixscale_ALMA,bmin_ALMA/pixscale_ALMA, bpa_ALMA + 90, 
# #                        facecolor = 'white', edgecolor = 'black', zorder = 2, lw = 1.5)

# #         ax.add_patch(beam)

# #         ax.text( data_ALMA.shape[0]/2 - (box_size - 2*bmaj_ALMA)/pixscale_ALMA , data_ALMA.shape[0]/2 - (box_size-bmaj_ALMA )/pixscale_ALMA,     str(np.round(bmaj_ALMA,2)) + '"x'+str(np.round(bmaj_ALMA,2))+'"', color = 'white', fontsize = 25,va='center' )    
    

def import_fits_HST(path_file):
    
    '''
    Import main parameters of HST cubes/images:
    
    Input: directory path of a single fits file. 
    
    Output parameters:
    
    1) data
    2) header 
    3) wcs 
    4) pixelscale 
    5) rms
    6) imagesize
    7) fwhm - arcsec
    8) fwhm angle - degree   
    
    '''
    f = fits.open(path_file)
    
    data = f[0].data
    header = f[0].header
    header['CTYPE1'] = header['CTYPE1']+str('-SIP')
    header['CTYPE2'] = header['CTYPE2']+str('-SIP')
    wcs = WCS(header, naxis =2)
    wcs.sip = None
    
    pixelscale = 3600 * header['CD2_2']

    channelsize = data.shape[0]    
    rms = sigma_clipped_stats(data, mask = np.isnan(data))[2]
    
    fhwm = 0.18
    fhwm_angle = 0.
    
    return data, header, wcs, pixelscale, rms,channelsize, fhwm, fhwm_angle


def import_fits_JWST(path_file):
    
    '''
    Import main parameters of JWST images:
    
    Input: directory path of a single fits file. 
    
    Output parameters:
    
    1) data
    2) header 
    3) wcs 
    4) pixelscale 
    5) rms
    6) imagesize
    7) fwhm - arcsec
    8) fwhm angle - degree   
    
    '''
    f = fits.open(path_file)

    data = f[0].data
    header = f[0].header
    header['CTYPE1'] = header['CTYPE1']+str('-SIP')
    header['CTYPE2'] = header['CTYPE2']+str('-SIP')
    wcs = WCS(header, naxis =2)
    wcs.sip = None

    pixelscale = 3600 * header['CDELT2']

    channelsize = data.shape[0]    
    rms = sigma_clipped_stats(data, mask = np.isnan(data))[2]

    fhwm = 0.269
    fhwm_angle = 0.

    data = data - sigma_clipped_stats(data, mask = np.isnan(data))[1]

    rms = sigma_clipped_stats(data, mask = np.isnan(data))[2]

    return data, header, wcs, pixelscale, rms,channelsize, fhwm, fhwm_angle

def convolution(data_HST_toconvolve,
                data_ALMA_toconvolve,
                beam_ALMA_toconvolve,
                pixscale_ALMA_toconvolve,
                fhwm_HST_toconvolve,
                pixscale_HST_toconvolve):

    '''
    
    data_HST_toconvolve - data HST to convolve
    
    data_ALMA_toconvolve - data ALMA to convolve
    
    beam_ALMA_toconvolve - (bmaj_ALMA, bmin_ALMA, bpa not corrected) arcsec,arcsec,degre
    
    fhwm_HST_toconvolve - (fwhm-maj_HST, fwhm-maj_HST, bpa not corrected) arcsec,arcsec,degre
    
    output ALMA convolved, HST convolved
    
    '''

    #Sigma ALMA with HST pixel scale
    sigma_maj_ALMA_in_HST = (beam_ALMA_toconvolve[0]/2.355)/(pixscale_HST_toconvolve)
    sigma_min_ALMA_in_HST = (beam_ALMA_toconvolve[1]/2.355)/(pixscale_HST_toconvolve)

    #Sigma HST with ALMA pixel scale
    sigma_maj_HST_in_ALMA = (fhwm_HST_toconvolve[0]/2.355)/(pixscale_ALMA_toconvolve)
    sigma_min_HST_in_ALMA = (fhwm_HST_toconvolve[1]/2.355)/(pixscale_ALMA_toconvolve)

    beam_ALMA_in_HST = Gaussian2DKernel(x_stddev = sigma_maj_ALMA_in_HST,
                                        y_stddev=  sigma_min_ALMA_in_HST,
                                        theta = ((beam_ALMA_toconvolve[2] +90) * np.pi)/180,
                                        x_size = data_HST_toconvolve.shape[0],
                                        y_size = data_HST_toconvolve.shape[0])


    psf_HST_in_ALMA = Gaussian2DKernel(x_stddev = sigma_maj_HST_in_ALMA,
                                       y_stddev=sigma_min_HST_in_ALMA,
                                       theta = ((fhwm_HST_toconvolve[2] +90) * np.pi)/180,
                                       x_size = data_ALMA_toconvolve.shape[0] ,
                                       y_size = data_ALMA_toconvolve.shape[0])


    conv_HST = convolve_fft(data_HST_toconvolve, beam_ALMA_in_HST)

    conv_ALMA = convolve_fft(data_ALMA_toconvolve, psf_HST_in_ALMA)

    return conv_ALMA, conv_HST

def centroid(data,
                 x_init, y_init,
                 box_size, centroid_func,
                 rms , name_cent_func,
                 verbose_output = False):


    x, y = centroid_sources(data, x_init,y_init, box_size = box_size,
                        centroid_func=centroid_func)

    if verbose_output == True:
        fig_centroid = plt.figure()
        axc = fig_centroid.add_subplot(111)
        axc.imshow(data, origin='lower', vmin = -3*rms, vmax = 20*rms)
        axc.set_xlim(x_init - box_size/2, x_init + box_size/2)
        axc.set_ylim(y_init - box_size/2, y_init + box_size/2)
        axc.scatter(x, y, marker='+', s=200, label = name_cent_func)
    return x,y


def eccentricity(bmaj_galaxy,bmin_galaxy):
    
    return np.sqrt(1 - (bmin_galaxy/bmaj_galaxy)**2)


def set_plot_param():
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size = 20)
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['xtick.major.size'] = 12
    plt.rcParams['xtick.minor.size'] = 8
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.minor.size'] = 7
    plt.rcParams['ytick.minor.size'] = 7
    plt.rcParams['xtick.minor.visible'] = True
    plt.rcParams['ytick.minor.visible'] = True
    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['ytick.right'] = True 
    plt.rcParams['xtick.major.width']  = 1
    plt.rcParams['ytick.major.width']  = 1

def plot_maps(pathfile_ALMA, box_size, vmin_rms, vmax_rms, cmap, rms_levels, title_image, saveas):
    
    data_ALMA, header_ALMA, wcs_ALMA, pixscale_ALMA, bmaj_ALMA, bmin_ALMA, bpa_ALMA, rms_ALMA, imagesize_ALMA = import_fits_ALMA(pathfile_ALMA)
    
    set_plot_param()
    
    plt.figure(figsize= (10,10))
    
    ax = plt.subplot(1,1,1, projection = wcs_ALMA)
    
    ax.imshow(data_ALMA, vmin = vmin_rms*rms_ALMA, vmax = vmax_rms*rms_ALMA, cmap = cmap)
    
    cs = ax.contour(data_ALMA, levels = rms_levels*rms_ALMA, colors = 'black' )
    #cnames = {rms_levels[i]*rms_ALMA: str(rms_levels[i]) for i in range(len(rms_levels))}
        
    
    #ax.clabel(cs, levels =  rms_levels, inline=True, fontsize=8)
    ax.text( data_ALMA.shape[0]/2 - (box_size*0.8)/pixscale_ALMA, data_ALMA.shape[0]/2 + (box_size*0.8)/pixscale_ALMA,    title_image, color = 'white', fontsize = 35,va='center' )
    
    beam = Ellipse((data_ALMA.shape[0]/2 - (box_size - bmaj_ALMA)/pixscale_ALMA, data_ALMA.shape[0]/2 - (box_size-bmaj_ALMA )/pixscale_ALMA), bmaj_ALMA/pixscale_ALMA,bmin_ALMA/pixscale_ALMA, bpa_ALMA + 90, 
                   facecolor = 'white', edgecolor = 'black', zorder = 2, lw = 1.5)
    
    ax.add_patch(beam)
    
    ax.text( data_ALMA.shape[0]/2 - (box_size - 2*bmaj_ALMA)/pixscale_ALMA , data_ALMA.shape[0]/2 - (box_size-bmaj_ALMA )/pixscale_ALMA,     str(np.round(bmaj_ALMA,2)) + '"x'+str(np.round(bmaj_ALMA,2))+'"', color = 'white', fontsize = 25,va='center' )
    
    ax.set_xlabel('Right Ascension', fontsize = 25)
    ax.set_ylabel('Declination', fontsize = 25)
    ax.tick_params(axis='both',  labelsize=25)
    ax.set_xlim(data_ALMA.shape[0]/2 - box_size/pixscale_ALMA, data_ALMA.shape[0]/2 + box_size/pixscale_ALMA)
    ax.set_ylim(data_ALMA.shape[0]/2 - box_size/pixscale_ALMA, data_ALMA.shape[0]/2 + box_size/pixscale_ALMA)

    if saveas != None:
        plt.savefig(saveas, dpi = 300)
    
