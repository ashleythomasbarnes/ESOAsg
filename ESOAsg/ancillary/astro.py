"""
Module to perform some quick astronomical calculations

The GW part is taken from a python script written by Martino Romaniello
http://www.eso.org/~mromanie/martino/Home.html
"""

import numpy as np
import requests
import json

from astropy import units as u
from astropy import constants
from astropy import wcs

import importlib

from os import path
from os import remove
from astroquery.mast import Observations
import healpy 

from ESOAsg import msgs
from ESOAsg.ancillary import cleaning_lists
from ESOAsg.ancillary import checks
from ESOAsg.ancillary import polygons  # ToDo to be removed
from ESOAsg.core import fitsfiles

import matplotlib
from matplotlib import pyplot as plt

# TODO - Also needs fixing... 
# STARTING_MATPLOTLIB_BACKEND = matplotlib.rcParams["backend"]
from ligo.skymap.tool.ligo_skymap_contour import main as ligo_skymap_contour
import ligo.skymap.plot 
# matplotlib.rcParams["backend"] = STARTING_MATPLOTLIB_BACKEND


def download_kepler_data(target_name):
    r"""Run a query to MAST to obtain the kepler data of an object

    The data is saved in the directories: `./mastDownload/Kepler/<target_name>_*/*.fits`

    Args:
        target_name (str): name of the target: e.g., 'kplr011446443'

    Returns:
        bool: it returns `True` if data are retrieved and `False` otherwise.

    """
    kepler_obs = Observations.query_criteria(target_name=target_name, obs_collection='Kepler')
    if len(kepler_obs) == 0:
        msgs.warning('No Kepler data for target {}'.format(target_name))
        return False
    for idx in np.arange(0, len(kepler_obs)):
        kepler_prods = Observations.get_product_list(kepler_obs[idx])
        kepler_fits_prods = Observations.filter_products(kepler_prods, extension='fits', mrp_only=False)
        Observations.download_products(kepler_fits_prods, mrp_only=False, cache=False)
    return True


def download_tess_data(target_name):
    r"""Run a query to MAST to obtain the TESS data of an object.

    Data is saved in the directories: `./mastDownload/TESS/<target_name>_*/*.fits`

    Args:
        target_name (str): name of the target: e.g., '231663901'

    Returns:
        bool: it returns `True` if data are retrieved and `False` otherwise.

    """
    tess_obs = Observations.query_criteria(target_name=target_name, obs_collection='TESS')
    if len(tess_obs) == 0:
        msgs.warning('No TESS data for target {}'.format(target_name))
        return False
    for idx in np.arange(0, len(tess_obs)):
        tess_prods = Observations.get_product_list(tess_obs[idx])
        tess_fits_prods = Observations.filter_products(tess_prods, extension='fits', mrp_only=False)
        Observations.download_products(tess_fits_prods, mrp_only=False, cache=False)
    return True


def download_gw_bayestar(superevent_name, file_name='bayestar.fits.gz', output_dir='./'):
    r"""Download the bayestar.fits.gz of a GW superevent

    This is simply checking the existence of: `https://gracedb.ligo.org/superevents/<superevent_name>/files/`
    and downloading the file: `https://gracedb.ligo.org/apiweb/superevents/<superevent_name>/iles/bayestar.fits.gz`

    The downloaded file will be renamed as `superevent_name`+`file_name`

    Args:
        superevent_name (str): name of the superevent: e.g., 'S191205ah'
        file_name (str): name of the file to download. Default is `bayestar.fits.gz`

    Returns:
        bool: it returns `True` if data are retrieved and `False` otherwise.

    """

    # Some checks
    assert isinstance(superevent_name, str), '{} is not a valid string'.format(superevent_name)
    assert isinstance(file_name, str), '{} is not a valid string'.format(file_name)
    
    # Checks if superevent exists
    gw_files_url = 'https://gracedb.ligo.org/superevents/' + superevent_name + '/files/'
    checks.connection_to_website(gw_files_url)
    
    # download the and save the superevent file
    gw_bayestar_url = 'https://gracedb.ligo.org/apiweb/superevents/' + superevent_name + '/files/' + file_name
    r = requests.get(gw_bayestar_url, allow_redirects=True)

    if r.status_code == 404:
        msgs.warning('Failed to access to: {}'.format(gw_bayestar_url))

    open(output_dir+superevent_name+'_bayestar.fits.gz', 'wb').write(r.content)

    # check that the file actually arrived on disk
    if path.isfile(superevent_name + '_bayestar.fits.gz'):
        _test_header = fitsfiles.header_from_fits_file(superevent_name + '_' + file_name)
        if len(_test_header) > 0:
            msgs.info('File {}_{} successfully downloaded'.format(superevent_name, file_name))
        else:
            msgs.warning('Not a valid fits file in: {}'.format(gw_bayestar_url))
    else:
        msgs.warning('Already downloaded: {}'.format(gw_bayestar_url))


def contours_from_gw_bayestar(file_name, credible_level=50.):
    r"""Given a bayestar.fits.gz HEALPix maps of a GW superevent it extracts contours at the given `credible_level`

    This is a wrapper around `ligo_skymap_contour`. The different contours are then counterclockwise oriented to be
    compatible with TAP and ASP queries.

    Args:
        file_name (str): Name of HEALPix maps
        credible_level (float): Probability level at which contours are returned

    Returns:
        list: [RA,Dec] list defining the contours location. These are counter-clockwise oriented as seen on the sky from
           inside the sphere. The length of contours represent the number of not connected regions identified.
    """

    # Some checks
    assert isinstance(file_name, (str)), '{} is not a valid string'.format(file_name)
    if not checks.fits_file_is_valid(file_name):
        msgs.error('{} not a valid fits file'.format(file_name))
    assert isinstance(credible_level, (int, float)), '`credible_level` is not a float or an int'

    # Create temporary file where to store output from ligo_skymap_contour
    contour_tmp_file = file_name + '.tmp.json'
    if path.isfile(contour_tmp_file):
        remove(contour_tmp_file)

    # Use ligo_skymap_contour to compute contours
    ligo_skymap_contour(args=[file_name, '--output', contour_tmp_file, '--contour', str(credible_level)])

    # Parse the content of the resulting json file into a dict
    with open(contour_tmp_file, 'r') as json_file:
        json_data = json_file.read()
    contours_dict = json.loads(json_data)
    json_file.close()
    
    # cleaning up
    remove(contour_tmp_file)

    # Parse the resulting json to obtain a list of coordinates defining the vertices of the contours of each peak
    # contours is a list, each element of which contains the vertices of one contour encircling the desired significance
    contours = contours_dict['features'][0]['geometry']['coordinates']
    # Make sure that the orientation of the polygons on the sky is the correct one for the TAP queries,
    # i.e. counter-clockwise as seen on the sky from inside the sphere.
    for iii, contour in enumerate(contours):
        contours[iii] = _ensure_orientation(contour)

    # Quick summary:
    if len(contours) > 0:
        msgs.info('Extracted the contours for {} regions at {} credible level'.format(len(contours), credible_level))
    else:
        msgs.info('No contours extracted at {} credible level'.format(credible_level))
        contours = None

    return contours


def _ensure_orientation(contour):
    r"""Set orientation of a polyong on the sky to be counterclockwise

    Take the input polygon and turn it into counterclockwise orientation as seen on the sky from inside the
    unitary sphere. This is done  through a `spherical_geometry.polygon.SphericalPolygon` object, which
    upon creation enforces the correct orientation independently of the one of the input (they call it clockwise,
    because it's see from outside).

    Args:
        contour (`list`):
            list of the vertices of a polygon

    Returns:
        contour_oriented (`list`)
            oriented list of the vertices of a polygon
    """

    # create the spherical polygon, which by construction has the right orientation on the sky
    contour_array = np.asarray(contour, dtype=np.float32)
    spherical_polygon = polygons.SphericalPolygon.from_lonlat(contour_array[:, 0], contour_array[:, 1], degrees=True)

    # reformat the correctly-oriented polygon to the same format as the input
    contour_oriented = []
    for radec in spherical_polygon.to_lonlat():
        for ra, dec in zip(radec[0], radec[1]):
            contour_oriented.append([ra, dec])

    return contour_oriented


def _array_split(array_in, thr):
    r"""Split an input array if two consecutive elements in x differ
    more than thr. The purpose is to deal with RAs folding at 360 degrees.
    """
    xx_in, yy_in = array_in[:, 0], array_in[:, 1]

    # Rearrange the input arrays to so that its first element is the most eastern one.
    ii = np.argmax(xx_in)
    xx_in = np.append(xx_in[ii:], xx_in[0:ii])
    yy_in = np.append(yy_in[ii:], yy_in[0:ii])

    dxx = np.diff(xx_in)
    split_indices = np.where(np.abs(dxx) > thr)
    return np.split(xx_in, split_indices[0] + 1), np.split(yy_in, split_indices[0] + 1)


def show_contours_from_gw_bayestar(file_name, contours=None, cmap='cylon', contours_color='black',
                                   credible_level=None, save_figure=None, return_fig=False):
    """
    Display a sky credibility map from a HEALPix map file and overlay contours indicating credible regions.

    Args:
        file_name (str): Path to the FITS file containing the HEALPix map to display.
        contours (list, optional): Precomputed contours to overlay. Defaults to None, in which case they are computed.
        cmap (str, optional): Colormap to use for the HEALPix map. Defaults to 'cylon'.
        contours_color (str, optional): Color of the contours. Defaults to 'black'.
        credible_level (int or None, optional): Credible level percentage for the contours if they are computed.
        save_figure (str or None, optional): File path to save the figure. Defaults to None (no saving).
        return_fig (bool, optional): Whether to return the figure object. Defaults to False.

    Returns:
        matplotlib.figure.Figure or None: The figure object if `return_fig` is True, otherwise None.
    """

    # Ensure the file_name argument is a valid string.
    assert isinstance(file_name, (str)), f'{file_name} is not a valid string'
    
    # Check if the FITS file is valid using a helper function.
    if not checks.fits_file_is_valid(file_name):
        msgs.error(f'{file_name} is not a valid FITS file')

    # Handle contours: if none are provided, generate them based on the credible level.
    if contours is None:
        msgs.warning("No contours defined, showing `credible_level=50` contours")
        credible_level = credible_level or 50  # Default to 50% credible level if not specified
        plot_contours = contours_from_gw_bayestar(file_name, credible_level=credible_level)
    else:
        plot_contours = contours

    # Read the HEALPix map data and header from the FITS file.
    map_data, map_header = healpy.read_map(file_name, h=True, verbose=False, dtype=None)
    map_header = dict(map_header)  # Convert header to a dictionary for easier access.

    # Extract object name and observation date/time from the FITS header.
    object_name = map_header['OBJECT']
    date = map_header['DATE'].split('T')[0]
    time = map_header['DATE'].split('T')[1]

    # Initialize a plot with a Mollweide projection.
    # Other nice examples given here - https://lscsoft.docs.ligo.org/ligo.skymap/plot/allsky.html
    fig = plt.figure(figsize=(10., 7.))
    ax = fig.add_subplot(111, projection='astro hours mollweide')

    # Display the HEALPix map on the plot using the specified colormap.
    ax.imshow_hpx(map_data, cmap=cmap, visible=True, zorder=1)
    
    # Add text annotations for the GW event name and credible level (if specified).
    ax.text(0.0, 0.99, f'GW event - {object_name}', transform=ax.transAxes)
    if credible_level is not None:
        ax.text(0.0, 0.95, f'Credible level: {credible_level}%', transform=ax.transAxes, fontsize=8)

    # Prepare for drawing contours: convert world coordinates to pixel coordinates.
    w = wcs.WCS(ax.header)  # Create a WCS object from the plot's header.

    # Determine the step size to split contours into segments for smoother plotting.
    split_step = (10 * (ax.get_xlim()[1] - ax.get_xlim()[0]) / 360.)  # 10-degree step in RA.

    # Loop through each contour and draw it on the plot.
    for contour_world in plot_contours:
        contour_pix = w.wcs_world2pix(contour_world, 0)  # Convert contour coordinates to pixels.
        x_contours, y_contours = _array_split(contour_pix, split_step)  # Split into manageable segments.
        for x_contour, y_contour in zip(x_contours, y_contours):
            ax.plot(x_contour, y_contour, linewidth=0.75, color=contours_color, zorder=5)

    # Add a grid to the plot for reference.
    ax.grid(True, ls=':', alpha=0.2, color='black')

    # Save the figure to the specified file if `save_figure` is provided.
    if save_figure is not None: 
        msgs.info(f'Saving figure to {save_figure}')
        # print(f'Saving figure to {save_figure}')
        fig.savefig(save_figure, bbox_inches='tight')

    # Return the figure object if `return_fig` is True.
    if return_fig: 
        return fig
    

def mag2f_nu(abmags):
    r"""Convert list of AB magnitudes in a list of spectral flux densities per unit frequency in jansky

    The input can be either a list of `astropy.units.Quantity` or a list of `float`. In the latter
    case it is assumed that the entry is in mag

    Args:
        abmags (any): input list of magnitudes in the AB system to be converted in flux density in jansky

    Returns:
        list: list of spectral flux densities per unit frequency. Each element of the list is an
            `astropy.units.Quantity` with jansky as unit

    """
    abmags_list = cleaning_lists.from_element_to_list_of_quantities(abmags, unit=u.mag)
    f_nus_list = []
    for abmag in abmags_list:
        f_nus_list.append(np.power(10., (8.90 - abmag.value) / 2.5) * u.jansky)
    return f_nus_list


def f_nu2mag(f_nus):
    r"""Convert list of spectral flux densities per unit frequency in a list of AB magnitudes

    The input can be either a list of `astropy.units.Quantity` or a list of `float`. In the latter
    case it is assumed that the entry is in jansky

    Args:
        f_nus (any): input list of spectral flux densities per unit frequency to be converted to AB magnitudes. If
            `f_nus` does not have units it is assumed to be in jansky

    Returns:
        list: list of AB magnitudes. Each element of the list is an `astropy.units.Quantity` with
            mag as unit

    """
    f_nus_list = cleaning_lists.from_element_to_list_of_quantities(f_nus, unit=u.jansky)
    abmag_list = []
    for f_nu in f_nus_list:
        abmag_list.append(((-2.5 * np.log10(f_nu.value)) + 8.90) * u.mag)
    return abmag_list


def f_nu2f_lambda(f_nus, wavelengths):
    r"""Convert list of spectral flux densities per unit frequency into spectral flux densities per unit  wavelength

    Args:
        f_nus (any): input list of spectral flux densities per unit frequency to be converted to AB magnitudes. If
            `f_nus` does not have units it is assumed to be in jansky
        wavelengths (any): specific wavelengths at which perform the conversion. If the input does not have units, Ang.
            is assumed. It must be either a single value or a list with the same length of `f_nus`
    Returns:
        list: list of spectral flux densities per unit wavelength

    """
    f_nus_list = cleaning_lists.from_element_to_list_of_quantities(f_nus, unit=u.jansky)
    wavelengths_list = cleaning_lists.from_element_to_list_of_quantities(wavelengths, unit=u.AA)
    f_lambdas_list = []
    if len(wavelengths_list) == 1:
        for f_nu in f_nus_list:
            f_lambdas_list.append(f_nu.to(u.erg*u.centimeter**-2.*u.second**-1.*u.AA**-1.,
                                          equivalencies=u.spectral_density(wavelengths_list[0])))
    elif len(wavelengths_list) == len(f_nus_list):
        for f_nu, wavelength in zip(f_nus_list, wavelengths_list):
            f_lambdas_list.append(f_nu.to(u.erg*u.centimeter**-2.*u.second**-1.*u.AA**-1.,
                                          equivalencies=u.spectral_density(wavelength)))
    else:
        msgs.error('Length of `wavelengths` not compatible with the one of `f_nus`')
        return
    return f_lambdas_list


def f_lambda2f_nu(f_lambdas, wavelengths):
    r"""Convert list of spectral flux densities per unit wavelength into spectral flux densities per unit frequency

    Args:
        f_lambdas (any): input list of spectral flux densities per unit wavelength to be converted to AB magnitudes. If
            `f_lambdas` does not have units it is assumed to be in erg/s/cm**2/Ang.
        wavelengths (any): specific wavelengths at which perform the conversion. If the input does not have units, Ang.
            is assumed. It must be either a single value or a list with the same length of `f_lambdas`

    Returns:
        list: list of spectral flux densities per unit frequency

    """
    f_lambdas_list = cleaning_lists.from_element_to_list_of_quantities(f_lambdas,
                                                                       unit=u.erg*u.centimeter**-2.*u.second**-1.*u.AA**-1.)
    wavelengths_list = cleaning_lists.from_element_to_list_of_quantities(wavelengths, unit=u.AA)
    f_nus_list = []
    if len(wavelengths_list) == 1:
        for f_lambda in f_lambdas_list:
            f_nus_list.append(f_lambda.to(u.jansky,
                                          equivalencies=u.spectral_density(wavelengths_list[0])))
    elif len(wavelengths_list) == len(f_lambdas_list):
        for f_lambda, wavelength in zip(f_lambdas_list, wavelengths_list):
            f_nus_list.append(f_lambda.to(u.jansky,
                                          equivalencies=u.spectral_density(wavelength)))
    else:
        msgs.error('Length of `wavelengths` not compatible with the one of `f_lambdas`')
        return
    return f_nus_list


def mag2f_lambda(abmags, wavelengths):
    r"""Convert list of AB magnitudes in a list of spectral flux densities per unit frequency in jansky

    The input can be either a list of `astropy.units.Quantity` or a list of `float`. In the latter
    case it is assumed that the entry is in mag

    Args:
        abmags (any): input list of magnitudes in the AB system to be converted in flux density in jansky
        wavelengths (any): specific wavelengths at which perform the conversion. If the input does not have units, Ang.
            is assumed. It must be either a single value or a list with the same length of `abmags`

    Returns:
        list: list of spectral flux densities per unit wavelength. Each element of the list is an
            `astropy.units.Quantity` with erg/s/cm**2/Ang as unit

    """
    abmags_list = cleaning_lists.from_element_to_list_of_quantities(abmags, unit=u.mag)
    f_nus_list = []
    for abmag in abmags_list:
        f_nus_list.append(np.power(10., (8.90 - abmag.value) / 2.5) * u.jansky)
    return f_nu2f_lambda(f_nus_list, wavelengths)


def f_lambda2mag(f_lambdas, wavelengths):
    r"""Convert list of spectral flux densities per unit wavelength in a list of AB magnitudes

    The input can be either a list of `astropy.units.Quantity` or a list of `float`. In the latter
    case it is assumed that the entry is in jansky

    Args:
        f_lambdas (any): input list of spectral flux densities per unit frequency to be converted to AB magnitudes. If
            `f_lambdas` does not have units it is assumed to be in erg/s/cm**2/Ang
        wavelengths (any): specific wavelengths at which perform the conversion. If the input does not have units, Ang.
            is assumed. It must be either a single value or a list with the same length of `f_lambdas`

    Returns:
        list: list of AB magnitudes. Each element of the list is an `astropy.units.Quantity` with
            mag as unit

    """
    f_lambdas_list = cleaning_lists.from_element_to_list_of_quantities(f_lambdas,
                                                                       unit=u.erg*u.centimeter**-2.*u.second**-1.*u.AA**-1.)
    wavelengths_list = cleaning_lists.from_element_to_list_of_quantities(wavelengths, unit=u.AA)
    f_nus_list = f_lambda2f_nu(f_lambdas_list, wavelengths_list)
    abmag_list = []
    for f_nu in f_nus_list:
        abmag_list.append(((-2.5 * np.log10(f_nu.value)) + 8.90) * u.mag)
    return abmag_list


def mag2sb(abmags, areas):
    r"""Convert list of AB magnitudes into surface brightness given an area

    .. math::

        SB = m + 2.5 log10(area/arcsec^2)

    Args:
        abmags (any): input list of magnitudes in the AB system to be converted in surface brightness in mag/arcsec**2
        areas (any): area over which calculate the surface brightness. If the input does not have units, arcsec**2
            is assumed. It must be either a single value or a list with the same length of `abmags`

    Returns:
        list: list of surface brightnesses. Each element of the list is an `astropy.units.Quantity` with
            mag/arcsec**2 as unit
    """
    abmags_list = cleaning_lists.from_element_to_list_of_quantities(abmags, unit=u.mag)
    areas_list = cleaning_lists.from_element_to_list_of_quantities(areas, unit=u.arcsec**2.)
    surface_brightnesses_list = []
    if len(areas_list) == 1:
        for abmag in abmags_list:
            surface_brightnesses_list.append((abmag + (2.5 * np.log10(areas_list[0]/u.arcsec**2.)))/u.arcsec**2.)
    elif len(areas_list) == len(abmags_list):
        for abmag, area in zip(abmags_list, areas_list):
            surface_brightnesses_list.append((abmag + (2.5 * np.log10(area/u.arcsec**2.)))/u.arcsec**2.)
    else:
        msgs.error('Length of `areas` not compatible with the one of `abmags`')
        return
    return surface_brightnesses_list


def sb2mag(surface_brightnesses, areas):
    r"""Convert list of surface brightness into AB magnitudes given an area

    .. math::

        m = SB - 2.5 log10(area/arcsec^2)

    Args:
        surface_brightnesses (any): input list of surface brightnesses to be converted in magnitudes in mag. If the
            input does not have units, mag/arcsec**2 is assumed.
        areas (any): area over which calculate the surface brightness. If the input does not have units, arcsec**2
            is assumed. It must be either a single value or a list with the same length of `sbs`

    Returns:
        list: list of magnitudes. Each element of the list is an `astropy.units.Quantity` with mag as unit

    """
    surface_brightnesses_list = cleaning_lists.from_element_to_list_of_quantities(surface_brightnesses,
                                                                                  unit=u.mag/u.arcsec**2)
    areas_list = cleaning_lists.from_element_to_list_of_quantities(areas, unit=u.arcsec**2.)
    abmag_list = []
    if len(areas_list) == 1:
        for surface_brightness in surface_brightnesses_list:
            abmag_list.append((surface_brightness*u.arcsec**2.) - (2.5 * np.log10(areas_list[0]/u.arcsec**2.) * u.mag))
    elif len(areas_list) == len(surface_brightnesses_list):
        for surface_brightness, area in zip(surface_brightnesses_list, areas_list):
            abmag_list.append((surface_brightness*u.arcsec**2.) - (2.5 * np.log10(area/u.arcsec**2.) * u.mag))
    else:
        msgs.error('Length of `areas` not compatible with the one of `surface_brightnesses`')
        return
    return abmag_list


# ToDo -> ABMAGLIM

def abmaglim(rms, seeing_fwhm, exptime=None, zero_point=None, sigma=5.):
    r"""Calculate the N-sigma magnituide limit.

    The code behaves in different ways depending if zero_point is set to None or not.
    If zero_points is None:
        rms needs to have units and exptime is not considered.
        abmaglim = -2.5 * Log10( sigma * rms * PI*(.5*seeing_fwhm)**2. / (3631.*u.jansky) )
    else:
        abmaglim is calculated as:
        abamglim = -2.5 * Log10( sigma * rms * PI*(.5*seeing_fwhm)**2. / exptime ) + zero_point

    Args:
        rms (`float`):
           Calculated rms on the image
        exptime (`float`):
            Exposure time in seconds
        zero_point (`float`):
            AB zero point of the image. If not set, the code will take the RMS with the associated
            astropy.units.
        seeing_fwhm (`float`):
            FWHM of the seeing in pixels
        sigma (`float`):
           Number of sigma to consider a detection significant

    Returns:
        abmaglim (`float`):
            AB magnitude limit of for an image
    """
    msgs.work('Calculating {}-sigma AB mag limit'.format(sigma))

    if zero_point is None:
        if exptime is not None:
            msgs.error('`exptime` needs to be `None` if `zero_point` is `None`')
        rms_in_jansky = rms.to(u.jansky, equivalencies=u.spectral())
        abmaglim = -2.5 * np.log10(sigma * rms_in_jansky * np.pi * np.power(seeing_fwhm / 2., 2.) / (3631. * u.jansky))
    else:
        abmaglim = -2.5 * np.log10(sigma * rms * np.pi * np.power(seeing_fwhm / 2., 2.) / exptime) + zero_point
    return abmaglim

