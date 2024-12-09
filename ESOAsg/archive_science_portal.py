from astropy import coordinates
from astropy.coordinates import ICRS
import webbrowser

from ESOAsg import msgs
from ESOAsg.ancillary import cleaning_lists
from ESOAsg.core import asp_queries


def query_from_radec(positions, radius=None, instruments=None, data_types=None, open_link=False, show_link=False):
    r"""Query the ESO `ASP service <http://archive.eso.org/scienceportal/home>`_ given a position

     The `positions` value (or list) needs to be given as an `astropy.coordinates.SkyCoord` object. For further detail
     see here: `astropy coordinates <https://docs.astropy.org/en/stable/coordinates/>`_

    Args:
        positions (astropy.coordinates.SkyCoord): coordinates (or list of coordinates) of the sky you want to query
        radius (float): search radius in arcseconds
        instruments (list): list of `str` (or single `str`) containing the instruments used to limit the search
        data_types (list): list of `str` (or single `str`) containing the data types used to limit the search
        open_link (bool): open a link to the ASP page
        show_link (bool): show the link on the terminal

    Returns:
        None

    """
    # Check input position
    positions_list = cleaning_lists.from_element_to_list(positions, element_type=coordinates.SkyCoord)
    # Working on instruments
    instruments_list = cleaning_lists.from_element_to_list(instruments, element_type=str)
    # Working on data_types
    data_types_list = cleaning_lists.from_element_to_list(data_types, element_type=str)
    # Run one query per position
    for iii, position in enumerate(positions_list):
        position.transform_to(ICRS)
        url = '{0}{1}{2}{3}{4}'.format(asp_queries.base_url(),
                                       asp_queries.condition_position(position.ra.degree, position.dec.degree, radius,
                                                                      connector=None),
                                       asp_queries.condition_instruments(instruments_list, connector='&'),
                                       asp_queries.condition_data_types(data_types_list, connector='&'),
                                       asp_queries.sort_by('-obs_date', connector='&'))
        asp_queries.run_query(url, show_link=show_link, open_link=open_link)

    return


def query_from_polygons(polygons=None, instruments=None, data_types=None, open_link=False, show_link=False):
    r"""Query the ESO `ASP service <http://archive.eso.org/scienceportal/home>`_ given a polygon

    The `polygons` value (or list) needs to be given as a string defining the location in the sky of the polygon
    with RA, Dec, separated by commas and with the first RA, Dec pair that matches the last one (to close the
    polygon)

    Args:
        polygons (list): list of `str` (or single `str`) containing the coordinates of the polygon in the sky you want
            to query
        instruments (list): list of `str` (or single `str`) containing the instruments used to limit the search
        data_types (list): list of `str` (or single `str`) containing the data types used to limit the search
        open_link (bool): open a link to the ASP page
        show_link (bool): show the link on the terminal

    Returns:
        None

    """
    # Check input position
    polygons_list = cleaning_lists.from_element_to_list(polygons, element_type=str)
    # Working on instruments
    instruments_list = cleaning_lists.from_element_to_list(instruments, element_type=str)
    # Working on data_types
    data_types_list = cleaning_lists.from_element_to_list(data_types, element_type=str)
    # Run one query per position
    for iii, polygon in enumerate(polygons_list):
        url = '{0}{1}{2}{3}{4}'.format(asp_queries.base_url(),
                                       asp_queries.condition_polygon(polygon, connector=None),
                                       asp_queries.condition_instruments(instruments_list, connector='&'),
                                       asp_queries.condition_data_types(data_types_list, connector='&'),
                                       asp_queries.sort_by('-obs_date', connector='&'))
        asp_queries.run_query(url, show_link=show_link, open_link=open_link)

    return


def open_url(url, show_link=True, open_link=True):
    r"""Open some url in the browser

    Args:
        url (str): url to open
        open_link (bool): open a link to page
        show_link (bool): show the link on the terminal

    Returns:
        None

    """
    if show_link:
        msgs.info('The link is:\n {}\n'.format(url))
    if open_link:
        webbrowser.open(url)
    return

def check_preview(dp_ids, show_link=False, open_link=True): 
    """Check the preview of the data product
        Note that this will download .pdf file to the browser download directory"""
    msgs.info('This will download a .pdf file from the browser directory')

    # Check if list
    if isinstance(dp_ids, str):
        dp_ids = [dp_ids]
    dp_ids = cleaning_lists.from_element_to_list(dp_ids, element_type=str)
    
    # Loop over the dp_ids
    for dp_id in dp_ids:
        url = f'https://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?{dp_id}&eso_download=preview'
        open_url(url, show_link, open_link)


def check_header(dp_ids, show_link=False, open_link=True): 
    """Check the header of the data product"""
    msgs.info('This will open a new tab with the header information in browser')

    # Check if list
    if isinstance(dp_ids, str):
        dp_ids = [dp_ids]
    dp_ids = cleaning_lists.from_element_to_list(dp_ids, element_type=str)

    # Loop over the dp_ids
    for dp_id in dp_ids:
        url = f'https://archive.eso.org/hdr?DpId={dp_id}'
        open_url(url, show_link, open_link)


def check_details(dp_ids, show_link=False, open_link=True): 
    """Check the details of the data product"""
    msgs.info('This will open a new tab with the details in browser')

    # Check if list
    if isinstance(dp_ids, str):
        dp_ids = [dp_ids]
    dp_ids = cleaning_lists.from_element_to_list(dp_ids, element_type=str)

    # Loop over the dp_ids
    for dp_id in dp_ids:
        url = f'https://archive.eso.org/dataset/{dp_id}'
        open_url(url, show_link, open_link)


def check_hips(dp_ids, show_link=False, open_link=True):
    """Check the HIPS of the data product"""
    msgs.info('This will open a new tab with the HIPS information in browser')

    # Check if list
    if isinstance(dp_ids, str):
        dp_ids = [dp_ids]
    dp_ids = cleaning_lists.from_element_to_list(dp_ids, element_type=str)

    # Loop over the dp_ids
    for dp_id in dp_ids:
        url = f'https://archive.eso.org/previews/v1/files/{dp_id}/hips'
        open_url(url, show_link, open_link)


def check_release_info(dp_ids, show_link=False, open_link=True): 
    """Check the release info of the data product"""
    msgs.info('This will open a new tab with the release description infomation in browser')

    # Check if list
    if isinstance(dp_ids, str):
        dp_ids = [dp_ids]
    dp_ids = cleaning_lists.from_element_to_list(dp_ids, element_type=str)

    # Loop over the dp_ids
    for dp_id in dp_ids:
        url = f'https://archive.eso.org/datalink/links?ID=ivo://eso.org/ID?{dp_id}&eso_download=data_documentation'
        open_url(url, show_link, open_link)


def check_metadata(dp_ids, show_link=False, open_link=True):
    """Check the HTML metadata of the data product"""
    msgs.info('This will open a new tab with the metadata information in browser')

    # Check if list
    if isinstance(dp_ids, str):
        dp_ids = [dp_ids]
    dp_ids = cleaning_lists.from_element_to_list(dp_ids, element_type=str)

    # Loop over the dp_ids
    for dp_id in dp_ids:
        url = f'https://archive.eso.org/wdb/wdb/adp/phase3_main/query?dp_id={dp_id}'
        open_url(url, show_link, open_link)