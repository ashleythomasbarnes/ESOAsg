from astropy import coordinates
from astropy.coordinates import ICRS

import numpy as np # Library for numerical operations
import urllib # Library for opening URLs
import requests # Library for making HTTP requests
import os # Library for interacting with the operating system
import cgi  # Library for parsing headers

from ESOAsg import msgs
from ESOAsg import default
from ESOAsg.core import tap_queries
from ESOAsg.queries import query_observations
from ESOAsg.ancillary import checks
from ESOAsg.ancillary import cleaning_lists


def query_from_radec(positions=None, radius=None, instruments=None, data_types=None, top=None, columns=None, verbose=False,
                     maxrec=None, em_min=None, em_max=None, enclosed=True, snr=None, em_res_power=None, order_by=None):
    r"""Query the ESO archive for data at a given position in RA and Dec

    The `positions` value (or list) needs to be given as an
    `astropy.coordinates.SkyCoord <https://docs.astropy.org/en/stable/coordinates/>`_ object.

    The output is in an (list of) `astropy.table` with columns defined in: `core.tap_queries.COLUMNS_FROM_OBSCORE`
    It is possible to change the columns to query by setting the value of `columns`


    .. note::
        In case you are querying radius=`None` is set, the query will performed with:
        `INTERSECT(POINT('',RA,Dec), s_region)`
        instead of:
        `INTERSECT(s_region,CIRCLE('',RA,Dec,radius/3600.))`.
        See here for further examples: `tap obs examples <http://archive.eso.org/tap_obs/examples>`_


    Args:
        positions (astropy.coordinates.SkyCoord): coordinates (or list of coordinates) of the sky you want to query
        radius (float, optional): search radius in arcseconds
        instruments (list): list of `str` (or single `str`) containing the instruments used to limit the search
        data_types (list): list of `str` (or single `str`) containing the data types used to limit the search
        columns (list): list of `str` (or single `str`) containing the columns to be queried
        verbose (bool): if set to `True` additional info will be displayed
        maxrec (int, optional): define the maximum number of entries that a single query can return. If it is `None` the
            value is set by the limit of the service.
        em_min (float, optional): minimum value of the wavelength to be queried
        em_max (float, optional): maximum value of the wavelength to be queried

    Returns:
        any: results from the queries

    """

    # Check inputs:
    # Working on positions
    positions_list = cleaning_lists.from_element_to_list(positions, element_type=coordinates.SkyCoord)
    # Working on radius
    if radius is not None:
        if isinstance(radius, int):
            radius = float(radius)
        else:
            assert isinstance(radius, float), r'Input radius is not a number'
    # Working on instruments
    instruments_list = cleaning_lists.from_element_to_list(instruments, element_type=str)
    # Working on data_types
    data_types_list = cleaning_lists.from_element_to_list(data_types, element_type=str)
    # Working on columns
    columns_list = _is_column_list_in_obscore(columns)

    if verbose:
        how_many_positions = len(positions_list)
        if how_many_positions > 1:
            msgs.work('Exploring ESO archive around {} locations in the sky'.format(how_many_positions))
        else:
            msgs.work('Exploring ESO archive around the input location in the sky')

    # Running over all positions
    results_from_query = []
    for idx, position in enumerate(positions_list):
        position.transform_to(ICRS)
        ra, dec = np.float_(position.ra.degree), np.float_(position.dec.degree)
        msgs.work('Running query {} to the ESO archive (out of {} total)'.format(idx + 1, len(positions_list)))

        # Define query
        query = "{0}{1}{2}{3}{4}{5}{6}{7}".format(tap_queries.create_query_obscore_base(top, columns_list),
                                                        tap_queries.condition_intersects_ra_dec(ra, dec, radius=radius),
                                                        tap_queries.condition_instruments_like(instruments_list),
                                                        tap_queries.condition_data_types_like(data_types_list),
                                                        tap_queries.condition_em_like(em_min, em_max, enclosed),
                                                        tap_queries.condition_snr_like(snr), 
                                                        tap_queries.condition_em_res_power_like(em_res_power),
                                                        tap_queries.condition_order_by_like(order_by))
        
        # instantiate ESOCatalogues
        query_for_observations = query_observations.ESOObservations(query=query, type_of_query='sync', maxrec=maxrec)

        # running query and append results to the list
        if verbose:
            query_for_observations.print_query()

        # Obtaining query results
        query_for_observations.run_query(to_string=True)
        result_from_query = query_for_observations.get_result_from_query()
        if len(result_from_query) < 1:
            msgs.warning('No data has been retrieved')
        else:
            msgs.info('A total of {} entries has been retrieved'.format(len(result_from_query)))
            if verbose:
                msgs.info('For the following instrument:')
                for inst_name in np.unique(result_from_query['instrument_name'].data):
                    msgs.info(' - {}'.format(inst_name))

        results_from_query.append(result_from_query)

    # Returning results
    return _return_results_from_query(results_from_query)


def query_from_polygons(polygons, instruments=None, data_types=None, verbose=False, columns=None, maxrec=None,
                        em_min=None, em_max=None):
    r"""Query the ESO archive for data at a area in the sky defined by a polygon

    The `polygons` value (or list) needs to be given as a string defining the location in the sky of the polygon
    with RA, Dec, separated by commas and with the first RA, Dec pair that matches the last one (to close the
    polygon)

    The output is in an (list of) `astropy.table` with columns defined in: `core.tap_queries.COLUMNS_FROM_OBSCORE`
    It is possible to change the columns to query by setting the value of `columns`


    Args:
        polygons (list): ist of `str` (or single `str`) containing the coordinates of the polygon in the sky you want
            to query
        instruments (list): list of `str` (or single `str`) containing the instruments used to limit the search
        data_types (list): list of `str` (or single `str`) containing the data types used to limit the search
        columns (list): list of `str` (or single `str`) containing the columns to be queried
        verbose (bool): if set to `True` additional info will be displayed
        maxrec (int, optional): define the maximum number of entries that a single query can return. If it is `None` the
            value is set by the limit of the service.
        em_min (float, optional): minimum value of the wavelength to be queried
        em_max (float, optional): maximum value of the wavelength to be queried

    Returns:
        any: results from the queries

    """
    # Check inputs:
    # Working on polygons
    polygons_list = cleaning_lists.from_element_to_list(polygons, element_type=str)
    # Working on instruments
    instruments_list = cleaning_lists.from_element_to_list(instruments, element_type=str)
    # Working on data_types
    data_types_list = cleaning_lists.from_element_to_list(data_types, element_type=str)
    # Working on columns
    columns_list = _is_column_list_in_obscore(columns)

    if verbose:
        how_many_polygons = len(polygons_list)
        if how_many_polygons > 1:
            msgs.work('Exploring ESO archive around {} locations in the sky'.format(how_many_polygons))
        else:
            msgs.work('Exploring ESO archive around the input location in the sky')

    # Running over all positions
    results_from_query = []
    for idx, polygon in enumerate(polygons_list):
        msgs.work('Running query {} to the ESO archive (out of {} total)'.format(idx + 1, len(polygons_list)))
        # Define query
        query = "{0}{1}{2}{3}{4}".format(tap_queries.create_query_obscore_base(columns_list),
                                      tap_queries.condition_intersects_polygon(polygon),
                                      tap_queries.condition_instruments_like(instruments_list),
                                      tap_queries.condition_data_types_like(data_types_list),
                                      tap_queries.condition_em_min_like(em_min, em_max))
        # instantiate ESOCatalogues
        query_for_observations = query_observations.ESOObservations(query=query, type_of_query='sync', maxrec=maxrec)
        # running query and append results to the list
        if verbose:
            query_for_observations.print_query()
        # Obtaining query results
        query_for_observations.run_query(to_string=True)
        result_from_query = query_for_observations.get_result_from_query()
        if len(result_from_query) < 1:
            msgs.warning('No data has been retrieved')
        else:
            msgs.info('A total of {} entries has been retrieved (with maxrec={})'.format(len(result_from_query),
                                                                                         maxrec))
            msgs.info('For the following instrument:')
            for inst_name in np.unique(result_from_query['instrument_name'].data):
               msgs.info(' - {}'.format(inst_name))

        results_from_query.append(result_from_query)

    # Returning results
    return _return_results_from_query(results_from_query)


def getDispositionFilename( response ): #get the filename from the header
    """Get the filename from the Content-Disposition in the response's http header"""
    contentdisposition = response.headers.get('Content-Disposition')
    if contentdisposition == None:
        return None
    value, params = cgi.parse_header(contentdisposition)
    filename = params["filename"]
    return filename


def writeFile( response, folder_path): #write content of file on disk
    """Write on disk the retrieved file specifying a folder path"""
    if response.status_code == 200:
        # The ESO filename can be found in the response header
        filename = getDispositionFilename( response )
        # Let's write on disk the downloaded FITS spectrum using the ESO filename:
        full_path = os.path.join(folder_path, filename)
        with open(full_path, 'wb') as f:
            f.write(response.content)
        return filename 


def download(dp_ids, output_dir='./', min_disk_space=float(default.get_value('min_disk_space')), 
             base_url='https://dataportal.eso.org/dataPortal/', 
             verbose=False):
    r"""Given a filename in the ADP format, the code download the file from the dataportal.eso.org

    Args:
        dp_ids (any): list data product ID (or single product ID) to be downloaded
        min_disk_space (float): the file will be downloaded only if there is this amount of space (in Gb) free on the
            disk
        base_url (str): base url to be used to download the file
        verbose (bool): if set to `True` additional info will be displayed

    Returns:
        None

    """
    # Check for disk space
    checks.check_disk_space(min_disk_space=min_disk_space)

    # Cleaning list
    dp_ids_list = cleaning_lists.from_element_to_list(cleaning_lists.from_bytes_to_string(dp_ids), element_type=str)

    for dp_id in dp_ids_list:

        msgs.work('Retrieving file {}.fits'.format(dp_id))
   
        download_url = f"{base_url}file/{dp_id}"
        if verbose:
                print(download_url)
        response = requests.get(download_url)

        if response.status_code != 200:
            msgs.info(f'No file downloaded - Response code: {response.status_code}')
        else: 
            output_filename = writeFile(response, output_dir) # Write the file to disk and get the filename
            msgs.info(f'File {output_filename} downloaded')


def download_cutout(dp_ids, output_dir='./', min_disk_space=float(default.get_value('min_disk_space')),
                    base_url = 'https://dataportal.eso.org/dataPortal/', 
                    em_min_cut=None, em_max_cut=None, 
                    verbose=False):
    r"""Given a filename in the ADP format, the code download the file from the dataportal.eso.org
    A portion of the spectrum can be cropped by specifying the minimum and maximum wavelength
    TODO: add the possibility to crop... 

    Args:
        dp_ids (any): list data product ID (or single product ID) to be downloaded
        min_disk_space (float): the file will be downloaded only if there is this amount of space (in Gb) free on the
            disk
        base_url (str): base url to be used to download the file
        em_min_cut (float, optional): minimum value of the wavelength to be cropped (in units of m)
        em_max_cut (float, optional): maximum value of the wavelength to be cropped (in units of m)
        verbose (bool): if set to `True` additional info will be displayed

    Returns:
        None

    """
    # Check for disk space
    checks.check_disk_space(min_disk_space=min_disk_space)

    # Cleaning list
    dp_ids_list = cleaning_lists.from_element_to_list(cleaning_lists.from_bytes_to_string(dp_ids), element_type=str)

    for dp_id in dp_ids_list:

        msgs.work('Retrieving file {}.fits'.format(dp_id))
  
        if em_min_cut is not None or em_max_cut is not None:

            em_min_int = int(em_min_cut*1e9) #convert to nm
            em_max_int = int(em_max_cut*1e9) #convert to nm
            cut_append = f"&PREFIX={em_min_int}-{em_max_int}&BAND={em_min_int}e-9+{em_max_int}e-9"
            download_url = f"{base_url}soda/sync?ID={dp_id}{cut_append}"
            if verbose:
                print(download_url)
            response = requests.get(download_url)

            if response.status_code != 200:
                msgs.info(f'No file downloaded - Response code: {response.status_code}')
            else: 
                output_filename = writeFile(response, output_dir) # Write the file to disk and get the filename
                msgs.info(f'File {output_filename} downloaded')

def columns_info(verbose=False):
    r"""Load a query that get names (and corresponding ucd) of the columns present in `ivoa.ObsCore`

    Args:
        verbose (bool): if set to `True` additional info will be displayed

    Returns:
        astropy.table: table of all columns present in a table/collection. Information are stored in `table_name`,
            `column_name`, `ucd`, `datatype`, `description`, and `unit`

    """
    # instantiate ESOCatalogues
    query_all_columns_info = query_observations.ESOObservations(query=tap_queries.create_query_obscore_all_columns())
    # Print query
    if verbose:
        query_all_columns_info.print_query()
    # Obtaining query results
    query_all_columns_info.run_query(to_string=True)
    all_columns_table = query_all_columns_info.get_result_from_query()
    return all_columns_table


def _is_column_list_in_obscore(columns):
    r"""Check if a given list of columns is present in `ivoa.ObsCore`

    Args:
        columns (any): list of string containing the column_name (or the single `str`) to be tested

    Returns:
        list: same of `columns` but columns not present at in the collections/tables are removed

    """
    assert columns is None or isinstance(columns, (str, list)), r'`columns` must be `None` or a `str` or a `list`'
    columns_list = cleaning_lists.from_element_to_list(columns, element_type=str)
    if columns is not None:
        # test if it is a valid column
        clean_columns = []
        for column in columns_list:
            if _is_column_in_obscore(column):
                clean_columns.append(column)
    else:
        clean_columns = None
    return clean_columns


def _is_column_in_obscore(column_name):
    r"""Check if a given column is present in `ivoa.ObsCore`

    Args:
        column_name (str): column to be tested

    Returns:
        bool: `True` if the column is present in `ivoa.ObsCore`. `False` and warning raised otherwise

    """
    is_in_obscore = True
    table_all_columns = columns_info(verbose=False)
    all_column_list = table_all_columns['column_name'].data.data.tolist()
    if column_name not in all_column_list:
        msgs.warning('Column: {} not recognized. Possible values are:\n{}'.format(column_name, all_column_list))
        is_in_obscore = False
    return is_in_obscore


def _return_results_from_query(results_from_query):
    r"""Return the result from a query

    This is ancillary to return the result from a query in different format depending on the input. There are three
    options:

    - if len(`results_from_query`) == 0: it returns `None`
    - if len(`results_from_query`) == 1: it returns a single astropy table
    - if len(`results_from_query`) > 1: it returns a list of astropy tables

    Args:
        results_from_query (list): input list of tables to be transformed (if necessary)

    Returns:
        any: `None`, a single astropy table, or a list of astropy tables depending on the input

    """
    if len(results_from_query) == 0:
        return None
    elif len(results_from_query) == 1:
        return results_from_query[0]
    else:
        return results_from_query
