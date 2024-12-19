r"""Module to create and run TAP queries

The Table Access Protocol (TAP) is a web-service protocol that gives access to collections of tabular data referred to
collectively as a `tableset`. TAP services accept queries posed against the `tableset` available via the service and
return the query response as another table, in accord with the relational model. Queries to the ESO TAP services are
submitted using the Astronomical Data Query Language ADQL [O08]_.

Currently, ESOAsg offers two TAP services:

* `eso_tap_obs`: to query both the database tables describing the observed raw and reduced data obtained at the La
  Silla Paranal Observatory, and the database table containing the ESO ambient conditions and meteorological
  measurements (seeing, isoplanatic angle, precipitable water, turbulence profiles, etc.)

* `eso_tap_cat`: to query the scientific catalogues provided by the principal investigators of ESO observing programmes

More information on the ESO TAP service are at `this web-page
<http://archive.eso.org/cms/eso-data/programmatic-access.html>`_ and some examples are given in `these notebooks
<http://archive.eso.org/programmatic/HOWTO/>`_

.. [O08] Oritz et al., (2008) `IVOA Astronomical Data Query Language <http://www.ivoa.net/documents/latest/ADQL.html>`_

"""

from pyvo import dal
from pyvo.dal import DALQueryError
from pyvo.dal import DALFormatError

from ESOAsg import default
from ESOAsg import msgs

import numpy as np

# currently supported tap services
TAP_SERVICES = ['eso_tap_cat', 'eso_tap_obs']
# type of queries currently allowed
TAP_QUERY_TYPES = ['sync', 'async']
# list of columns that will be queried in
COLUMNS_FROM_OBSCORE = ['target_name', 'dp_id', 's_ra', 's_dec', 't_exptime', 'em_min', 'em_max', 'dataproduct_type',
                        'instrument_name', 'obstech', 'abmaglim', 'proposal_id', 'obs_collection']


def define_tap_service(which_tap_service):
    r"""Load a Table Access Protocol (TAP) service from defaults

    Currently the supported `TAP services <http://archive.eso.org/programmatic/#TAP>`_ are:
    * `eso_tap_cat`: TAP service for scientific catalogues generated by ESO observing teams
    * `eso_tap_obs`: TAP service for raw, reduced, and ambient data

    See `pyvo docs <https://pyvo.readthedocs.io/en/latest/api/pyvo.dal.TAPService.html#>`_ for further details

    Args:
        which_tap_service (str): Select the `TAP services <http://archive.eso.org/programmatic/#TAP>`_ to be queried

    Returns:
        pyvo.dal.tap.TAPService: TAP service used for the queries

    """
    if which_tap_service not in TAP_SERVICES:
        msgs.error('{} not a valid entry for TAP services. Possibilities are: {}'.format(which_tap_service,
                                                                                         TAP_SERVICES))
    tap_service = dal.tap.TAPService(default.get_value(which_tap_service))
    return tap_service


def print_query(query):
    r"""Print the query on the terminal

    In case the `query` is empty, a warning is raised

    Args:
        query (str): String containing the query

    Returns:
          None

    """
    if query is None:
        msgs.warning('The query is empty')
    else:
        msgs.info('The query is:')
        msgs.info('{}'.format(query))
    return


def which_service(tap_service):
    r"""Print a summary description of the TAP service used

    Args:
        tap_service (pyvo.dal.tap.TAPService): TAP service used for the queries

    Returns:
        None

    """
    msgs.info('The TAP service used is:')
    tap_service.describe()
    return


def run_query(tap_service, query, type_of_query, maxrec=default.get_value('maxrec')):
    r"""Run query to TAP service and return result as an `astropy.Table`

    If the job requires to much time to run, the code will move to an asynchronous query.

    Args:
        tap_service (pyvo.dal.tap.TAPService): TAP service that will be used for the query
        query (str): query to be run
        type_of_query (str): type of query to be run
        maxrec (int): define the maximum number of entries that a single query can return. Default is set
            by default.get_value('maxrec')

    Returns:
        astropy.table: result from the query to the TAP service

    """
    if type_of_query not in TAP_QUERY_TYPES:
        msgs.error('{} not a valid entry for the type of TAP query. Possibilities are: {}'.format(type_of_query,
                                                                                                  TAP_QUERY_TYPES))
    # Obtaining query results and convert it to an astropy table
    if query is not None:
        if type_of_query == 'sync':
            result_from_query = run_query_sync(tap_service, query, maxrec=maxrec)
        else:
            result_from_query = run_query_async(tap_service, query, maxrec=maxrec)
    else:
        msgs.warning('Empty query provided')
        result_from_query = None
    return result_from_query


def run_query_sync(tap_service, query, maxrec=default.get_value('maxrec')):
    r"""Run a synchronous query to TAP service and return result as an `astropy.Table`

    If the synchronous query fails, the code automatically tries to run the same query asynchronously

    Args:
        tap_service (pyvo.dal.tap.TAPService): TAP service that will be used for the query
        query (str): query to be run
        maxrec (int): define the maximum number of entries that a single query can return. Default is set
            by default.get_value('maxrec')

    Returns:
        astropy.table: result from the query to the TAP service

    """
    try:
        if maxrec is None:
            result_from_query = tap_service.search(query=query, maxrec=maxrec).to_table()
        else: 
            result_from_query = tap_service.search(query=query, maxrec=int(maxrec)).to_table()

    except (ValueError, DALQueryError, DALFormatError):
        msgs.warning('The query timed out. Trying with maxrec=100, but consider using `async` instead.')
        result_from_query = tap_service.search(query=query, maxrec=100).to_table()
        # msgs.error('The query failed. Trying `async` instead.')
        # msgs.warning('The query timed out. Trying `async` instead')
        # result_from_query = run_query_async(tap_service=tap_service, query=query, maxrec=maxrec)
    return result_from_query


def run_query_async(tap_service, query, maxrec=default.get_value('maxrec')):
    r"""Run an asynchronous query to TAP service and return result as an `astropy.Table`

    Args:
        tap_service (pyvo.dal.tap.TAPService): TAP service that will be used for the query
        query (str): query to be run
        maxrec (int): define the maximum number of entries that a single query can return. Default is set
            by default.get_value('maxrec')

    Returns:
        astropy.table: result from the query to the TAP service

    """
    tap_job = tap_service.submit_job(query=query, maxrec=maxrec)
    tap_job.run()
    # Wait for Executing
    tap_job.wait(phases=["EXECUTING", "ERROR", "ABORTED"], timeout=10.)
    msgs.info('The query to the tap_service is in the status: {}'.format(str(tap_job.phase)))
    # Wait for Completed
    tap_job.wait(phases=["COMPLETED", "ERROR", "ABORTED"], timeout=10.)
    msgs.info('The query to the tap_service is in the status: {}'.format(str(tap_job.phase)))
    # Fetch the results
    tap_job.raise_if_error()
    return tap_job.fetch_result().to_table()

# Query builders:

# General


def _create_comma_separated_list(list_of_strings):
    r"""Given a list of strings, returns them in a single string using comma as separator

    If `list_of_string is None, a `*` character is returned

    Args:
        list_of_strings (list, optional): list of `str`

    Returns:
        str: single string containing all the entries of the input list separated by a comma

    """
    if list_of_strings is None:
        final_string = '*'
    else:
        final_string = ' '
        for single_string in list_of_strings:
            final_string = final_string + ', ' + single_string
        # remove initial `,`
        final_string = final_string.strip()[2:]
    return final_string


# Catalogues:


def create_query_all_catalogues(all_versions=False, collections=None, tables=None):
    r"""Create TAP query that returns info on all catalogues in the ESO archive

    If `collections` or `tables` are not `None` only the query for the selected collections/tables is returned

    Args:
        all_versions (bool): if set to `True` also obsolete versions of the catalogues are listed
        collections (list, optional): list of `str` containing the names of the collections that are going to be
            selected
        tables (list, optional): list of `str` containing the table_name of the tables that are going to be selected

    Returns:
        str: string containing the query to obtain all catalogues present in the ESO archive

    """
    query_all_catalogues = '''
        SELECT 
            collection, title, version, table_name, filter, instrument, telescope, publication_date, 
            ref.description as description, number_rows, number_columns, rel_descr_url, acknowledgment,
            cat_id, mjd_obs, mjd_end, skysqdeg, bibliography, document_id, kc.from_column as from_column,
            k.target_table as target_table, kc.target_column as target_column, schema_name
        FROM
            TAP_SCHEMA.tables as ref
        LEFT OUTER JOIN 
            TAP_SCHEMA.keys as k on ref.table_name = k.from_table 
        AND
            k.target_table in (SELECT
                                   T.table_name
                               FROM 
                                   TAP_SCHEMA.tables as T
                               WHERE 3 in (SELECT 
                                               count(*) 
                                           FROM
                                               TAP_SCHEMA.columns
                                           WHERE
                                               table_name=T.table_name
                                           AND
                                               (ucd = 'pos.eq.ra;meta.main' OR
                                                ucd = 'pos.eq.dec;meta.main' OR
                                                ucd = 'meta.id;meta.main')
                                           )
                               )
        LEFT OUTER JOIN 
            TAP_SCHEMA.key_columns as kc on k.key_id=kc.key_id
        WHERE
            schema_name = 'safcat' '''

    if not all_versions:
        query_last_version_only = '''
        AND
            cat_id in (SELECT
                           t1.cat_id cat_id
                       FROM
                           TAP_SCHEMA.tables t1
                       LEFT OUTER JOIN
                           TAP_SCHEMA.tables t2 on (t1.title = t2.title AND
                                                    t1.version < t2.version)
                       WHERE
                           t2.title is null
                       )'''
        query_all_catalogues = query_all_catalogues + query_last_version_only

    if collections is not None:
        query_select_collections = '''
        AND
            ({}
            )'''.format(condition_collections_like(collections))
        query_all_catalogues = query_all_catalogues + query_select_collections

    if tables is not None:
        query_select_tables = '''
        AND
            ({}
            )'''.format(condition_tables_like(tables))
        query_all_catalogues = query_all_catalogues + query_select_tables

    return query_all_catalogues


def create_query_all_columns(collections=None, tables=None):
    r"""Create a query that get names (and corresponding info) of the columns present in a collection (or in a table)

    If `collections` and `tables` are `None` the query for the column of all collections in the ESO archive is
    returned.

    Args:
        collections (list, optional): list of `str` containing the names of the collections for which columns
            information will be returned
        tables (list, optional): list of `str` containing the table_name for which columns information will be returned

    Returns:
        str: string containing the query to obtain information on all columns present in the ESO archive

    """
    query_all_columns = '''
        SELECT 
            table_name, column_name, ucd, datatype, description, unit
        FROM 
            TAP_SCHEMA.columns
        WHERE 
            table_name in
                (SELECT
                    table_name
                 FROM
                    TAP_SCHEMA.tables
                 WHERE {}
                )
        AND
            ({}
            )'''.format(condition_collections_like(collections), condition_tables_like(tables))
    return query_all_columns


def create_query_table_base(table_name, columns=None, top=None):
    r"""Create a query to return selected columns from a table

    Args:
        columns (list, optional): list of `str` containing the columns that will be queried
        table_name (str): name of the table that will be queried
        top (int, optional): number of rows that will be returned. If set to `None` all rows will be returned

    Returns:
        str: query to obtain data from a table

    """

    if top is not None:
        query = '''
            SELECT TOP {}
                {} 
            FROM 
                {}'''.format(int(top), _create_comma_separated_list(columns), table_name)

    else:
        query = '''
                SELECT 
                    {} 
                FROM 
                    {}'''.format(_create_comma_separated_list(columns), table_name)
    
    return query


def condition_source_ids_like(source_ids, source_id_name='SOURCEID'):
    r"""Create the WHERE `source_id_name` = condition string for a query

    Args:
        source_ids (list): list of `str` containing the source_ids  for which columns
            information will be returned
        source_id_name (str): name of the source id column in a catalogue

    Returns:
        str: String containing the WHERE `source_id_name` = condition for a query

    """
    if source_ids is None:
        return ''
    else:
        condition_source_ids = '''
            WHERE
                ('''
        list_source_ids = source_ids.copy()
        for source_id in list_source_ids:
            condition_source_id = '''
                {} = {} OR'''.format(source_id_name, source_id)
            condition_source_ids = condition_source_ids + condition_source_id
            # remove the last OR
            condition_source_ids = condition_source_ids[:-3] + '''
                )'''
    return condition_source_ids

# Observations


def create_query_obscore_base(top=None, columns=None):
    r"""Create the base string for a query to `ivoa.ObsCore`

    If `columns` is set to None, default list of columns is set by the global variable `COLUMNS_FROM_OBSCORE`

    Args:
        columns (list, optional): list of `str` containing the columns that will be queried

    Returns:
        str: Base for the `ivoa.ObsCore` queries

    """

    if columns is None:
        columns = COLUMNS_FROM_OBSCORE

    if top is not None:
        query_base = '''
            SELECT TOP {}
                {}
            FROM
                ivoa.ObsCore'''.format(top, _create_comma_separated_list(columns))

    else: 
        query_base = '''
                SELECT
                    {}
                FROM
                    ivoa.ObsCore'''.format(_create_comma_separated_list(columns))
    return query_base


def create_query_obscore_all_columns():
    r"""Create a query that get names (and corresponding info) of the columns present `ivoa.ObsCore`

    Returns:
        str: string containing the query to obtain information on all columns present in `ivoa.ObsCore`

    """
    query_obscore_all_columns = '''
        SELECT 
            table_name, column_name, ucd, datatype, description, unit, utype
        FROM 
            TAP_SCHEMA.columns
        WHERE 
            table_name='ivoa.ObsCore' '''
    return query_obscore_all_columns


def condition_intersects_ra_dec(ra, dec, radius=None):
    r"""Create the WHERE INTERSECTS condition string for a query

    Args:
        ra (float): RA of the target in degrees and in the ICRS system
        dec (float): Dec of the target in degrees and in the ICRS system
        radius (float, optional): Search radius in arcsec. If set to `None` no radius will be considered in the
            INTERSECT condition

    Returns:
        str: String containing the WHERE INTERSECT condition for a query

    """
    if radius is None:
        query_intersects_ra_dec = '''
            WHERE
                INTERSECTS(POINT('ICRS',{0},{1}), s_region) = 1'''.format(str(ra),
                                                                          str(dec))
    else:
        query_intersects_ra_dec = '''
            WHERE
                INTERSECTS(s_region,CIRCLE('ICRS',{0},{1},{2}/3600.)) = 1'''.format(str(ra),
                                                                                    str(dec),
                                                                                    str(radius))
    return query_intersects_ra_dec


def condition_intersects_polygon(polygon):
    r"""Create the WHERE INTERSECTS polygon condition string for a query

    Args:
        polygon (str): coordinates of the vertices of the polygon in the sky

    Returns:
        str: String containing the WHERE INTERSECT condition for a query

    """
    query_intersects_polygon = '''
            WHERE
                INTERSECTS(s_region,POLYGON('', {0} )) = 1'''.format(polygon)
    return query_intersects_polygon


# Conditions

# ToDo, these can be compacted in one single function


def condition_tables_like(tables=None):
    r"""Create a `LIKE` - `OR` condition over a list of table_names

    If `tables` is `None` the query the conditions will be substitute with '%' meaning that will have no effect

    Args:
        tables (list): list of `str` containing the table_name for which columns information will be returned

    Returns:
        str: set of conditions in the format: table_name LIKE --- OR

    """
    condition_tables = ''''''
    if tables is None:
        list_tables = ['%']
    else:
        list_tables = tables.copy()
    for table in list_tables:
        condition_table = '''
                    table_name LIKE '{}' OR'''.format(table)
        condition_tables = condition_tables + condition_table
    # remove the last OR
    condition_tables = condition_tables[:-3]
    return condition_tables


def condition_instruments_like(instruments=None):
    r"""Create condition string to select only specific instruments in `ivoa.ObsCore`

    Args:
        instruments (list): limit the search to the selected list of
            instruments (e.g., `XSHOOTER`)

    Returns:
        str: string containing the `instrument_name` condition for a query

    """
    condition_instruments = '''
            AND
                ('''
    if instruments is None:
        list_instruments = ['%']
    else:
        list_instruments = instruments.copy()
    for instrument in list_instruments:
        condition_instrument = '''
                instrument_name LIKE '{}' OR'''.format(instrument)
        condition_instruments = condition_instruments + condition_instrument
    # remove the last OR
    condition_instruments = condition_instruments[:-3] + '''
                )'''
    return condition_instruments


def condition_data_types_like(data_types=None):
    r"""Create condition string to select only specific data product types in `ivoa.ObsCore`

    Args:
        data_types (list): limit the search to the selected list of dataproduct types (e.g., `spectrum`)

    Returns:
        str: string containing the `dataproduct_type` condition for a query

    """
    condition_data_types = '''
            AND
                ('''
    if data_types is None:
        list_data_types = ['%']
    else:
        list_data_types = data_types.copy()
    for data_type in list_data_types:
        condition_data_type = '''
                dataproduct_type LIKE '{}' OR'''.format(data_type)
        condition_data_types = condition_data_types + condition_data_type
    # remove the last OR
    condition_data_types = condition_data_types[:-3] + '''
                )'''
    return condition_data_types


def condition_collections_like(collections=None):
    r"""Create a `LIKE` - `OR` condition over a list of collections

    If `collections` is `None` the query the conditions will be substitute with '%' meaning that will have no effect

    Args:
        collections (list, optional): list of `str` containing the names of the collections for which columns
            information will be returned

    Returns:
        str: set of conditions in the format: collection LIKE --- OR

    """
    condition_collections = ''''''
    if collections is None:
        list_collections = ['%']
    else:
        list_collections = collections.copy()
    for collection in list_collections:
        condition_collection = '''
                    collection LIKE '{}' OR'''.format(collection)
        condition_collections = condition_collections + condition_collection
    # remove the last OR
    condition_collections = condition_collections[:-3]
    return condition_collections

def condition_em_like(em_min=None, em_max=None, enclosed=True):
    r"""Create condition string to select only specific em_min and em_max in `ivoa.ObsCore`

    Args:
        em_min (float): limit the search to the selected list of em_min (e.g., `>1000`)
        em_max (float): limit the search to the selected list of em_max (e.g., `<10000`)
        enclosed (bool): if set to `True` the condition will be `em_min > em_min` and `em_max < em_max`, otherwise
            the condition will be `em_min < em_min` and `em_max > em_max`

    Returns:
        str: string containing the `em_min` and `em_max` condition for a query

    """

    condition_em = ''''''
    if enclosed: 

        if em_min is not None:
            condition_em += f'''
                AND 
                    em_min>{em_min}'''
            
        if em_max is not None:
            condition_em += f'''
                AND 
                    em_max<{em_max}'''
    else:
        if em_min is not None:
            condition_em += f'''
                AND 
                    em_min<{em_min}'''
            
        if em_max is not None:
            condition_em += f'''
                AND 
                    em_max>{em_max}'''
        
    return condition_em

def condition_snr_like(snr=None):
    r"""Create condition string to select only specific snr in `ivoa.ObsCore`

    Args:
        snr (list): limit the search to the selected list of snr (e.g., `>10`)

    Returns:
        str: string containing the `snr` condition for a query

    """
    condition_snr = ''''''
    if snr is not None:
        condition_snr = f'''
            AND 
                snr>{snr}'''
        
    return condition_snr

def condition_em_res_power_like(em_res_power=None):
    r"""Create condition string to select only specific em_res_power in `ivoa.ObsCore`

    Args:
        em_res_power (list): limit the search to the selected list of em_res_power (e.g., `>1000`)

    Returns:
        str: string containing the `em_res_power` condition for a query

    """
    condition_em_res_power = ''''''
    if em_res_power is not None:
        condition_em_res_power = f'''
            AND 
                em_res_power>{em_res_power}'''
        
    return condition_em_res_power

def condition_order_by_like(order_by=None, order='ascending'):
    r"""Create condition string to select only specific order_by in `ivoa.ObsCore`

    Args:
        order_by (list): limit the search to the selected list of order_by (e.g., `t_exptime`)

    Returns:
        str: string containing the `order_by` condition for a query

    """
    condition_order_by = ''''''
    if order_by is not None:

        if order == 'descending':
            condition_order_by = f'''
                ORDER BY 
                    {order_by} DESC'''

        else: 
            condition_order_by = f'''
                ORDER BY 
                    {order_by} ASC'''
        
    return condition_order_by

def conditions_dict_like(conditions_dict=None):
    
    conditions_query = ''''''
    if conditions_dict is not None:
        for key, value in conditions_dict.items():
            conditions_query += f'''
                WHERE
                    {key}={value}''' 
            
    return conditions_query
