"""
Library Features:

Name:          reflex_tools_utils
Author(s):     Mauro Arcorace (mauro.arcorace@cimafoundation.org)
               Fabio Delogu (fabio.delogu@cimafoundation.org)
               Andrea Libertino (andrea.libertino@cimafoundation.org)
Date:          '20220404'
Version:       '2.0.0'
"""
########################################################################################################################

import binascii, math, json, logging
import os
import subprocess
import sys

# ---------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------
# Return proper elapsed time from a given total execution time in seconds
def Give_Elapsed_Time(time):
    if (time > 60):
        elapsed_time = time / 60.0
        time_units = "mins"
    elif (time > 3600):
        elapsed_time = time / 3600.0
        time_units = "hours"
    elif (time > 86400):
        elapsed_time = time / 86400.0
        time_units = "days"
    else:
        elapsed_time = time
        time_units = "seconds"

    return elapsed_time, time_units
# ---------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------
def Start_GRASS_py3(grass_bin, epsg_code, grass_temp_db_name, out_folder):
    # Check GRASS config path
    startcmd = grass_bin + ' --config path'
    p = subprocess.Popen(startcmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    # Set GISBASE environment variable
    gisbase = out.decode("utf-8").strip('\n')  #######
    os.environ['GISBASE'] = gisbase

    # Set GRASS-Python environment
    gpydir = os.path.join(gisbase, "etc", "python")
    sys.path.append(gpydir)

    # Define temp GRASS database
    grass_temp_db = os.path.join(out_folder, grass_temp_db_name)
    if not os.path.exists(grass_temp_db):
        os.makedirs(grass_temp_db)

    # location/mapset: use random names for batch jobs
    location = binascii.hexlify(os.urandom(16)).decode("utf-8")
    mapset = 'PERMANENT'
    location_path = os.path.join(grass_temp_db, location)

    # Create new location for mapset
    startcmd = grass_bin + ' -c epsg:' + epsg_code + ' -e ' + location_path
    # startcmd = grass_bin + ' -c ' + input_dem_geotiff + ' -e ' + location_path

    p = subprocess.Popen(startcmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    # Set GISDBASE environment variable
    os.environ['GISDBASE'] = grass_temp_db

    # Set path to GRASS libs
    path = os.getenv('LD_LIBRARY_PATH')
    dir = os.path.join(gisbase, 'lib')
    if path:
        path = dir + os.pathsep + path
    else:
        path = dir
    print("--> LD_LIBRARY_PATH environment variable: %s" % str(path))

    # Set os environments
    os.environ['LD_LIBRARY_PATH'] = path
    os.environ['LANG'] = 'en_US'
    os.environ['LOCALE'] = 'C'

    return gisbase, grass_temp_db, location, mapset

# ---------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------
# Based on lat and lng, return best utm epsg-code
# Credits to https://gis.stackexchange.com/users/33092/alex?tab=profile
def convert_wgs_to_utm(lon, lat):
    utm_band = str((math.floor((lon + 180) / 6 ) % 60) + 1)
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
    if lat >= 0:
        epsg_code = '326' + utm_band
        return epsg_code
    epsg_code = '327' + utm_band
    return epsg_code
# ---------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Method to read file json
def read_file_json(file_name):
    env_ws = {}
    for env_item, env_value in os.environ.items():
        env_ws[env_item] = env_value

    with open(file_name, "r") as file_handle:
        json_block = []
        for file_row in file_handle:

            for env_key, env_value in env_ws.items():
                env_tag = '$' + env_key
                if env_tag in file_row:
                    env_value = env_value.strip("'\\'")
                    file_row = file_row.replace(env_tag, env_value)
                    file_row = file_row.replace('//', '/')

            # Add the line to our JSON block
            json_block.append(file_row)

            # Check whether we closed our JSON block
            if file_row.startswith('}'):
                # Do something with the JSON dictionary
                json_dict = json.loads(''.join(json_block))
                # Start a new block
                json_block = []

    return json_dict
# ----------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to set logging information
def set_logging(logger_file='log.txt', logger_format=None):
    if logger_format is None:
        logger_format = '%(asctime)s %(name)-12s %(levelname)-8s ' \
                        '%(filename)s:[%(lineno)-6s - %(funcName)20s()] %(message)s'

    # Remove old logging file
    if os.path.exists(logger_file):
        os.remove(logger_file)

    # Set level of root debugger
    logging.root.setLevel(logging.INFO)

    # Open logging basic configuration
    logging.basicConfig(level=logging.INFO, format=logger_format, filename=logger_file, filemode='w')

    # Set logger handle
    logger_handle_1 = logging.FileHandler(logger_file, 'w')
    logger_handle_2 = logging.StreamHandler()
    # Set logger level
    logger_handle_1.setLevel(logging.INFO)
    logger_handle_2.setLevel(logging.INFO)
    # Set logger formatter
    logger_formatter = logging.Formatter(logger_format)
    logger_handle_1.setFormatter(logger_formatter)
    logger_handle_2.setFormatter(logger_formatter)
    # Add handle to logging
    logging.getLogger('').addHandler(logger_handle_1)
    logging.getLogger('').addHandler(logger_handle_2)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
