#!/usr/bin/env python
"""
REFlEx - Utils - Assign discharge hydro-hydra
__date__ = '20230315'
__version__ = '1.0.1'
__author__ =
        'Andrea Libertino (andrea.libertino@cimafoundation.org)'
        'Lorenzo Campo (lorenzo.campo@cimafoundation.org)'
__library__ = 'REFlEx'
General command line:
### python reflex_step1_hydro_derivatives.py -settings_file "settings.json"
Version(s):
20230315 (1.0.1) --> Add check on available neighbours when system does not converge
20220530 (1.0.0) --> Beta release
"""

# -------------------------------------------------------------------------------------
# Import python libraries
from rasterio.features import shapes
from shapely.geometry import shape
import geopandas as gpd
import rioxarray as rxr
import numpy as np
import os
import pandas as pd
import logging, json
from argparse import ArgumentParser
import time
from sklearn.neighbors import BallTree

# -------------------------------------------------------------------------------------
# Script Main
def main():
    # -------------------------------------------------------------------------------------
    # Version and algorithm information
    alg_name = 'REFlEx - Utils - Assign discharge hydro-hydra '
    alg_version = '1.0.1'
    alg_release = '2023-03-15'
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    alg_settings = get_args()

    # Set algorithm settings
    data_settings = read_file_json(alg_settings)
    settings_path = data_settings["paths"]
    domain = data_settings["domain"]["name"]

    # Fill paths
    settings_path_filled = {}
    for key in settings_path.keys():
        settings_path_filled[key] = settings_path[key].replace("{domain}",domain)

    os.makedirs(settings_path_filled["output_folder"], exist_ok = True)
    os.makedirs(settings_path_filled["ancillary_folder"], exist_ok = True)
    os.makedirs(settings_path_filled["log_folder"], exist_ok=True)

    # Set algorithm logging
    os.makedirs(settings_path_filled["log_folder"], exist_ok=True)
    set_logging(
        logger_file=os.path.join(settings_path_filled["log_folder"], "assignation_" + domain + "_log.txt"))
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    logging.info(' ============================================================================ ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> START ... ')
    logging.info(' ')

    logging.info(" --> Domain : " + domain)

    # Time algorithm information
    start_time = time.time()

    # Read discharge values
    col_names = ["Q" + str(i) for i in data_settings["discharge"]["T_vals"]]
    discharge_vals = pd.read_csv(data_settings["discharge"]["file"], sep=data_settings["discharge"]["sep"], index_col=0,
                                 usecols=[data_settings["discharge"]["index_col"]] + data_settings["discharge"][
                                     "q_cols"])

    # Manage hydro features
    logging.info(" --> Manage hydrological basins and network...")
    with rxr.open_rasterio(settings_path_filled["basins_hydro_file"]) as raster:
        image = raster.values
        crs = raster.rio.crs
        if crs is None:
            logging.warning(" ---> WARNING! The hydrological basin shapefile is not properly georeferenced! EPSG:4326 is considered!")
            crs = "EPSG:4326"
        list_pop = [
            {'bas_cd': value, 'area_skm':-9999.0, 'geometry': shape(shp)}
            for i, (shp, value)
            in enumerate(shapes(image.astype('uint16'), transform=raster.rio.transform()))
        ]
        basin_hydro = gpd.GeoDataFrame(list_pop, crs=crs).to_crs(epsg=4326)
        basin_hydro = basin_hydro.dissolve(by='bas_cd', as_index=False)
        accum = rxr.open_rasterio(settings_path_filled["basin_hydro_area"]).reindex_like(raster, method="nearest")

        if data_settings["flags"]["drop_all_zeroes_quantiles"]:
            allzeroes_areas = discharge_vals.index[np.nansum(discharge_vals.values, axis=1) == 0]
            basin_hydro = basin_hydro.loc[basin_hydro["bas_cd"].isin(allzeroes_areas.tolist()) == False]

        basin_hydro_centroid = gpd.GeoDataFrame(geometry=basin_hydro.centroid)
        basin_hydro_centroid["area_skm"] = -9999.0

        for ind_bas, bas in zip(basin_hydro.index, basin_hydro["bas_cd"]):
            try:
                basin_hydro.loc[ind_bas, "area_skm"] = np.max(accum.values[image == bas])
                basin_hydro_centroid.loc[ind_bas, "area_skm"] = np.max(accum.values[image == bas])
                basin_hydro_centroid.loc[ind_bas, "bas_cd"] = bas
            except:
                basin_hydro = basin_hydro.drop(ind_bas, axis=0)
                basin_hydro_centroid = basin_hydro_centroid.drop(ind_bas, axis=0)

        basin_hydro_centroid = basin_hydro_centroid.set_index(np.arange(0, len(basin_hydro_centroid)))
        basin_hydro.to_file(os.path.join(settings_path_filled["ancillary_folder"], "basin_dissolved.shp"))
        basin_hydro_centroid.to_file(os.path.join(settings_path_filled["ancillary_folder"], "basin_hydro_centroid_dissolved.shp"))
    logging.info(" --> Manage hydrological basins and network...DONE")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Manage reflex features
    logging.info(" --> Manage reflex basins and network...")
    basin_reflex = gpd.read_file(settings_path_filled["basin_reflex_file"])
    stream_reflex = gpd.read_file(settings_path_filled["network_reflex_file"])

    basin_reflex_centroid = gpd.GeoDataFrame(geometry=basin_reflex.centroid).set_index(basin_reflex["stream"].values)

    for bas in stream_reflex["stream"].values:
        try:
            basin_reflex_centroid.loc[bas,"area_skm"] = stream_reflex.loc[stream_reflex["stream"].values==bas, "flowAcc_sk"].values
        except:
            pass

    basin_reflex_centroid = basin_reflex_centroid.dropna(axis=0, how='any')
    logging.info(" --> Manage reflex basins and network...DONE")
    # -------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------
    from copy import deepcopy
    input_values = deepcopy(basin_reflex_centroid.index)
    out_df = pd.DataFrame(index=input_values, columns=["near_bas", "diff_facc", "dist_km"])

    max_attempts = data_settings["algorithm"]["max_attempts"]
    k_neighbors = data_settings["algorithm"]["number_of_neighbours"]
    attempt = 1
    last_run = False

    # Assign discharge
    while len(input_values) > 0 and attempt <= max_attempts:
        logging.info(" --> Assign hydrological and hydraulic network... Attempt: " + str(attempt))
        in_pts = [(x,y) for x,y in zip(basin_reflex_centroid.loc[input_values].geometry.x , basin_reflex_centroid.loc[input_values].geometry.y)]
        qry_pts = [(x,y) for x,y in zip(basin_hydro_centroid.geometry.x , basin_hydro_centroid.geometry.y)]

        tab_close, tab_dist = get_nearest(in_pts, qry_pts, k_neighbors=k_neighbors)

        for i_row, i in enumerate(input_values):
            logging.info(" ---> Compute stream " + str(i))
            nearest_indeces = tab_close[i_row,:]
            nearest_distances = tab_dist[i_row,:]
            diff_flow_acc_nearest = np.abs(basin_reflex_centroid.loc[i,"area_skm"] - basin_hydro_centroid.loc[nearest_indeces, "area_skm"])/basin_reflex_centroid.loc[i,"area_skm"]
            out_df.loc[i,"near_bas"]= basin_hydro_centroid.loc[diff_flow_acc_nearest.idxmin(), "bas_cd"]
            out_df.loc[i,"diff_facc"] = diff_flow_acc_nearest.min()
            out_df.loc[i,"dist_km"] = nearest_distances[nearest_indeces==diff_flow_acc_nearest.idxmin()][0]
        logging.info(" --> Assign hydrological and hydraulic network... Attempt: " + str(attempt) + " DONE")

        input_values = out_df.loc[out_df["diff_facc"].values > data_settings["algorithm"]["diff_flowacc_accept"]].index
        attempt = attempt + 1
        k_neighbors = k_neighbors * 2

        if last_run == True:
            break

        # If search for a larger than the available basin, colalpse to the maximum and do the last try
        if k_neighbors > len(qry_pts):
            k_neighbors = len(qry_pts)
            last_run = True

    if len(input_values) > 0:
        out_df.loc[input_values] = np.nan
        logging.warning("WARNING! " + str(len(input_values)) + " streams have not been assigned")

    # Manage reflex features
    #logging.info(" --> Assign hydrological and hydraulic network...")
    #in_pts = [(x, y) for x, y in zip(basin_reflex_centroid.geometry.x, basin_reflex_centroid.geometry.y)]
    #qry_pts = [(x, y) for x, y in zip(basin_hydro_centroid.geometry.x, basin_hydro_centroid.geometry.y)]

    #tab_close, tab_dist = get_nearest(in_pts, qry_pts, k_neighbors=data_settings["algorithm"]["number_of_neighbours"])

    #out_df = pd.DataFrame(index=basin_reflex_centroid.index, columns=["near_bas", "diff_facc", "dist_km"])

    #for i_row, i in enumerate(out_df.index):
    #    logging.info(" ---> Compute stream " + str(i))
    #    nearest_indeces = tab_close[i_row, :]
    #    nearest_distances = tab_dist[i_row, :]
    #    diff_flow_acc_nearest = np.abs(
    #        basin_reflex_centroid.loc[i, "area_skm"] - basin_hydro_centroid.loc[nearest_indeces, "area_skm"]) / \
    #                            basin_reflex_centroid.loc[i, "area_skm"]
    #    out_df.loc[i, "near_bas"] = basin_hydro_centroid.loc[diff_flow_acc_nearest.idxmin(), "bas_cd"]
    #    out_df.loc[i, "diff_facc"] = diff_flow_acc_nearest.min()
    #    out_df.loc[i, "dist_km"] = nearest_distances[nearest_indeces == diff_flow_acc_nearest.idxmin()][0]
    #logging.info(" --> Assign hydrological and hydraulic network...DONE")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Prepare output
    out_stream_reflex = stream_reflex.copy()

    for i_row in out_stream_reflex.index.values:
        stream = out_stream_reflex.loc[i_row, "stream"]
        near_stream = out_df.loc[stream, "near_bas"]
        if not np.isnan(near_stream):
            out_stream_reflex.loc[i_row, col_names] = discharge_vals.loc[near_stream].values
        else:
            out_stream_reflex.loc[i_row, col_names] = -9999

    if not "Q0" in col_names:
        out_stream_reflex["Q0"] = 0
        col_names = ["Q0"] + col_names
    out_stream_reflex.to_file(os.path.join(settings_path_filled["output_folder"], domain + "_assigned.shp"))
    out_df.to_csv(os.path.join(settings_path_filled["output_folder"], domain + "_assigned_stats.csv"))
    out_stream_reflex[["stream"] + col_names].to_csv(os.path.join(settings_path_filled["output_folder"], domain + "_assigned.csv"), index=False)
    logging.info(" --> Assign hydrological and hydraulic network...DONE")

    if data_settings["flags"]["clear_ancillary"]:
        os.system("rm " + os.path.join(settings_path_filled["ancillary_folder"],'*'))
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    time_elapsed = round(time.time() - start_time, 1)

    logging.info(' ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    logging.info(' ==> ... END')
    logging.info(' ==> Bye, Bye')
    logging.info(' ============================================================================ ')
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
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


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-settings_file', action="store", dest="alg_settings")
    parser_values = parser_handle.parse_args()

    if parser_values.alg_settings:
        alg_settings = parser_values.alg_settings
    else:
        alg_settings = 'configuration.json'

    return alg_settings

# -------------------------------------------------------------------------------------

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
def get_nearest(src_points, candidates, k_neighbors=5):
    """
    Find nearest neighbors for all source points from a set of candidate points
    modified from: https://automating-gis-processes.github.io/site/notebooks/L3/nearest-neighbor-faster.html
    """

    # Create tree from the candidate points
    tree = BallTree(candidates, leaf_size=15, metric='euclidean')

    # Find closest points and distances
    distances, indices = tree.query(src_points, k=k_neighbors)

    # Transpose to get distances and indices into arrays
    distances = distances.transpose()
    indices = indices.transpose()

    # Get closest indices and distances (i.e. array at index 0)
    # note: for the second closest points, you would take index 1, etc.
    closest = np.stack([indices[i] for i in np.arange(0, k_neighbors, 1)], axis=0).T
    dist = np.stack([distances[i] for i in np.arange(0, k_neighbors, 1)], axis=0).T

    # Return indices and distances
    return closest, dist * 111
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------

