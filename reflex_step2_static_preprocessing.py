#!/usr/bin/env python
"""
REFlEx - Step2 - Static data preprocessing
__date__ = '20230101'
__version__ = '2.0.3'
__author__ =
        'Mauro Arcorace' (mauro.arcorace@cimafoundation.org',
        'Alessandro Masoero (alessandro.masoero@cimafoundation.org',
        'Valerio Basso',
        'Alessia Matanò',
        'Giulia Bruno (giulia.bruno@cimafoundation.org)',
        'Andrea Libertino (andrea.libertino@cimafoundation.org)'
__library__ = 'REFlEx'
General command line:
### python reflex_step2_static_preprocessing.py -log_file "/path/to/log.txt" -settings_file "settings.json" -base_path "/path/to/base_folder"
Version(s):
20190220 (1.0.0) --> Beta release
20220406 (2.0.0) --> Full revision - Mask file produced as unique shapefile for each stream
                                     Add the possibility of choose which concentration time equations to include
                                     Revised stream management (volume extimation moved to step 4)
                                     Automatic selection of best epsg for proj
                                     Parallel implementation
20220726 (2.0.1) --> Fix basin delineation procedure
20230101 (2.0.3) --> Optimized multiprocessing
                     Fixed pfafstetter codification
                     Optimized singleprocessing for big domains
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'REFlEx - STEP 2 - Static Data Processing'
alg_version = '2.0.3'
alg_release = '2023-01-01'
# Algorithm parameter(s)
time_format = '%Y%m%d%H%M'
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Import python libraries
import os
import logging
from argparse import ArgumentParser
import time
from lib.reflex_tools_utils import Give_Elapsed_Time, convert_wgs_to_utm, set_logging, read_file_json
from lib.reflex_tools_basins import create_masks, calculate_basins_stats
import geopandas as gpd
import numpy as np
import pandas as pd
import sys
from multiprocessing import Pool, Manager, cpu_count, set_start_method, get_context
from numba import config
from copy import deepcopy


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Script Main
def main():
    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    config_file, alg_log, base_path = get_args()

    # Set algorithm settings
    data_settings = read_file_json(config_file)

    os.makedirs(os.path.dirname(alg_log), exist_ok=True)
    set_logging(logger_file=alg_log)

    overwrite_mode = True
    quiet_mode = True
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    logging.info(' ============================================================================ ')
    logging.info(' ================ REFlEx - Rapid Estimation od FLood EXtent ================= ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> START ... ')
    logging.info(' ')
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # GET PARAMETERS FROM CONFIG FILE
    # -------------------------------

    # Path to executable and version
    script_ver = data_settings["algorithm"]["version"]
    grass_bin = data_settings["algorithm"]["grass_bin"]

    # Domain and DEM info
    acp_x = data_settings["domain"]["name"].upper()
    domain_name = data_settings["domain"]["name"]
    rrs_km = data_settings["domain"]["res_km"]
    in_res_DD = data_settings["domain"]["res_deg"]
    rrs_m = rrs_km * 1000
    rrs = data_settings["domain"]["res_str"]
    source_epsg = data_settings["domain"]["source_epsg"]
    target_epsg = data_settings["domain"]["target_epsg"]

    # Previous step settings
    drain_method_streams = data_settings["step_1"]["stream_definition_method"]

    # Step settings
    step0_dir_name = data_settings["step_0"]["dir_name"].replace("{base_path}", base_path)
    step1_dir_name = data_settings["step_1"]["dir_name"].replace("{base_path}", base_path)
    step2_dir_name = data_settings["step_2"]["dir_name"].replace("{base_path}", base_path)

    buffer_watershed_distance_cell = data_settings["step_2"]["buffer_for_coastal_expansion_cells"]
    c_time = data_settings["step_2"]["concentration_time"]
    try:
        max_attempts = data_settings["step_2"]["max_attempts_make_mask"]
    except:
        max_attempts = 3

    if data_settings["step_2"]["multiprocessing"]["enable"]:
        process_max = data_settings["step_2"]["multiprocessing"]["max_cores"]
        if process_max is None:
            process_max = cpu_count() - 1
    else:
        process_max = 1

    chunk_size = data_settings["step_2"]["multiprocessing"]["chunk_size"]

    # Set numba execution
    config.THREADING_LAYER_PRIORITY = ["tbb", "omp", "workqueue"]
    ################################################################################

    ################################################################################

    # CREATE OUTPUT DIRECTORIES
    # -------------------------

    # Define output MAIN directory
    grass_output_db = step2_dir_name
    logging.info('--> Output main directory: %s' % grass_output_db)
    if not os.path.exists(grass_output_db):
        os.makedirs(grass_output_db)

    # Define output vector sub-directory
    grass_output_db_vct_dirname = 'vct'
    grass_output_db_vct = os.path.join(grass_output_db, grass_output_db_vct_dirname)
    logging.info('--> Output vector sub-directory: %s' % grass_output_db_vct)
    if not os.path.exists(grass_output_db_vct):
        os.makedirs(grass_output_db_vct)

    # Define output document sub-directory
    grass_output_db_txt_dirname = 'txt'
    grass_output_db_txt = os.path.join(grass_output_db, grass_output_db_txt_dirname)
    logging.info('--> Output txt or cvs sub-directory: %s' % grass_output_db_txt)
    if not os.path.exists(grass_output_db_txt):
        os.makedirs(grass_output_db_txt)

    # Define output document sub-directory
    grass_output_db_tmp_dirname = 'tmp'
    grass_output_db_tmp = os.path.join(grass_output_db, grass_output_db_txt_dirname)
    logging.info('--> Output txt or cvs sub-directory: %s' % grass_output_db_txt)
    if not os.path.exists(grass_output_db_tmp):
        os.makedirs(grass_output_db_tmp)

    # -------------------------

    ################################################################################

    ################################################################################

    code_start_time = time.time()

    # Define output GRASS database
    grass_step1_db = step1_dir_name
    grass_step0_db = step0_dir_name

    # CALCULATE STREAMS PFAFSTETTER HIERARCHY
    # ---------------------------------------

    # ---------------------------------------
    shape_in = {}
    shape_in["shp_stream_pfaf"] = os.path.join(grass_step1_db, "vct",
                                               'v_' + domain_name + '_' + rrs + '_streams_mod_pf.shp')
    shape_in["shp_subbasins"] = os.path.join(grass_step1_db, "vct",
                                             'v_' + domain_name + '_' + rrs + '_subbasins_dissolved.shp')
    # -------------------------------

    # -------------------------------
    # read basin geodataframe
    basin_gdf = gpd.read_file(shape_in["shp_subbasins"]).copy()
    additional_cols = ["lat_min", "lat_max", "lon_min", "lon_max"]
    for field in additional_cols:
        basin_gdf[field] = np.nan
    basin_gdf["mask_type"] = ""

    logging.info(' --> Create masks for ' + str(len(basin_gdf)) + " basins...")

    centroid_x = basin_gdf.total_bounds[0] + (basin_gdf.total_bounds[2] - basin_gdf.total_bounds[0]) / 2
    centroid_y = basin_gdf.total_bounds[1] + (basin_gdf.total_bounds[3] - basin_gdf.total_bounds[1]) / 2

    prj_epsg = convert_wgs_to_utm(centroid_x, centroid_y)

    logging.info(' --> Selected EPSG for area calculation: ' + str(prj_epsg))

    area_km = pd.DataFrame(data=basin_gdf.to_crs('epsg:%s' % str(prj_epsg)).area / (10 ** 6)).set_index(
        basin_gdf["stream"].values)

    # read stream geodtaframe
    streams_gdf = gpd.read_file(shape_in["shp_stream_pfaf"]).copy()
    additional_cols = ["tconc", "tpeak", "treces", "str_len_km", "flowAcc_skm"]
    for field in additional_cols:
        streams_gdf[field] = np.nan
    streams_gdf["str_len_km"] = streams_gdf.to_crs('epsg:%s' % str(prj_epsg)).geometry.length.values / (10 ** 3)

    stream_ids = streams_gdf["stream"].values

    for stream in stream_ids:
        streams_gdf.loc[streams_gdf["stream"] == stream, "area_km2"] = area_km.loc[stream].values
        streams_gdf.loc[streams_gdf["stream"] == stream, "distance"] = min(15000,
                                                                           area_km.loc[stream].values * 0.13 + 2000)

    # calculate min and max pfaf per macrobasin and initialise settings
    pfaf_limits = {}
    for iMBID in streams_gdf['MBID_new'].unique():
        pd_tmp = streams_gdf.loc[streams_gdf['MBID_new'] == iMBID].copy()
        maxpfA = np.nanmax(pd_tmp['pfA'])
        minpfA = np.nanmin(pd_tmp['pfA'])
        pfaf_limits[iMBID] = [minpfA, maxpfA]

    masks_settings = {}
    masks_settings["buffer_distance"] = float(in_res_DD) * float(buffer_watershed_distance_cell)
    masks_settings["pfaf_limits"] = pfaf_limits
    masks_settings["masks_folder"] = os.path.join(grass_output_db_vct, 'masks')
    masks_settings["input_epsg"] = target_epsg
    masks_settings["tmp_folder"] = grass_output_db_tmp

    os.makedirs(masks_settings["masks_folder"], exist_ok=True)

    # Find upstream basins
    shreve_unique = np.unique(streams_gdf["shreve"].values)

    upst_basin = {}
    stream_to_delete = []

    for shreve in shreve_unique:
        sub_lev = streams_gdf.loc[streams_gdf["shreve"].values == shreve]
        for _, row in sub_lev.iterrows():
            upst_basin[row["stream"]] = [row["stream"]]
            if row["prev_str01"] > 0:
                try:
                    upst_basin[row["stream"]] = upst_basin[row["stream"]] + upst_basin[row["prev_str01"]]
                except:
                    stream_to_delete = stream_to_delete + [row["stream"]]
            if row["prev_str02"] > 0:
                try:
                    upst_basin[row["stream"]] = upst_basin[row["stream"]] + upst_basin[row["prev_str02"]]
                except:
                    if not row["stream"] in stream_to_delete:
                        stream_to_delete = stream_to_delete + [row["stream"]]

    masks_settings["upstream_basins"] = upst_basin

    logging.warning(" --> WARNING! Streams " + ", ".join(
        [str(i) for i in stream_to_delete]) + " have been deleted because upstream branches are missing!")

    if len(stream_to_delete) > 0:
        streams_gdf = streams_gdf[~streams_gdf['stream'].isin(stream_to_delete)]
        basin_gdf = basin_gdf[~basin_gdf['stream'].isin(stream_to_delete)]

    manager = Manager()
    d = manager.dict()
    d["streams_gdf"] = streams_gdf
    d["basins_gdf"] = basin_gdf

    chunks = [streams_gdf["stream"].values[i:i + chunk_size] for i in
              range(0, len(streams_gdf["stream"].values), chunk_size)]

    for chunk in chunks:
        logging.info(" ---> Launching chunk from " + str(min(chunk)) + " to " + str(max(chunk)))
        exec_pool = get_context('spawn').Pool(process_max)
        for stream in chunk:  # stream_ids:
            exec_pool.apply_async(create_masks, args=(stream, masks_settings, d))
        exec_pool.close()
        exec_pool.join()

    missing_masks = []
    for stream in streams_gdf["stream"].values:
        if not os.path.isfile(os.path.join(masks_settings["masks_folder"], 'masks_shp_{}.shp'.format(stream))):
            missing_masks += [stream]
    attempt_no = 0

    while len(missing_masks) > 0 and attempt_no <= max_attempts:
        logging.info(str(len(missing_masks)) + " masks are missing! Compute...")
        logging.info(', '.join([str(i) for i in missing_masks]))
        logging.info("Attempt " + str(attempt_no))
        exec_pool = get_context('spawn').Pool(process_max)
        for stream in missing_masks:
            exec_pool.apply_async(create_masks, args=(stream, masks_settings, d))
        exec_pool.close()
        exec_pool.join()

        missing_masks_out = deepcopy(missing_masks)
        for stream in missing_masks:
            if os.path.isfile(os.path.join(masks_settings["masks_folder"], 'masks_shp_{}.shp'.format(stream))):
                missing_masks_out.remove(stream)
        missing_masks = missing_masks_out

        attempt_no = attempt_no + 1

    if len(missing_masks) > 0:
        logging.error(" --> Some masks has not been produced after " + str(
            attempt_no) + " attempts! Verify problems in streams: " + ",".join([str(i) for i in missing_masks]))
        raise RuntimeError
    # ------------------------------

    # -------------------------------
    logging.info(' --> Calculate stream statistics...')

    rst_in = {}
    rst_in["slope"] = os.path.join(grass_step1_db, 'rst', 'r_' + domain_name + "_" + rrs + "_dem_slope.tif")
    rst_in["channel"] = os.path.join(grass_step1_db, 'rst', 'r_' + domain_name + "_" + rrs + "_" + drain_method_streams.lower() + "_streams.tif")
    rst_in["flowacc_skm"] = os.path.join(grass_step1_db, 'rst', 'r_' + domain_name + "_" + rrs + "_" + drain_method_streams.lower() + "_flow_acc_skm.tif")
    rst_in["dem"] = os.path.join(grass_step0_db, 'r_' + domain_name + "_" + rrs + "_filled_dem.tif")
    rst_in["slope_channel"] = os.path.join(grass_step1_db, 'rst', 'r_' + domain_name + "_" + rrs + "_slope_channel.tif")

    input_data = {}
    input_data["masks_folder"] = masks_settings["masks_folder"]
    input_data["maps_in"] = rst_in

    conc_time_in = [auth for auth in c_time.keys() if c_time[auth] is True]
    if len(conc_time_in) == 0:
        logging.error(" --> ERROR! Choose at least a concentration time formula!")
        raise ValueError

    input_data["conc_time_in"] = conc_time_in
    streams_gdf_out = d["streams_gdf"].copy()
    time_df_out = pd.DataFrame(index=streams_gdf_out["stream"].values, columns=conc_time_in)
    results = []

    if process_max > 1:
        for chunk in chunks:
            logging.info(" ---> Launching chunk from " + str(min(chunk)) + " to " + str(max(chunk)))
            exec_pool = get_context('spawn').Pool(process_max)
            for stream in chunk:
                results.append(exec_pool.apply_async(calculate_basins_stats, args=(stream, input_data, d)))
            exec_pool.close()
            exec_pool.join()
        logging.info(" --> Collecting output..")
        for result in results:
            res = result.get()
            streams_gdf_out.loc[streams_gdf_out["stream"] == res[0], ["tconc", "tpeak", "treces", "flowAcc_skm"]] = res[1:5]
            time_df_out.loc[time_df_out.index == res[0], conc_time_in] = res[5:]
        logging.info(" --> Collecting output..DONE")
    else:
        for chunk in chunks:
            logging.info(" ---> Launching chunk from " + str(min(chunk)) + " to " + str(max(chunk)))
            for stream in chunk:
                streams_gdf_out.loc[streams_gdf_out["stream"] == stream, ["tconc", "tpeak", "treces",
                                                                          "flowAcc_skm"]] = calculate_basins_stats(
                    stream, input_data, d)[1:5]

    streams_gdf_out.to_file(os.path.join(grass_output_db_vct, 'v_' + domain_name + '_' + rrs + '_streams_features.shp'), driver='ESRI Shapefile')
    time_df_out.to_csv(os.path.join(grass_output_db_txt, 'tab_' + domain_name + '_corr_time_estimation.csv'))
    # -------------------------------

    # ------------------------------------------
    # Estimate total execution time
    tot_elapsed_time_sec = float(time.time() - code_start_time)
    tot_elapsed_time, time_units = Give_Elapsed_Time(tot_elapsed_time_sec)

    # -------------------------------------------------------------------------------------
    # Info algorithm

    logging.info(' ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME ELAPSED: ' + str(tot_elapsed_time) + time_units)
    logging.info(' ==> ... END')
    logging.info(' ==> Bye, Bye')
    logging.info(' ============================================================================ ')
    # -------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-log_file', action="store", dest="alg_log")
    parser_handle.add_argument('-settings_file', action="store", dest="alg_settings")
    parser_handle.add_argument('-base_path', action="store", dest="alg_output")

    parser_values = parser_handle.parse_args()

    if parser_values.alg_settings:
        alg_settings = parser_values.alg_settings
    else:
        alg_settings = 'configuration.json'

    if parser_values.alg_log:
        alg_log = parser_values.alg_log
    else:
        alg_log = 'log.txt'

    if parser_values.alg_output:
        alg_output = parser_values.alg_output
    else:
        alg_output = None

    return alg_settings, alg_log, alg_output


# -------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == "__main__":
    sys.stdout.flush()
    set_start_method('spawn', force=True)
    main()
# ------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------