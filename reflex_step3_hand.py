#!/usr/bin/env python
"""
REFlEx - Step3 - HAND map computation
__date__ = '20230101'
__version__ = '2.0.3'
__author__ =
        'Mauro Arcorace' (mauro.arcorace@cimafoundation.org',
        'Alessandro Masoero (alessandro.masoero@cimafoundation.org)',
        'Valerio Basso',
        'Alessia Matanò',
        'Andrea Libertino (andrea.libertino@cimafoundation.org)'
__library__ = 'REFlEx'
General command line:
### python reflex_step3_hand.py -log_file "/path/to/log.txt" -settings_file "settings.json" -base_path "/path/to/base_folder"
Version(s):
20190220 (1.0.0) --> Beta release
20220406 (2.0.0) --> Full revision - Use pysheds for HAND definition
                                     Parallel implementation
20220808 (2.0.2) --> Optimized multiprocessing
20220808 (2.0.3) --> Optimized singleprocessing for big domains
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Import full libraries
import logging
from argparse import ArgumentParser
import time
import os, sys
import geopandas as gpd
from multiprocessing import Pool, Manager, set_start_method, get_context, cpu_count
from lib.reflex_tools_basins import compute_hand, compute_hand_ce
from lib.reflex_tools_utils import Give_Elapsed_Time, set_logging, read_file_json
from numba import config
from pysheds.grid import Grid
from copy import deepcopy

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'REFlEx - STEP 3 - Hand Map Computation'
alg_version = '2.0.3'
alg_release = '2023-01-01'
# Algorithm parameter(s)
time_format = '%Y%m%d%H%M'


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

    # Previous steps settings
    drain_method_streams = data_settings["step_1"]["stream_definition_method"]
    drain_method_hand = data_settings["step_1"]["hand_definition_method"]
    buffer_watershed_distance_cell = data_settings["step_2"]["buffer_for_coastal_expansion_cells"]

    # Step settings
    step0_dir_name = data_settings["step_0"]["dir_name"].replace("{base_path}", base_path)
    step1_dir_name = data_settings["step_1"]["dir_name"].replace("{base_path}", base_path)
    step2_dir_name = data_settings["step_2"]["dir_name"].replace("{base_path}", base_path)
    step3_dir_name = data_settings["step_3"]["dir_name"].replace("{base_path}", base_path)
    try:
        max_attempts = data_settings["step_3"]["max_attempts_make_hand"]
    except:
        max_attempts = 5

    coastal_expansion_active = data_settings["step_3"]["coastal_expansion"]["enable"]
    gradient_limit = data_settings["step_3"]["coastal_expansion"]["gradient_limit"]
    JL = data_settings["step_3"]["coastal_expansion"]["head_loss_cm_m"]
    in_str_elev_perc_value = data_settings["step_3"]["coastal_expansion"]["slope_percentile"]
    try:
        produce_only_used_hand_ce = data_settings["step_3"]["coastal_expansion"]["produce_only_used_hand_ce"]
    except:
        produce_only_used_hand_ce = True
    try:
        manual_ce = data_settings["step_3"]["coastal_expansion"]["manual_activation"]
    except:
        manual_ce = False

    if data_settings["step_3"]["multiprocessing"]["enable"]:
        process_max = data_settings["step_3"]["multiprocessing"]["max_cores"]
        if process_max is None:
            process_max = cpu_count() - 1
    else:
        process_max = 1

    chunk_size = data_settings["step_3"]["multiprocessing"]["chunk_size"]

    # Set numba execution
    config.THREADING_LAYER_PRIORITY = ["tbb", "omp", "workqueue"]
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Record starting simulation time
    code_start_time = time.time()

    output_folder = os.path.join(step3_dir_name, 'rst', "")
    for mask_type in ["ext_no_ce", "ext_ce"]:
        output_folder_type = os.path.join(output_folder, mask_type, "")
        os.makedirs(output_folder_type, exist_ok=True)

    input_dict = {"domain": domain_name, "res": rrs, "dir_method": drain_method_hand,
                  "dir_method_stream": drain_method_streams}
    input_files_raw = {
        "shp_streams": os.path.join(step2_dir_name, 'vct', 'v_{domain}_{res}_streams_features.shp'),
        "shp_masks": os.path.join(step2_dir_name, 'vct', 'masks', 'masks_shp_{stream}.shp'),
        "rst_pnt": os.path.join(step1_dir_name, 'rst', 'r_{domain}_{res}_{dir_method}_flow_dir.tif'),
        "rst_stream": os.path.join(step1_dir_name, 'rst', 'r_{domain}_{res}_{dir_method_stream}_streams.tif'),
        "rst_dem": os.path.join(step1_dir_name, 'rst', 'r_{domain}_{res}_dem.tif')
    }
    input_files = fill_input_files(input_files_raw, input_dict)
    streams_gdf = gpd.read_file(input_files["shp_streams"])

    if manual_ce is True and "manual_ce" not in streams_gdf.columns:
        raise ValueError(
            "ERROR! In coastal expansion manual mode a column named 'manual_ce' should be included in the stream shapefile!")

    hand_settings = {"input_files": input_files, "hand_method": drain_method_hand, "output_folder": output_folder,
                     "coastal_expansion_gradient_limit": gradient_limit, "head_loss": JL,
                     "z_percentile": in_str_elev_perc_value,
                     "buffer_n_cells": buffer_watershed_distance_cell}

    manager = Manager()
    d = manager.dict()
    d["streams_gdf"] = streams_gdf
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Compute hands no ce
    logging.info("--> Compute hand maps...")
    grid = Grid.from_raster(hand_settings["input_files"]["rst_dem"])
    hand_settings["grid"] = grid
    chunks = [streams_gdf["stream"].values[i:i + chunk_size] for i in
              range(0, len(streams_gdf["stream"].values), chunk_size)]
    if process_max > 1:
        for chunk in chunks:
            logging.info(" ---> Launching chunk from " + str(min(chunk)) + " to " + str(max(chunk)))
            exec_pool = get_context('spawn').Pool(process_max)
            for stream in chunk:
                exec_pool.apply_async(compute_hand, args=(stream, hand_settings, d))
            exec_pool.close()
            exec_pool.join()
        logging.info("--> Compute hand maps...DONE")
    else:
        for chunk in chunks:
            logging.info(" ---> Launching chunk from " + str(min(chunk)) + " to " + str(max(chunk)))
            for stream in chunk:
                compute_hand(stream, hand_settings, d)

    logging.info("--> Verify presence of missing hand maps...")
    missing_hand = []
    attempt_no = 0
    mask_type = "ext_no_ce"
    missing_hand = check_missing_hands(streams_gdf["stream"].values, missing_hand, "ext_no_ce", hand_settings)

    if len(missing_hand) == 0:
        logging.info(" --> All the maps have been correctly produced...")
    else:
        while len(missing_hand) > 0 and attempt_no <= max_attempts:
            logging.info(str(len(missing_hand)) + " hand maps are missing! Compute...")
            logging.info(', '.join([str(i) for i in missing_hand]))
            exec_pool = get_context('spawn').Pool(process_max)
            for stream in missing_hand:
                exec_pool.apply_async(compute_hand, args=(stream, hand_settings, d))
            exec_pool.close()
            exec_pool.join()

            missing_hand_out = deepcopy(missing_hand)
            for stream in missing_hand:
                if os.path.isfile(os.path.join(hand_settings["output_folder"], mask_type,
                                               "hand_" + hand_settings["hand_method"] + "_" + mask_type + "_" + str(
                                                       stream) + ".tif")):
                    missing_hand_out.remove(stream)
            missing_hand = missing_hand_out

            attempt_no = attempt_no + 1

        if len(missing_hand) > 0:
            logging.error(" --> Some hand maps has not been produced after " + str(
                attempt_no) + " attempts! Verify problems in streams: " + ', '.join([str(i) for i in missing_hand]))
            raise RuntimeError
    logging.info("--> Verify presence of missing hand maps...DONE")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Compute hands ce
    if coastal_expansion_active:
        if produce_only_used_hand_ce:
            if not manual_ce:
                streams_in = streams_gdf.loc[
                    (streams_gdf["next_strea"] <= 0) | (streams_gdf["gradient"] <= gradient_limit)]
            else:
                streams_in = streams_gdf.loc[streams_gdf["manual_ce"] == 1]
        else:
            streams_in = streams_gdf
        logging.info("--> Compute hand maps with coastal expansion...")
        chunks = [streams_in["stream"].values[i:i + chunk_size] for i in
                  range(0, len(streams_in["stream"].values), chunk_size)]
        if process_max > 1:
            for chunk in chunks:
                logging.info(" ---> Launching chunk from " + str(min(chunk)) + " to " + str(max(chunk)))
                exec_pool = get_context('spawn').Pool(process_max)
                for stream in chunk:
                    exec_pool.apply_async(compute_hand_ce, args=(stream, hand_settings, d))
                exec_pool.close()
                exec_pool.join()
                logging.info("--> Compute hand maps with coastal expansion...DONE")
        else:
            for chunk in chunks:
                logging.info(" ---> Launching chunk from " + str(min(chunk)) + " to " + str(max(chunk)))
                for stream in chunk:
                    compute_hand_ce(stream, hand_settings, d)

        logging.info("--> Verify presence of missing hand maps...")
        missing_hand = []
        attempt_no = 0
        mask_type = "ext_ce"
        missing_hand = check_missing_hands(streams_in["stream"].values, missing_hand, mask_type, hand_settings)

        if len(missing_hand) == 0:
            logging.info(" --> All the maps have been correctly produced...")
        else:
            while len(missing_hand) > 0 and attempt_no <= max_attempts:
                logging.info(str(len(missing_hand)) + " masks are missing! Compute...")
                logging.info(', '.join([str(i) for i in missing_hand]))
                exec_pool = get_context('spawn').Pool(process_max)
                for stream in missing_hand:
                    exec_pool.apply_async(compute_hand_ce, args=(stream, hand_settings, d))
                exec_pool.close()
                exec_pool.join()

                missing_hand_out = deepcopy(missing_hand)
                for stream in missing_hand:
                    if os.path.isfile(os.path.join(hand_settings["output_folder"], mask_type,
                                                   "hand_" + hand_settings["hand_method"] + "_" + mask_type + "_" + str(
                                                           stream) + ".tif")):
                        missing_hand_out.remove(stream)
                missing_hand = missing_hand_out

                attempt_no = attempt_no + 1
            if len(missing_hand) > 0:
                logging.error(" --> Some hand maps has not been produced after " + str(
                    attempt_no) + " attempts! Verify problems in streams: " + ",".join([str(i) for i in missing_hand]))
                raise RuntimeError
        logging.info("--> Verify presence of missing hand maps...DONE")
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Estimate total execution time
    tot_elapsed_time_sec = float(time.time() - code_start_time)
    tot_elapsed_time, time_units = Give_Elapsed_Time(tot_elapsed_time_sec)

    # Info algorithm

    logging.info(' ')
    logging.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    logging.info(' ==> TIME ELAPSED: ' + str(tot_elapsed_time) + time_units)
    logging.info(' ==> ... END')
    logging.info(' ==> Bye, Bye')
    logging.info(' ============================================================================ ')
    # -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
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

def fill_input_files(input_files_raw, input_dict):
    filled_input = {}
    input_dict["stream"] = "{stream}"
    for key in input_files_raw.keys():
        filled_input[key] = input_files_raw[key].format(**input_dict)
    return filled_input


# -------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def check_missing_hands(streams, missing_hand, mask_type, hand_settings):
    for stream in streams:
        if not os.path.isfile(os.path.join(hand_settings["output_folder"], mask_type,
                                           "hand_" + hand_settings["hand_method"] + "_" + mask_type + "_" + str(
                                               stream) + ".tif")):
            missing_hand += [stream]
    return missing_hand


# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == "__main__":
    sys.stdout.flush()
    set_start_method('spawn', force=True)
    main()
# ----------------------------------------------------------------------------