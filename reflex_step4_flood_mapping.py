#!/usr/bin/env python
"""
REFlEx - Step4 - Flood mapping
__date__ = '20220406'
__version__ = '2.0.0'
__author__ =
        'Mauro Arcorace' (mauro.arcorace@cimafoundation.org',
        'Alessandro Masoero (alessandro.masoero@cimafoundation.org',
        'Valerio Basso,
        'Alessia Matanò,
        'Andrea Libertino (andrea.libertino@cimafoundation.org',
        'Lorenzo Alfieri (lorenzo.alfieri @cimafoundation.org'
__library__ = 'REFlEx'
General command line:
### python reflex_step4_flood_mapping.py -log_file "/path/to/log.txt" -settings_file "settings.json" -base_path "/path/to/base_folder" (-volumes_file "/path/to/file.csv" -rp "10")
Version(s):
20190220 (1.0.0) --> Beta release
20220406 (2.0.0) --> Full revision
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------

# Import python libraries
import logging
from argparse import ArgumentParser
import geopandas as gpd
import numpy as np
import os
import sys
import time
import rioxarray as rxr
import xarray as xr
import shutil
import rasterio as rio
import rasterio.merge as rio_merge

import pandas as pd
from lib.reflex_tools_flooding import optimise_volume
from lib.reflex_tools_utils import Give_Elapsed_Time, Start_GRASS_py3, set_logging, read_file_json
from multiprocessing import Pool, Manager, set_start_method, get_context, cpu_count
from numba import config
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'REFlEx - STEP 4 - Flood mapping'
alg_version = '2.0.0'
alg_release = '2022-04-06'
# Algorithm parameter(s)
time_format = '%Y%m%d%H%M'

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Script Main
def main():
    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    config_file, alg_log, base_path, discharge_file_name, return_period = get_args()

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
    coastal_expansion_active = data_settings["step_3"]["coastal_expansion"]["enable"]
    gradient_limit = data_settings["step_3"]["coastal_expansion"]["gradient_limit"]

    # Step settings
    step0_dir_name = data_settings["step_0"]["dir_name"].replace("{base_path}", base_path)
    step1_dir_name = data_settings["step_1"]["dir_name"].replace("{base_path}", base_path)
    step2_dir_name = data_settings["step_2"]["dir_name"].replace("{base_path}", base_path)
    step3_dir_name = data_settings["step_3"]["dir_name"].replace("{base_path}", base_path)
    step4_dir_name = data_settings["step_4"]["dir_name"].replace("{base_path}", base_path)

    runoff_volume_min = data_settings["step_4"]["minimum_volume_included"]
    optimise_setting = data_settings["step_4"]["optimization"]
    roughness_coeff = data_settings["step_4"]["roughness_coefficient"]

    # Ancillary info (can be provided also by command line)
    if return_period is None:
        return_period = data_settings["step_4"]["return_period"]
    if discharge_file_name is None:
        discharge_file_name = data_settings["step_4"]["volume_file_path"]

    step4_dir_name = os.path.join(step4_dir_name, "T" + str(return_period),"")

    if data_settings["step_4"]["multiprocessing"]["enable"]:
        process_max = data_settings["step_3"]["multiprocessing"]["max_cores"]
        if process_max is None:
            process_max = cpu_count() - 1
    else:
        process_max = 1

    max_maps = data_settings["step_4"]["merging_maps"]["max_maps_opened_together"]

    config.THREADING_LAYER_PRIORITY = ["tbb", "omp", "workqueue"]
    ################################################################################

    # Record starting simulation time
    code_start_time = time.time()

    # Make output folders
    os.makedirs(os.path.join(step4_dir_name, 'rst', ''), exist_ok=True)
    os.makedirs(os.path.join(step4_dir_name, 'tmp', ''), exist_ok=True)
    os.makedirs(os.path.join(step4_dir_name, 'txt', ''), exist_ok=True)
    # --------------------------------------------
    # Calculate volume
    stream_gdf = gpd.read_file(os.path.join(step2_dir_name, 'vct', 'v_' + domain_name + '_' + rrs + '_streams_features.shp')).copy()
    volume_table = pd.read_csv(discharge_file_name, sep=',', header=0, index_col='stream')
    stream_gdf = stream_gdf.join(volume_table, on='stream', rsuffix='_temp')

    stream_gdf["V0"] = ((stream_gdf['Q' + str(return_period)] - stream_gdf['Q0']) * (stream_gdf['tconc']) / 2.).round(0)

    stream_gdf["Q_p_t"] = stream_gdf["Q" + str(return_period)] - stream_gdf["Q0"]

    # Calculate angle of triangular hydrograph
    stream_gdf["slope_hydrograph"] = np.arctan2(stream_gdf["Q_p_t"].values,stream_gdf["tconc"].values)

    # Remove streams with volume below a treshold
    stream_gdf["V0"].values[stream_gdf["V0"].values is None] = -1
    stream_gdf["V0"].values[np.isnan(stream_gdf["V0"].values)] = -1
    stream_gdf = stream_gdf.drop(stream_gdf[stream_gdf["V0"].values < runoff_volume_min].index)

    stream_ids = stream_gdf["stream"]

    optimise_setting["roughness_coeff"] = roughness_coeff
    optimise_setting["coastal_expansion_gradient_limit"] = gradient_limit
    optimise_setting["coastal_expansion_active"] = coastal_expansion_active
    optimise_setting["rst_hand"] = os.path.join(step3_dir_name,'rst','{mask_type}','hand_' + drain_method_hand + '_{mask_type}_{stream_id}.tif')
    optimise_setting["rst_stream"] = os.path.join(step3_dir_name, 'rst', '{mask_type}', 'stream_{mask_type}_{stream_id}.tif')
    optimise_setting["out_path"] = step4_dir_name
    optimise_setting["areacell_m"] = (rrs_km * 1000) ** 2

    # Optimization loop
    logging.info(" --> Launch optimization...")
    manager = Manager()
    d = manager.dict()
    d["streams_gdf"] = stream_gdf

    streams_gdf_out = stream_gdf.copy().drop(columns=["geometry"])
    streams_gdf_out["V_cut"] = -9999
    streams_gdf_out["v_i"] = -9999
    streams_gdf_out["diff_final"] = -9999

    results = []
    exec_pool = get_context('spawn').Pool(process_max)
    for stream in stream_ids:
        results.append(exec_pool.apply_async(optimise_volume, args=(stream, optimise_setting, d)))
    exec_pool.close()
    exec_pool.join()

    logging.info(" --> Collecting output..")
    for result in results:
        res = result.get()
        streams_gdf_out.loc[streams_gdf_out["stream"] == res[0], ["V_cut", "v_i", "diff_final"]] = res[1:]
    logging.info(" --> Collecting output...DONE!")

    # Postprocessing
    logging.info(" --> Merging flood maps...")
    grid = rxr.open_rasterio(os.path.join(step1_dir_name,'rst','r_' + domain_name  + "_" + rrs + "_" + drain_method_streams.lower() + "_streams.tif"))
    grid.values = np.zeros(grid.values.shape).astype("float32")
    base_map = os.path.join(optimise_setting["out_path"], "tmp", "base_map.tif")
    grid.rio.to_raster(base_map, compress="DEFLATE")

    src_grid = rio.open(os.path.join(optimise_setting["out_path"], "tmp", "base_map.tif"))
    optimise_setting["meta"] = src_grid.meta.copy()
    optimise_setting["meta"]["compress"] = "DEFLATE"
    chunks = [stream_ids.values[i:i + max_maps] for i in range(0, len(stream_ids.values), max_maps)]

    if data_settings["step_4"]["merging_maps"]["multiprocessing_enabled"]:
        results =[]
        exec_pool = get_context('spawn').Pool(process_max)

        for chunk in chunks:
            results.append(exec_pool.apply_async(merge_maps_mp, args=(chunk, optimise_setting)))
        exec_pool.close()
        exec_pool.join()

        logging.info(" --> Merge final map")
        out_merged_array, out_trans = rio_merge.merge([result.get() for result in results], method=custom_merge)
    else:
        l_0 = [base_map]
        out_merged_array = grid.copy().values
        for chunk in chunks:
            logging.info(" --> Merge maps for chunk " + str(chunk[0]) + "_" + str(chunk[-1]))
            lista = l_0 + [os.path.join(optimise_setting["out_path"], "tmp", "flood_m_bas_" + str(stream_id) + ".tif") for stream_id in chunk]
            chunk_max, out_trans = rio_merge.merge(lista, method=custom_merge)
            out_merged_array = np.maximum(chunk_max, out_merged_array)
            logging.info(" --> Merge maps for chunk " + str(chunk[0]) + "_" + str(chunk[-1]) + "...DONE")

    with rio.open(os.path.join(optimise_setting["out_path"], 'rst', 'flood_map_T' + str(return_period) + '_m.tif'), "w", **optimise_setting["meta"]) as dest:
        dest.write(out_merged_array)

    streams_gdf_out.to_csv(os.path.join(step4_dir_name, 'txt', 'flood_stats_T' + str(return_period) + '_m.csv'))

    shutil.rmtree(os.path.join(optimise_setting["out_path"], "tmp"))
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Estimate total execution time
    tot_elapsed_time_sec = float(time.time() - code_start_time)
    tot_elapsed_time,time_units = Give_Elapsed_Time(tot_elapsed_time_sec)

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
    parser_handle.add_argument('-volumes_file', action="store", dest="alg_volumes")
    parser_handle.add_argument('-rp', action="store", dest="alg_rp")

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

    if parser_values.alg_volumes:
        alg_volumes = parser_values.alg_volumes
    else:
        alg_volumes = None

    if parser_values.alg_rp:
        alg_rp = parser_values.alg_rp
    else:
        alg_rp = None

    return alg_settings, alg_log, alg_output, alg_volumes, alg_rp

# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def fill_input_files(input_files_raw, input_dict):
    filled_input = {}
    input_dict["stream"] = "{stream}"
    for key in input_files_raw.keys():
        filled_input[key] = input_files_raw[key].format(**input_dict)
    return filled_input
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def custom_merge(old_data, new_data, old_nodata, new_nodata, index=None, roff=None, coff=None):
    old_data[:] = np.fmax(old_data, new_data)  # <== NOTE old_data[:] updates the old data array *in place*
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
def merge_maps_mp(chunk, optimise_setting):
    print(" --> Merge maps for chunk " + str(chunk[0]) + "_" + str(chunk[-1]))
    lista = [os.path.join(optimise_setting["out_path"], "tmp", "base_map.tif")] + [os.path.join(optimise_setting["out_path"], "tmp", "flood_m_bas_" + str(stream_id) + ".tif") for stream_id in chunk]
    merged_array, out_trans = rio_merge.merge(lista, method=custom_merge)
    out_name = os.path.join(optimise_setting["out_path"], "tmp", "merged_" + str(chunk[0]) + "_" + str(chunk[-1]) + ".tif")
    with rio.open(out_name, "w", **optimise_setting["meta"]) as dest:
        dest.write(merged_array)
    print(" --> Merge maps for chunk " + str(chunk[0]) + "_" + str(chunk[-1]) + "...DONE")
    return out_name
# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == "__main__":
    sys.stdout.flush()
    set_start_method('spawn', force=True)
    main()
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------