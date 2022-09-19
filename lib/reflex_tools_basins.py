"""
Library Features:

Name:          reflex_tools_basins
Author(s):     Mauro Arcorace (mauro.arcorace@cimafoundation.org)
               Alessandro Masoero (alessandro.masoero@cimafoundation.org)
               Valerio Basso
               Giulia Bruno (giulia.bruno@cimafoundation.org)
               Alessia Matanò
               Andrea Libertino (andrea.libertino@cimafoundation.org)
Date:          '20220916'
Version:       '2.0.3'
"""
########################################################################################################################
# Import libraries
from pysheds.grid import Grid
import fiona
from copy import deepcopy
import fiona.crs
import geopandas as gpd
import numpy as np
import os, sys
import pandas as pd
import warnings
import statistics
import rioxarray as rxr
from lib.reflex_tools_pfafstetter import Pfaf_find_upstream, Pfaf_find_Nextupstreams, Pfaf_find_downstream
import logging
from osgeo import gdal
import rioxarray as rxr
import xrspatial
import xarray as xr
########################################################################################################################

# Function for create basin masks with and without coastal expansion
def create_masks(stream, masks_settings, d):
    sys.stdout.flush()
    print(' ---> Create mask for stream ' + str(stream))
    streams_gdf = d["streams_gdf"]
    basin_gdf = d["basins_gdf"]
    head_basins = []

    row = streams_gdf.loc[streams_gdf["stream"] == stream]
    pfaf = row['pfA'].values[0]
    iMBID = row['MBID_new'].values[0]
    streams_gdf_mbid = streams_gdf.loc[streams_gdf['MBID_new'] == iMBID].copy()

    # find all upstreams
    iUpstream = Pfaf_find_upstream(pfaf, masks_settings["pfaf_limits"][iMBID][1])
    subset_Upstream = streams_gdf_mbid.loc[(streams_gdf_mbid['pfA'] >= pfaf) & (streams_gdf_mbid['pfA'] <= iUpstream)]


    # find up to 2 next upstreams
    tmp_Upstream_new, head_basins = Pfaf_find_Nextupstreams(pfaf, subset_Upstream, row, stream, head_basins, streams_gdf_mbid)

    # find 2 downstream pfaf(s)
    DownStreamLow, DownStreamUp = Pfaf_find_downstream(pfaf, masks_settings["pfaf_limits"][iMBID][0])

    DownStreamLow = streams_gdf_mbid.loc[(streams_gdf_mbid['pfA'] >= DownStreamLow) & (streams_gdf_mbid['pfA'] < DownStreamUp)]
    DownStreamUp = streams_gdf_mbid.loc[(streams_gdf_mbid['pfA'] >= DownStreamUp) & (streams_gdf_mbid['pfA'] < pfaf)]

    DownStreamLow_sel = DownStreamLow.loc[DownStreamLow['pfA'] == np.max(DownStreamLow['pfA'])]
    DownStreamUp_sel = DownStreamUp.loc[DownStreamUp['pfA'] == np.min(DownStreamUp['pfA'])]

    # list of pfaf(s) of interest
    streams = {}
    streams["all_up"] = masks_settings["upstream_basins"][stream]   #[i for i in subset_Upstream["stream"].values]
    streams["up"] = [i for i in tmp_Upstream_new["stream"].values]
    streams["down"] = [i for i in pd.concat((DownStreamLow_sel, DownStreamUp_sel))["stream"] if i not in streams["up"]]

    df_to_merge = {}

    for mask_type in ['basin', 'ext_no_ce', 'ext_ce']:

        if mask_type == 'ext_no_ce':
            data = streams["up"]
        elif mask_type == 'ext_ce':
            data = streams["up"] + streams["down"]
        #elif mask_type == 'upst_basin':
        #    data = streams["all_up"]
        elif mask_type == 'basin':
            data = streams["all_up"]
            # data = [stream]

        tmp_polygons = []
        for stream_code in data:
            tmp_polygons.append(basin_gdf.loc[(basin_gdf['stream'] == stream_code)])

        df = gpd.GeoDataFrame(pd.concat(tmp_polygons, ignore_index=True))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            df.crs = fiona.crs.from_epsg(int(masks_settings["input_epsg"]))

        temp = df.copy(deep=True)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if mask_type == 'ext_ce':
                temp['geometry'] = temp.buffer(masks_settings["buffer_distance"])
            else:
                temp['geometry'] = temp.buffer(0.000000001)

        # dissolve all polygons
        #if not mask_type == 'basin':
        temp['dissolvefi'] = 1
        df_to_merge[mask_type] = temp.dissolve(by='dissolvefi')
        #else:
        #    df_to_merge[mask_type] = temp

        df_to_merge[mask_type]["mask_type"] = mask_type
        df_to_merge[mask_type]['lat_min'] = df_to_merge[mask_type].bounds.miny.values
        df_to_merge[mask_type]['lat_max'] = df_to_merge[mask_type].bounds.maxy.values
        df_to_merge[mask_type]['lon_min'] = df_to_merge[mask_type].bounds.minx.values
        df_to_merge[mask_type]['lon_max'] = df_to_merge[mask_type].bounds.maxx.values

    df_out = pd.concat((df_to_merge['basin'], df_to_merge['ext_no_ce'], df_to_merge['ext_ce']))
    outfile = os.path.join(masks_settings["masks_folder"], 'masks_shp_{}.shp'.format(stream))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df_out.to_file(outfile)

    return

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to calculate basin main features
def calculate_basins_stats(stream, input_data, d):
    sys.stdout.flush()
    print(" ---> Analyse basin " + str(stream))

    # Import and clip inputs
    mask_FID_shp_i = os.path.join(input_data["masks_folder"], "masks_shp_" + str(stream) + ".shp")
    mask = gpd.read_file(mask_FID_shp_i)
    mask_basin = mask.loc[mask["mask_type"] == 'basin']
    maps_basin = {}
    for map in input_data["maps_in"].keys():
        maps_basin[map] = rxr.open_rasterio(input_data["maps_in"][map], masked=True).rio.clip(mask_basin.geometry, from_disk=True)

    # compute basin features derived from static data
    avg_basSlope = np.nanmean(maps_basin["slope"].values)
    h_basAvg = np.nanmean(maps_basin["dem"].values)
    avg_strSlope = np.nanmean(maps_basin["slope_channel"].values)
    h_basMin = np.nanpercentile(maps_basin["dem"].values, 2)
    h_basMax = np.nanpercentile(maps_basin["dem"].values, 98)
    flow_accum_skm = np.nanmax(maps_basin["flowacc_skm"].values)

    str_len_km = d["streams_gdf"].loc[d["streams_gdf"]["stream"] == stream, "cum_length"].values[0]*(10**(-3))
    #bas_area_km2 = np.nanmax(maps_basin["flowacc_skm"].values) #d["streams_gdf"].loc[d["streams_gdf"]["stream"] == stream, "area_km2"].values

    conc_time_in = input_data["conc_time_in"]
    conc_time = {}

    if 'kirpich' in conc_time_in:
        # Kirpich
        # stream_length[m] - basinslope [m/m]
        conc_time['kirpich'] = ((0.000325 * np.power(str_len_km * 1000, 0.77) * np.power(avg_basSlope / 100.0, -0.385)) * 3600)  # in seconds!!!

    if 'california_culvert_practice' in conc_time_in:
        # California Culvert Practice
        trav_len_miles = str_len_km * 0.621371
        H_diff_feet = (h_basMax - h_basMin) * 3.28084
        conc_time['california_culvert_practice'] = (3600 * (np.power((11.9 * np.power(trav_len_miles, 3) / H_diff_feet), 0.385)))  # in seconds!!!

    if 'ventura' in conc_time_in:
        # Ventura
        # area [km2] - stream_slope [m/m]
        conc_time['ventura'] = (0.1272 * np.sqrt(flow_accum_skm / (avg_strSlope / 100.0)) * 3600)  # in seconds!!!

    if 'pezzoli' in conc_time_in:
        # Pezzoli
        # StreamSlope [m/m] - StreamLength [Km]
        conc_time['pezzoli'] = ((0.055 * str_len_km / np.sqrt(avg_strSlope / 100.0)) * 3600)  # in seconds!!!

    if 'giandotti' in conc_time_in:
        # Giandotti
        # StreamSlope [m/m] - StreamLength [Km]
        conc_time['giandotti'] = 3600 * (((4 * np.sqrt(flow_accum_skm)) + (1.5 * str_len_km)) / (0.8 * np.sqrt((h_basAvg - h_basMin))))  # in seconds!!!

    if 'pasini' in conc_time_in:
        # Pasini
        # StreamSlope [m/m] - StreamLength [Km]
        conc_time['pasini'] = ((0.108 * (np.power((str_len_km * flow_accum_skm), 0.3333)) / (np.power((avg_strSlope / 100.0), 0.5))) * 3600)  # in seconds!!!

    if 'siccardi' in conc_time_in:
        # Siccardi
        # area [km2]
        conc_time['siccardi'] = (3600 * (0.27 * np.sqrt(flow_accum_skm) + 0.25))  # in seconds!!!

        # compute median tc

        # conc_time['tc_avg'] = conc_time[['tc_P', 'tc_K', 'tc_CCP', 'tc_V', 'tc_G', 'tc_Pas', 'tc_S']].mean(axis=1)
        # conc_time['tc_min'] = conc_time[['tc_P', 'tc_K', 'tc_CCP', 'tc_V', 'tc_G', 'tc_Pas', 'tc_S']].min(axis=1)

        # compute centered mean (e.g. excluding min and max of the series
        # sum_tc_tmp = conc_time[['tc_P', 'tc_K', 'tc_CCP', 'tc_V', 'tc_G', 'tc_Pas', 'tc_S']].sum(axis=1)
        # max_tc_tmp = conc_time[['tc_P', 'tc_K', 'tc_CCP', 'tc_V', 'tc_G', 'tc_Pas', 'tc_S']].max(axis=1)
        # min_tc_tmp = conc_time[['tc_P', 'tc_K', 'tc_CCP', 'tc_V', 'tc_G', 'tc_Pas', 'tc_S']].min(axis=1)
        # length = len(['tc_P', 'tc_K', 'tc_CCP', 'tc_V', 'tc_G', 'tc_Pas', 'tc_S'])

        # conc_time['tc_cent_avg'] = (sum_tc_tmp - max_tc_tmp - min_tc_tmp) / (length - 2)

    conc_time['out_value'] = statistics.median([conc_time[i] for i in conc_time_in])

    # If corrivation time degenerate to +Inf (difference of elevation in the basin is neglectible) use Viparelli
    if np.isinf(conc_time['out_value']):
        conc_time['out_value'] = 3600 * str_len_km / 5.4

    # Time to peak and recession time is computed according to
    # https: // www.wcc.nrcs.usda.gov / ftpref / wntsc / H & H / NEHhydrology / ch16.pdf

    # compute peak time
    conc_time['tpeak'] = ((0.133 * conc_time['out_value']) / 2 + 0.6 * conc_time['out_value'])

    # compute recession time
    conc_time['treces'] = (1.67 * conc_time['tpeak'])

    out_conc_time = np.array([conc_time[i] for i in conc_time_in])

    return np.concatenate((np.array((stream, conc_time['out_value'], conc_time['tpeak'], conc_time['treces'], flow_accum_skm)), out_conc_time))
    # np.array((stream, conc_time['tc_med'].values[0], conc_time['tpeak'].values[0], conc_time['treces'].values[0], flow_accum_skm))

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to calculate hand map
def compute_hand(stream_id, hand_settings, d):
    sys.stdout.flush()
    print(" --> Compute stream " + str(stream_id))
    basin_mask = gpd.read_file(hand_settings["input_files"]["shp_masks"].format(stream=str(stream_id)))

    mask_type = "ext_no_ce"
    hand_out_file = os.path.join(hand_settings["output_folder"], mask_type, "hand_" + hand_settings["hand_method"] + "_" + mask_type + "_" + str(stream_id) + ".tif")
    stream_out_file = os.path.join(hand_settings["output_folder"], mask_type, "stream_" + mask_type + "_" + str(stream_id) + ".tif")

    logging.info(" ---> Compute HANDS map")

    # Import data and mask with proper mask file
    basin_mask_shp = basin_mask.loc[basin_mask["mask_type"] == mask_type]
    basin_bbox = tuple(basin_mask_shp[["lon_min", "lat_min", "lon_max", "lat_max"]].values.squeeze())
    grid = deepcopy(hand_settings["grid"])
    #    Grid.from_raster(hand_settings["input_files"]["rst_dem"])
    dem = grid.read_raster(hand_settings["input_files"]["rst_dem"], window = basin_bbox, window_crs = grid.crs)
    basin_mask_raster = grid.rasterize(basin_mask_shp.geometry.values, fill=np.nan)
    grid.clip_to(basin_mask_raster)

    fdir = grid.read_raster(hand_settings["input_files"]["rst_pnt"], window = basin_bbox, window_crs = grid.crs, nodata=-1)
    stream = grid.read_raster(hand_settings["input_files"]["rst_stream"], window = basin_bbox, window_crs = grid.crs, nodata=255)

    if hand_settings["hand_method"] == 'dinf':
        hand = grid.compute_hand(fdir=fdir, dem=dem*100, mask=stream==stream_id, routing='dinf')
    elif hand_settings["hand_method"] == 'd8':
        hand = grid.compute_hand(fdir=fdir, dem=dem * 100, mask=stream == stream_id, dirmap=(2,1,8,7,6,5,4,3), routing='d8')

    grid.to_raster(hand, hand_out_file, apply_output_mask=True, dtype=np.float32, tiled=True)
    stream.mask=stream==stream_id
    grid.to_raster(stream, stream_out_file, apply_input_mask=True, apply_output_mask=True, dtype=np.int16, tiled=True)

    print(" --> Compute stream " + str(stream_id) + " ...DONE!")
    return

# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Function to calculate hand map
def compute_hand_ce(stream_id, hand_settings, d):
    sys.stdout.flush()
    streams_gdf = d["streams_gdf"]
    stream_row = streams_gdf[streams_gdf["stream"] == stream_id]

    mask_type = "ext_ce"
    # if the stream is an outlet or the slope below a slope limit
    if stream_row.next_strea.values[0] == -1 or stream_row.gradient.values[0] <= hand_settings["coastal_expansion_gradient_limit"]:
        print(" --> Compute stream " + str(stream_id))
        # Declare file names
        hand_out_file = os.path.join(hand_settings["output_folder"], mask_type, "hand_" + hand_settings["hand_method"] + "_" + mask_type + "_" + str(stream_id) + ".tif")
        stream_out_file = os.path.join(hand_settings["output_folder"], mask_type, "stream_" + mask_type + "_" + str(stream_id) + ".tif")

        # Import extended mask and use if for import clipped DEM and stream
        basin_mask = gpd.read_file(hand_settings["input_files"]["shp_masks"].format(stream=str(stream_id)))
        basin_mask_shp = basin_mask.loc[basin_mask["mask_type"] == mask_type]

        dem = rxr.open_rasterio(hand_settings["input_files"]["rst_dem"], masked=True).rio.clip(basin_mask_shp.geometry, from_disk=True)
        stream = rxr.open_rasterio(hand_settings["input_files"]["rst_stream"], masked=True).rio.clip(basin_mask_shp.geometry, from_disk=True)
        hand = rxr.open_rasterio(os.path.join(hand_settings["output_folder"], "ext_no_ce", "hand_" + hand_settings["hand_method"] + "_ext_no_ce_" + str(stream_id) + ".tif"), masked=True).rio.clip(basin_mask_shp.geometry, from_disk=True)

        # Calculate cost matrix: distance from current stream limited to the distance set in the stream file
        # First option buffer only up to the limit distance from the reference basin
        # cost = xrspatial.proximity(raster = stream.squeeze(), target_values = [stream_id], distance_metric='GREAT_CIRCLE', max_distance=stream_row["distance"].values[0])
        cost = xrspatial.proximity(raster=stream.squeeze(), target_values=[stream_id], distance_metric='GREAT_CIRCLE')
        cost = cost.reindex({"x": dem.x.values, "y": dem.y.values}, method='nearest', fill_value=np.nan, tolerance=np.abs(hand.rio.resolution()[0]))

        # Calculate elevation percentile
        z_p = np.nanpercentile(dem.values[stream.values==stream_id], hand_settings["z_percentile"])

        # Compute extended DH: DEM - z_p + loss_matrix
        hand_ext = 100 * dem.squeeze() - (100 * z_p) + cost * hand_settings["head_loss"]

        hand_ext = xr.where(hand_ext<0,0,hand_ext)
        if np.min(hand_ext) < 0:
            min_hand_ext = np.min(hand_ext)
        else:
            min_hand_ext = 0
        hand_ext = xr.where(hand_ext<0,hand_ext-min_hand_ext,hand_ext)

        # Regrid hand on extended domain and merge the 2 maps
        hand_ce = hand.reindex({"x": dem.x.values, "y": dem.y.values}, method='nearest', fill_value=np.nan, tolerance=np.abs(hand.rio.resolution()[0])).squeeze()

        # Raw merge
        hand_ce = xr.where(hand_ce.isnull(), hand_ext, hand_ce)

        # Write output
        xr.where(dem.squeeze().isnull(),np.nan,hand_ce).rio.write_nodata(-9999, inplace=True).rio.to_raster(hand_out_file)
        cost.values[cost.values>0]=np.nan
        cost = cost+1
        cost.rio.write_nodata(-9999, inplace=True).rio.to_raster(stream_out_file)

        print(" --> Compute stream " + str(stream_id) + " ...DONE!")
    else:
        print(" --> Stream " + str(stream_id) + "... SKIPPED!")

    return
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
