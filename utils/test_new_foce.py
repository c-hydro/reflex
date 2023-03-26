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

stream_id = 12
mask_type = "no_ce"
hand_settings = {}

hand_settings["output_folder"] = "/home/andrea/Desktop/nile_tests/output/"
hand_settings["hand_method"] = "d8"
hand_settings["input_files"] = {}
hand_settings["input_files"]["shp_masks"] = "/home/andrea/Desktop/nile_tests/mask_basin.shp"
hand_settings["input_files"]["rst_dem"] = "/home/andrea/Desktop/nile_tests/foce_fill.tif"
hand_settings["input_files"]["rst_stream"] = "/home/andrea/Desktop/nile_tests/streams_foce.tif"
hand_settings["input_files"]["rst_pnt"] = "/home/andrea/Desktop/nile_tests/foce_fdir.tif"
hand_settings["z_percentile"] = 0.001
hand_settings["head_loss"] = 0.02


print(" --> Compute stream " + str(stream_id))
basin_mask_shp = gpd.read_file(hand_settings["input_files"]["shp_masks"].format(stream=str(stream_id)))

mask_type = "no_ce"
hand_out_file = os.path.join(hand_settings["output_folder"], mask_type,
                             "hand_" + hand_settings["hand_method"] + "_" + mask_type + "_" + str(stream_id) + ".tif")
stream_out_file = os.path.join(hand_settings["output_folder"], mask_type,
                               "stream_" + mask_type + "_" + str(stream_id) + ".tif")

logging.info(" ---> Compute HANDS map")

# Import data and mask with proper mask file
basin_bbox = tuple([29,32,29,32])
grid = Grid.from_raster(hand_settings["input_files"]["rst_dem"])

#    Grid.from_raster(hand_settings["input_files"]["rst_dem"])
dem = grid.read_raster(hand_settings["input_files"]["rst_dem"], window_crs=grid.crs)
basin_mask_raster = grid.rasterize(basin_mask_shp.geometry.values, fill=np.nan)
grid.clip_to(basin_mask_raster)

fdir = grid.read_raster(hand_settings["input_files"]["rst_pnt"], window_crs=grid.crs, nodata=-1)
stream = grid.read_raster(hand_settings["input_files"]["rst_stream"], window_crs=grid.crs, nodata=255)

if hand_settings["hand_method"] == 'dinf':
    hand = grid.compute_hand(fdir=fdir, dem=dem * 100, mask=stream == stream_id, routing='dinf')
elif hand_settings["hand_method"] == 'd8':
    hand = grid.compute_hand(fdir=fdir, dem=dem * 100, mask=stream == stream_id, dirmap=(2, 1, 8, 7, 6, 5, 4, 3), routing='d8')

grid.to_raster(hand, hand_out_file, apply_output_mask=True, dtype=np.float32, tiled=True)
stream.mask = stream == stream_id
grid.to_raster(stream, stream_out_file, apply_input_mask=True, apply_output_mask=True, dtype=np.int16, tiled=True)

print(" --> Compute stream " + str(stream_id) + " ...DONE!")

mask_type = "ce"
# if the stream is an outlet or the slope below a slope limit

print(" --> Compute stream " + str(stream_id))
# Declare file names
hand_out_file = os.path.join(hand_settings["output_folder"], mask_type, "hand_" + hand_settings["hand_method"] + "_" + mask_type + "_" + str(stream_id) + ".tif")
stream_out_file = os.path.join(hand_settings["output_folder"], mask_type, "stream_" + mask_type + "_" + str(stream_id) + ".tif")

# Import extended mask and use if for import clipped DEM and stream
basin_mask_shp = gpd.read_file(hand_settings["input_files"]["shp_masks"].format(stream=str(stream_id)))
#basin_mask_shp = basin_mask.loc[basin_mask["mask_type"] == mask_type]

dem = rxr.open_rasterio(hand_settings["input_files"]["rst_dem"], masked=True).rio.clip(basin_mask_shp.geometry, from_disk=True)
stream = rxr.open_rasterio(hand_settings["input_files"]["rst_stream"], masked=True).rio.clip(basin_mask_shp.geometry, from_disk=True)
hand = rxr.open_rasterio(os.path.join(hand_settings["output_folder"], "no_ce", "hand_" + hand_settings["hand_method"] + "_no_ce_" + str(stream_id) + ".tif"), masked=True).rio.clip(basin_mask_shp.geometry, from_disk=True)

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