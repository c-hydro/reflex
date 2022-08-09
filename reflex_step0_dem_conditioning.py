#!/usr/bin/env python
"""
REFlEx - Step0 - Hydrogical Conditioning
__date__ = '20220406'
__version__ = '2.0.0'
__author__ =
        'Mauro Arcorace' (mauro.arcorace@cimafoundation.org',
        'Andrea Libertino (andrea.libertino@cimafoundation.org)'
__library__ = 'REFlEx'
General command line:
### python reflex_step0_dem_conditioning.py -log_file "/path/to/log.txt" -settings_file "settings.json" -dem "/path/to/dem.tif" -base_path "/path/to/base_folder" (-stream "/path/to/stream.shp")
Version(s):
20190220 (1.0.0) --> Beta release
20210126 (1.5.0) --> Official release
20220118 (1.6.0) --> Added simplified carving
20220406 (2.0.0) --> Full revision
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Import full libraries
import geopandas as gpd
from lib.reflex_tools_utils import Give_Elapsed_Time, Start_GRASS_py3, set_logging, read_file_json
import logging
from matplotlib import pyplot as plt
import numpy as np
import os
from osgeo import gdal, osr, gdalconst
from pysheds.grid import Grid
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import richdem as rd
from scipy import interpolate
from scipy.signal import savgol_filter
import shutil
import subprocess
import time
from argparse import ArgumentParser
import json

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'REFlEx - STEP 0 - Dem Conditioning'
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
    config_file, alg_log, input_dem_geotiff, input_stream_shapefile, base_path = get_args()

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
    min_threshold = data_settings["domain"]["min_nonnull_value"]

    # Step settings
    step0_dir_name = data_settings["step_0"]["dir_name"].replace("{base_path}", base_path)

    # Preprocessing
    remove_sinks_enabled = data_settings["step_0"]["preprocessing"]["hydrodem_remove_sink"]
    epsilon_filling_enabled = data_settings["step_0"]["preprocessing"]["richdem_eps_filling"]
    pysheds_preprocessing_enabled = data_settings["step_0"]["preprocessing"]["pysheds_pit_remove_and_solve_flats"]

    # Stream burning
    method_dem_burning = data_settings["step_0"]["stream_burning"]["burning_method"]
    if method_dem_burning is None:
        stream_burning_enabled = False
    else:
        stream_burning_enabled = True
    carve_depth_m = data_settings["step_0"]["stream_burning"]["A"]["carve_depth_m"]
    filtering_magnitude = data_settings["step_0"]["stream_burning"]["B"]["filtering_magnitude"]
    buffer_distance_cells = data_settings["step_0"]["stream_burning"]["B"]["buffer_distance_cells"]
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------

    # CREATE OUTPUT DIRECTORIES
    # -------------------------

    # Define output MAIN directory
    # out_folder='%s_%s_%s_%s'%(str(out_folder),str(buffer_cells),str(int(filtering_magnitude)),str(target_epsg))
    main_output_dir = os.path.join(base_path, step0_dir_name)
    logging.info('--> Output main directory: %s' % main_output_dir)
    if not os.path.exists(main_output_dir):
        os.makedirs(main_output_dir)

    # Define output tmp sub-directory
    out_tmp_dirname = 'tmp'
    tmp_output_dir = os.path.join(main_output_dir, out_tmp_dirname)
    logging.info('--> Output tmp sub-directory: %s' % tmp_output_dir)
    if not os.path.exists(tmp_output_dir):
        os.makedirs(tmp_output_dir)

    # Define output tmp sub-directory
    if method_dem_burning == 'B':
        out_segment_dirname = 'segments'
        out_segment_dir = os.path.join(main_output_dir, out_segment_dirname)
        logging.info('--> Output segment sub-directory: %s' % out_segment_dir)
        if not os.path.exists(out_segment_dir):
            os.makedirs(out_segment_dir)

    ################################################################################

    ################################################################################

    # CHECK EPSG CODE
    # -------------------------

    if target_epsg != source_epsg:

        logging.info("Projecting raster to EPSG to %s..." % (str(target_epsg)))
        output_DTM_file_prj = Project_geotiff(input_dem_geotiff, source_epsg, target_epsg, tmp_output_dir)

        logging.info("Projecting vector from EPSG %s to %s..." % (str(source_epsg), str(target_epsg)))
        output_prj_stream_shpfile = Project_shapefile(input_stream_shapefile, source_epsg, target_epsg, tmp_output_dir)

    else:
        logging.info("--> Not required to projecting input raster or vector files.")
        output_DTM_file_prj = input_dem_geotiff
        output_prj_stream_shpfile = input_stream_shapefile

        ################################################################################

    ################################################################################

    # INITIALIZE GRASS
    # ----------------

    # Start GRASS
    logging.info("--> Starting GRASS..")
    grass_temp_db_name = 'grass_temp'
    gisbase, grass_temp_db, location, mapset = Start_GRASS_py3(grass_bin, str(target_epsg), grass_temp_db_name,
                                                               main_output_dir)
    logging.info('--> GRASS temporary database path: %s' % grass_temp_db)
    logging.info('--> GRASS temporary database directory: %s' % location)
    logging.info('--> GRASS mapset: %s' % mapset)

    # Import GRASS Python bindings
    import grass.script as grass
    import grass.script.setup as gsetup

    # Launch session and do something
    gsetup.init(gisbase, grass_temp_db, location, mapset)
    logging.info('--> GRASS GIS environment: %s' % grass.gisenv())

    # Import GRASS modules
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.pygrass.modules.shortcuts import vector as v
    from grass.pygrass.modules.shortcuts import database as db
    from grass.pygrass.modules import Module

    # Define alias for GRASS modules
    r.out_gdal = Module('r.out.gdal')
    r.stream_basins = Module('r.stream.basins')
    r.stream_extract = Module('r.stream.extract')
    r.stream_order = Module('r.stream.order')
    r.stream_distance = Module('r.stream.distance')
    r.grow_distance = Module('r.grow.distance')
    r.to_vect = Module('r.to.vect')
    v.to_db = Module('v.to.db')
    v.db_dropcolumn = Module('v.db.dropcolumn')
    v.db_addcolumn = Module('v.db.addcolumn')
    v.what_rast = Module('v.what.rast')
    v.out_ogr = Module('v.out.ogr')
    v.db_renamecolum = Module('v.db.renamecolumn')
    v.db_select = Module('v.db.select')
    v.db_addtable = Module('v.db.addtable')
    v.in_ogr = Module('v.in.ogr')
    v.db_update = Module('v.db.update')
    v.db_select = Module('v.db.select')
    v.to_rast = Module('v.to.rast')
    db.out_ogr = Module('db.out.ogr')
    v.to_points_new = Module('v.to.points')

    # ----------------

    ################################################################################

    ################################################################################

    # DEFINE VARIABLE NAMES
    # ---------------------

    progress_dem = 'r_%s_%s_dem_in_progress' % (domain_name, str(rrs))

    dem_coverage = 'r_%s_%s_dem_coverage' % (domain_name, str(rrs))
    hydrodem_dem = 'r_%s_%s_dem_hydrodem' % (domain_name, str(rrs))

    conditioned_dem = 'r_%s_%s_conddem' % (domain_name, rrs)
    filled_dem = 'r_%s_%s_filled_dem' % (domain_name, rrs)

    # ---------------------
    code_start_time = time.time()
    ################################################################################

    ################################################################################

    # IMPORT DEM INTO GRASS AND MAKE HILLSHADE
    # ----------------------------------------
    dem_processed = os.path.join(tmp_output_dir, progress_dem)

    # Import DEM file into GRASS
    r.external(input=output_DTM_file_prj, output='temp_in')
    g.region(raster='temp_in')

    r.mapcalc("raw_dem=if(temp_in<" + str(min_threshold) + ",null(),temp_in)", overwrite=overwrite_mode, quiet=quiet_mode)
    r.out_gdal(input="raw_dem", output=dem_processed, type="Float32", createopt="COMPRESS=DEFLATE", nodata=-9999, overwrite=overwrite_mode, quiet=quiet_mode)
    ############################

    # Compute DEM coverage
    r.mapcalc("%s=not(if(isnull(%s)))" % (dem_coverage, 'raw_dem'), overwrite=overwrite_mode, quiet=quiet_mode)

    # ----------------------------------------

    ################################################################################

    ################################################################################

    # FILTERING DEM
    # -------------
    # Smooth input DEM using mdenoise (http://carlosgrohmann.com/blog/denoise-srtm-grass/)
    # /application/pi/Documents/Reflex_model/bin/mdenoise -i /application/pi/Documents/LIDAR/MATTM_LiDAR_MAGRA_resampled_bln_10m_3003_masked_trnslt.asc -o /application/pi/Documents/LIDAR/MATTM_LiDAR_MAGRA_resampled_bln_10m_3003_masked_trnslt_denoise5.asc -z -t 0.9 -n 5

    # Removing sinks by executing hydrodem
    if (remove_sinks_enabled is True):

        logging.info('-> Removing sinks by executing hydrodem...')
        r.hydrodem(flags='a', input='raw_dem', output=hydrodem_dem, mod='4', size='4', overwrite=overwrite_mode,
                   quiet=quiet_mode)
        r.out_gdal(input=hydrodem_dem, output=dem_processed, createopt="COMPRESS=DEFLATE",
                   overwrite=overwrite_mode, quiet=quiet_mode)

    else:
        pass

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # DEM preprocessing

    # Compute RichDEM epsilon filling
    if epsilon_filling_enabled:

        # Read raw dem
        input_E = rd.LoadGDAL(dem_processed, no_data=-9999)
        logging.info('-> Converting numpy to richdem array...')

        # DEM Epsilon filling
        logging.info('-> Epsilon filling raw DEM...')
        rd.FillDepressions(input_E, epsilon=True, in_place=True)
        rd.SaveGDAL(dem_processed, input_E)

    else:
        pass

    # Solve flats using pysheds
    if pysheds_preprocessing_enabled:
        logging.info('-> Pit filling and solving flats with pysheds...')
        grid = Grid.from_raster('%s' % str(dem_processed))
        dem = grid.read_raster('%s' % str(dem_processed), window=grid.bbox, window_crs=grid.crs)
        dem.nodata = -9999
        pit_filled_dem = grid.fill_pits(dem, nodata_out=np.nan)
        inflated_dem = grid.resolve_flats(pit_filled_dem, nodata_out=np.nan)
        grid.to_raster(inflated_dem, dem_processed, nodata=-9999, apply_output_mask=True)
    else:
        pass

    out_dem_filled_file = '%s/%s.tif' % (main_output_dir, filled_dem)
    gdal.Translate(out_dem_filled_file, dem_processed, format="GTiff", creationOptions=["COMPRESS=DEFLATE"], outputType=gdalconst.GDT_Float32)
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # Import filtered DEM raster file and overwrite the raw dem
    r.external(input=out_dem_filled_file, output='raw_dem', overwrite=overwrite_mode, quiet=quiet_mode)
    g.region(raster='raw_dem')
    r.mask(raster=dem_coverage)

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # Stream burning
    if stream_burning_enabled is True:

        v.in_ogr(input=output_prj_stream_shpfile, output='streams_vect', quiet=quiet_mode)

        if method_dem_burning == 'A':
            logging.info('-> Simplified stream burning...')
            v.to_rast(input='streams_vect', type='line', output='streams_rast', use='val', value=carve_depth_m,
                      quiet=quiet_mode)
            r.mapcalc(expression="out_conditioned_dem=if(isnull(streams_rast),raw_dem,raw_dem-streams_rast)",
                      quiet=quiet_mode)

            # Export and define output variable
            r.out_gdal(input='out_conditioned_dem', output='%s/%s.tif' % (main_output_dir, conditioned_dem), type='Float32',
                       overwrite=overwrite_mode, nodata=-9999, createopt="COMPRESS=DEFLATE", quiet=quiet_mode)
            # dem_input_cond_file='%s/%s.tif'%(main_output_dir,conditioned_dem)
            out_dem_cond_file = '%s/%s.tif' % (main_output_dir, conditioned_dem)
            logging.info('-> Simplified stream burning...DONE')

        elif method_dem_burning == 'B':
            # Get stream ID list
            stream_ID_list_orig = []
            p = grass.pipe_command('v.db.select', flags="c", map='streams_vect', columns='stream', quiet=quiet_mode)
            for line in p.stdout:
                stream_ID_list_orig.append(line.decode("utf-8").rstrip('\r\n'))
            p.wait()
            logging.info("--> stream_ID_list_orig: %s" % stream_ID_list_orig)

            # Get stream ID level and depth list from input vector file
            shape_input_stream_network = gpd.read_file(output_prj_stream_shpfile)
            shape_input_stream_network_sorted = shape_input_stream_network.iloc[
                shape_input_stream_network['stream'].sort_values(ascending=[True]).index.values]
            stream_ID_array = shape_input_stream_network_sorted['stream'].values
            stream_ID_level_array = shape_input_stream_network_sorted['Level'].values
            stream_ID_outlet_depth_array = shape_input_stream_network_sorted['Depth'].values
            stream_ID_inlet_depth_array = shape_input_stream_network_sorted['Depth_beg'].values
            stream_ID_list = stream_ID_array.tolist()
            stream_ID_level_list = stream_ID_level_array.tolist()
            stream_ID_outlet_depth_list = stream_ID_outlet_depth_array.tolist()
            stream_ID_inlet_depth_list = stream_ID_inlet_depth_array.tolist()
            del shape_input_stream_network, shape_input_stream_network_sorted

            # Extract maximum area value and identify location of sub-basin outlet for %s streams..
            tot_streams = len(stream_ID_list)
            logging.info(
                "--> Burning streams into DEM for %s branches (this process might take several minutes).." % tot_streams)

            logging.info("--> stream_ID_list: %s" % stream_ID_list)
            logging.info("--> stream_ID_level_list: %s" % stream_ID_level_list)
            logging.info("--> stream_ID_outlet_depth_list: %s" % stream_ID_outlet_depth_list)
            logging.info("--> stream_ID_inlet_depth_list: %s" % stream_ID_inlet_depth_list)

            segm_count = 0

            # for ID in stream_ID_list:
            for k in range(len(stream_ID_list)):

                ID = str(stream_ID_list[k])
                burning_level = int(stream_ID_level_list[k])
                burning_depth = float(stream_ID_outlet_depth_list[k])
                already_burned_depth = float(stream_ID_inlet_depth_list[k])

                logging.info("--> Burning stream no.: %s" % ID)
                logging.info("--> Burning stream no. %s with level %s and carve stream from %s to %s m of depth" % (
                    ID, str(burning_level), str(already_burned_depth), str(burning_depth)))
                flat_burning_mode_flag = False
                if already_burned_depth == burning_depth:
                    flat_burning_mode_flag = True
                else:
                    pass

                g.region(raster='raw_dem')
                g.region(zoom='raw_dem')

                v.extract(input='streams_vect', type='line', where='stream=%s' % (str(ID)),
                          output='stream_i_vect_%s' % str(ID), overwrite=overwrite_mode, quiet=quiet_mode)
                v.out_ogr(input='stream_i_vect_%s' % str(ID),
                          output='%s/stream_i_vect_%s.gpkg' % (tmp_output_dir, str(ID)),
                          format="GPKG", type="line", overwrite=overwrite_mode, quiet=quiet_mode)

                v.to_rast(input='stream_i_vect_%s' % str(ID), output='streams_rast_%s' % str(ID), type='line',
                          use='attr',
                          attribute_column='cat_', overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='streams_rast_%s' % str(ID), output='%s/streams_rast_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)

                r.mapcalc("stream_i_tmp_%s=if(streams_rast_%s==%s,%s,null())" % (str(ID), str(ID), str(ID), str(ID)),
                          overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='stream_i_tmp_%s' % str(ID), output='%s/stream_i_tmp_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)
                r.mapcalc("dem_nan_mask_i_%s=if(!isnull(raw_dem),1,0)" % (str(ID)), overwrite=overwrite_mode,
                          quiet=quiet_mode)
                r.mapcalc("stream_i_%s=if(dem_nan_mask_i_%s==1,stream_i_tmp_%s,null())" % (str(ID), str(ID), str(ID)),
                          overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='stream_i_%s' % str(ID), output='%s/stream_i_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)

                in_region = grass.region()
                logging.info("--> GRASS computational region: '%s'" % in_region)
                try:
                    r.thin(input='stream_i_%s' % str(ID), output='stream_i_thinned_%s' % str(ID))
                except:
                    print('SKIPPED stream_i_%s' % str(ID))
                    continue
                r.out_gdal(input='stream_i_thinned_%s' % str(ID),
                           output='%s/stream_i_thinned_%s.tif' % (tmp_output_dir, ID), createopt="COMPRESS=DEFLATE",
                           overwrite=overwrite_mode, quiet=quiet_mode)
                r.to_vect(input='stream_i_thinned_%s' % str(ID), output='stream_i_rst_vct_%s' % str(ID), type='line',
                          overwrite=overwrite_mode, quiet=quiet_mode)
                v.out_ogr(input='stream_i_rst_vct_%s' % str(ID), type="line",
                          output='%s/stream_i_rst_vct_%s' % (main_output_dir, str(ID)), format="ESRI_Shapefile",
                          overwrite=overwrite_mode, quiet=quiet_mode)

                input_stream_vct_shp = '%s/stream_i_rst_vct_%s/stream_i_rst_vct_%s.shp' % (
                    tmp_output_dir, str(ID), str(ID))
                v.generalize(input='stream_i_rst_vct_%s' % str(ID), output='stream_i_rst_vct_gen_%s' % str(ID),
                             method='reduction', threshold='1000000000', overwrite=overwrite_mode, quiet=quiet_mode)
                # v.out_ogr(input='stream_i_rst_vct_gen_%s'%str(ID), type="line", output='%s/stream_i_rst_vct_gen_%s'%(tmp_output_dir,str(ID)), format="ESRI_Shapefile", overwrite = overwrite_mode, quiet = quiet_mode)

                # Buffer streams
                r.buffer(input='stream_i_%s' % str(ID), output='stream_buff_i_%s' % str(ID),
                         distances='%s' % (str(buffer_distance_cells)), overwrite=overwrite_mode, quiet=quiet_mode)

                # # First method clip over stream buffer center line boundary
                # r.mapcalc("stream_buff_i_fc_%s=if(stream_buff_i_%s==1,1,null())"%(str(ID),str(ID)), overwrite = overwrite_mode, quiet = quiet_mode)
                # r.mapcalc("stream_buff_sv_i_%s=if(isnull(stream_buff_i_fc_%s),null(),1)"%(str(ID),str(ID)), overwrite = overwrite_mode, quiet = quiet_mode)
                # g.region(raster='stream_buff_sv_i_%s'%str(ID), quiet = quiet_mode)
                # g.region(zoom='stream_buff_sv_i_%s'%str(ID), quiet = quiet_mode)

                # Clip over stream buffer
                r.mapcalc("stream_buff_sv_i_%s=if(isnull(stream_buff_i_%s),null(),1)" % (str(ID), str(ID)),
                          overwrite=overwrite_mode, quiet=quiet_mode)
                g.region(raster='stream_buff_sv_i_%s' % str(ID), quiet=quiet_mode)
                g.region(zoom='stream_buff_sv_i_%s' % str(ID), quiet=quiet_mode)
                r.mask(raster='stream_buff_sv_i_%s' % str(ID))

                # Clip DEM
                r.mapcalc("dem_i_%s=raw_dem" % str(ID), overwrite=overwrite_mode, quiet=quiet_mode)
                # r.out_gdal(input='dem_i_%s'%str(ID), output='%s/dem_i_%s.tif'%(tmp_output_dir,ID), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)

                r.mask(flags='r')

                # Extract points to process
                v.to_points_new(input='stream_i_rst_vct_gen_%s' % str(ID),
                                output='stream_i_vect_start_point_%s' % str(ID),
                                use='start', overwrite=overwrite_mode, quiet=quiet_mode)
                v.to_points_new(input='stream_i_rst_vct_%s' % str(ID), output='stream_i_vect_all_point_%s' % str(ID),
                                use='vertex', overwrite=overwrite_mode, quiet=quiet_mode)
                v.to_points_new(input='stream_i_rst_vct_gen_%s' % str(ID),
                                output='stream_i_vect_end_point_%s' % str(ID),
                                use='end', overwrite=overwrite_mode, quiet=quiet_mode)

                # Make mask start + end point
                v.to_rast(input='stream_i_vect_start_point_%s' % str(ID),
                          output='stream_i_vect_start_point_rst_%s' % str(ID), type='point', use='attr',
                          attribute_column='value', overwrite=overwrite_mode, quiet=quiet_mode)
                v.to_rast(input='stream_i_vect_end_point_%s' % str(ID),
                          output='stream_i_vect_end_point_rst_%s' % str(ID),
                          type='point', use='attr', attribute_column='value', overwrite=overwrite_mode,
                          quiet=quiet_mode)
                r.mapcalc("start_point_mask_%s=if(!isnull(stream_i_vect_start_point_rst_%s),1,0)" % (str(ID), str(ID)),
                          overwrite=overwrite_mode, quiet=quiet_mode)
                r.mapcalc("end_point_mask_%s=if(!isnull(stream_i_vect_end_point_rst_%s),1,0)" % (str(ID), str(ID)),
                          overwrite=overwrite_mode, quiet=quiet_mode)
                r.mapcalc("startend_mask_%s=start_point_mask_%s+end_point_mask_%s" % (str(ID), str(ID), str(ID)),
                          overwrite=overwrite_mode, quiet=quiet_mode)

                # v.out_ogr(input='stream_i_vect_start_point_%s'%str(ID), type="point", output='%s/stream_i_vect_start_point_%s'%(tmp_output_dir,str(ID)), format="ESRI_Shapefile", overwrite = overwrite_mode, quiet = quiet_mode)
                # v.out_ogr(input='stream_i_vect_all_point_%s'%str(ID), type="point", output='%s/stream_i_vect_all_point_%s'%(tmp_output_dir,str(ID)), format="ESRI_Shapefile", overwrite = overwrite_mode, quiet = quiet_mode)
                # v.out_ogr(input='stream_i_vect_end_point_%s'%str(ID), type="point", output='%s/stream_i_vect_end_point_%s'%(tmp_output_dir,str(ID)), format="ESRI_Shapefile", overwrite = overwrite_mode, quiet = quiet_mode)

                r.mapcalc("start_point_elev_%s=if(!isnull(stream_i_vect_start_point_rst_%s),raw_dem,null())" % (
                    str(ID), str(ID)), overwrite=overwrite_mode, quiet=quiet_mode)
                r.mapcalc(
                    "end_point_elev_%s=if(!isnull(stream_i_vect_end_point_rst_%s),raw_dem,null())" % (
                        str(ID), str(ID)),
                    overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='start_point_elev_%s' % str(ID),
                           output='%s/start_point_elev_%s.tif' % (tmp_output_dir, ID), createopt="COMPRESS=DEFLATE",
                           overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='end_point_elev_%s' % str(ID),
                           output='%s/end_point_elev_%s.tif' % (tmp_output_dir, ID), createopt="COMPRESS=DEFLATE",
                           overwrite=overwrite_mode, quiet=quiet_mode)

                r.to_vect(input='start_point_elev_%s' % str(ID), output='start_point_elev_point_%s' % str(ID),
                          type='point', overwrite=overwrite_mode, quiet=quiet_mode)
                r.to_vect(input='end_point_elev_%s' % str(ID), output='end_point_elev_point_%s' % str(ID),
                          type='point',
                          overwrite=overwrite_mode, quiet=quiet_mode)
                v.out_ogr(input='start_point_elev_point_%s' % str(ID), type="point",
                          output='%s/start_point_elev_point_%s' % (tmp_output_dir, str(ID)),
                          format="ESRI_Shapefile",
                          overwrite=overwrite_mode, quiet=quiet_mode)
                v.out_ogr(input='end_point_elev_point_%s' % str(ID), type="point",
                          output='%s/end_point_elev_point_%s' % (tmp_output_dir, str(ID)), format="ESRI_Shapefile",
                          overwrite=overwrite_mode, quiet=quiet_mode)

                v.db_dropcolumn(map='stream_i_vect_start_point_%s' % str(ID), columns='value', quiet=quiet_mode)
                v.db_dropcolumn(map='stream_i_vect_start_point_%s' % str(ID), columns='label', quiet=quiet_mode)
                v.db_addcolumn(map='stream_i_vect_start_point_%s' % str(ID), columns="POS integer")
                v.db_addcolumn(map='stream_i_vect_start_point_%s' % str(ID), columns="value double precision")
                v.what_rast(map='stream_i_vect_start_point_%s' % str(ID), raster='raw_dem', column='value',
                            quiet=quiet_mode)
                v.out_ogr(input='stream_i_vect_start_point_%s' % str(ID), type="point",
                          output='%s/stream_i_vect_start_point_v2_%s' % (tmp_output_dir, str(ID)),
                          format="ESRI_Shapefile", overwrite=overwrite_mode, quiet=quiet_mode)
                stream_i_vect_start_point_shp = '%s/stream_i_vect_start_point_v2_%s/%s.shp' % (
                    tmp_output_dir, str(ID), 'stream_i_vect_start_point_%s' % str(ID))
                stream_i_vect_start_point_shape = gpd.read_file(stream_i_vect_start_point_shp)
                stream_i_vect_start_point_shape_ed = stream_i_vect_start_point_shape.copy()
                stream_i_vect_start_point_shape_ed = stream_i_vect_start_point_shape_ed.assign(POS=1)
                stream_i_vect_start_point_shape_ed = stream_i_vect_start_point_shape_ed.assign(cat=1)
                value_tmp = stream_i_vect_start_point_shape_ed.loc[
                    stream_i_vect_start_point_shape_ed['POS'] == 1, 'value']
                start_point_elevation = value_tmp.values[0]
                stream_i_vect_start_point_shape_ed = stream_i_vect_start_point_shape_ed[
                    ['cat', 'value', 'POS', 'geometry']]
                stream_i_vect_start_point_shp_v2 = '%s/%s_v2.shp' % (
                    tmp_output_dir, 'stream_i_vect_start_point_%s' % str(ID))
                stream_i_vect_start_point_shape_ed.to_file(driver='ESRI Shapefile',
                                                           filename=stream_i_vect_start_point_shp_v2)
                del stream_i_vect_start_point_shape, stream_i_vect_start_point_shape_ed, value_tmp

                v.db_dropcolumn(map='stream_i_vect_end_point_%s' % str(ID), columns='value', quiet=quiet_mode)
                v.db_dropcolumn(map='stream_i_vect_end_point_%s' % str(ID), columns='label', quiet=quiet_mode)
                v.db_addcolumn(map='stream_i_vect_end_point_%s' % str(ID), columns="POS integer")
                v.db_addcolumn(map='stream_i_vect_end_point_%s' % str(ID), columns="value double precision")
                v.what_rast(map='stream_i_vect_end_point_%s' % str(ID), raster='raw_dem', column='value',
                            quiet=quiet_mode)
                v.out_ogr(input='stream_i_vect_end_point_%s' % str(ID), type="point",
                          output='%s/stream_i_vect_end_point_v2_%s' % (tmp_output_dir, str(ID)),
                          format="ESRI_Shapefile", overwrite=overwrite_mode, quiet=quiet_mode)
                stream_i_vect_end_point_shp = '%s/stream_i_vect_end_point_v2_%s/%s.shp' % (
                    tmp_output_dir, str(ID), 'stream_i_vect_end_point_%s' % str(ID))
                stream_i_vect_end_point_shape = gpd.read_file(stream_i_vect_end_point_shp)
                stream_i_vect_end_point_shape_ed = stream_i_vect_end_point_shape.copy()
                stream_i_vect_end_point_shape_ed = stream_i_vect_end_point_shape_ed.assign(POS=1)
                stream_i_vect_end_point_shape_ed = stream_i_vect_end_point_shape_ed.assign(cat=1)
                value_tmp = stream_i_vect_end_point_shape_ed.loc[
                    stream_i_vect_end_point_shape_ed['POS'] == 1, 'value']
                end_point_elevation = value_tmp.values[0]
                stream_i_vect_end_point_shape_ed = stream_i_vect_end_point_shape_ed[
                    ['cat', 'value', 'POS', 'geometry']]
                stream_i_vect_end_point_shp_v2 = '%s/%s_v2.shp' % (
                    tmp_output_dir, 'stream_i_vect_end_point_%s' % str(ID))
                stream_i_vect_end_point_shape_ed.to_file(driver='ESRI Shapefile',
                                                         filename=stream_i_vect_end_point_shp_v2)
                del stream_i_vect_end_point_shape, stream_i_vect_end_point_shape_ed, value_tmp

                if float(start_point_elevation) >= float(end_point_elevation):

                    v.in_ogr(flags="o", input=stream_i_vect_start_point_shp_v2,
                             output='real_inlet_point_%s' % str(ID),
                             quiet=quiet_mode)
                    v.in_ogr(flags="o", input=stream_i_vect_end_point_shp_v2,
                             output='real_outlet_point_%s' % str(ID),
                             quiet=quiet_mode)
                    v.out_ogr(input='real_inlet_point_%s' % str(ID), type="point",
                              output='%s/real_inlet_point_%s' % (tmp_output_dir, str(ID)), format="ESRI_Shapefile",
                              overwrite=overwrite_mode, quiet=quiet_mode)
                    v.out_ogr(input='real_outlet_point_%s' % str(ID), type="point",
                              output='%s/real_outlet_point_%s' % (tmp_output_dir, str(ID)), format="ESRI_Shapefile",
                              overwrite=overwrite_mode, quiet=quiet_mode)
                else:
                    v.in_ogr(flags="o", input=stream_i_vect_start_point_shp_v2,
                             output='real_outlet_point_%s' % str(ID),
                             quiet=quiet_mode)
                    v.in_ogr(flags="o", input=stream_i_vect_end_point_shp_v2,
                             output='real_inlet_point_%s' % str(ID),
                             quiet=quiet_mode)
                    v.out_ogr(input='real_inlet_point_%s' % str(ID), type="point",
                              output='%s/real_inlet_point_%s' % (tmp_output_dir, str(ID)), format="ESRI_Shapefile",
                              overwrite=overwrite_mode, quiet=quiet_mode)
                    v.out_ogr(input='real_outlet_point_%s' % str(ID), type="point",
                              output='%s/real_outlet_point_%s' % (tmp_output_dir, str(ID)), format="ESRI_Shapefile",
                              overwrite=overwrite_mode, quiet=quiet_mode)

                # Create masks
                v.to_rast(input='real_inlet_point_%s' % str(ID), output='inlet_mask_%s' % str(ID), type='point',
                          use='attr', attribute_column='POS', overwrite=overwrite_mode, quiet=quiet_mode)
                v.to_rast(input='real_outlet_point_%s' % str(ID), output='outlet_mask_%s' % str(ID), type='point',
                          use='attr', attribute_column='POS', overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='inlet_mask_%s' % str(ID), output='%s/inlet_mask_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='outlet_mask_%s' % str(ID), output='%s/outlet_mask_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)
                r.mapcalc("stream_mask_%s=if(!isnull(stream_i_thinned_%s),0,null())" % (str(ID), str(ID)),
                          overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='stream_mask_%s' % str(ID), output='%s/stream_mask_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)
                r.patch(input='inlet_mask_%s,outlet_mask_%s,stream_mask_%s' % (str(ID), str(ID), str(ID)),
                        output='inoutlet_mask_%s' % str(ID), overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='inoutlet_mask_%s' % str(ID),
                           output='%s/inoutlet_mask_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)

                r.mapcalc("stream_i_thinned_elev_rst_%s=if(stream_i_thinned_%s==%s,raw_dem,null())" % (
                    str(ID), str(ID), str(ID)), overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='stream_i_thinned_elev_rst_%s' % str(ID),
                           output='%s/stream_i_thinned_elev_rst_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)

                # Extract profile along stream (elevation VS progressive distance from inlet)
                r.cost(input='stream_i_thinned_%s' % str(ID), output='progressive_%s' % str(ID),
                       start_raster='inlet_mask_%s' % str(ID), overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='progressive_%s' % str(ID), output='%s/progressive_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)

                # Calculate statistics
                progressive_i_stats = grass.parse_command('r.univar', map='progressive_%s' % str(ID), flags='g')
                progressive_i_count = int(progressive_i_stats['n'])  # To verify how to get non-null cells
                progressive_i_min = float(progressive_i_stats['min'])
                progressive_i_max = float(progressive_i_stats['max'])

                # Calculate standardized flow_acc
                r.mapcalc("progressive_adim_%s=(progressive_%s-%s)/(%s-%s)" % (
                    str(ID), str(ID), str(progressive_i_min), str(progressive_i_max), str(progressive_i_min)),
                          overwrite=overwrite_mode, quiet=quiet_mode)

                # Multiply flow_acc_i_adim by large number to avoid having problems reading small values E^x
                r.mapcalc("progressive_adim_%s=progressive_adim_%s*1000" % (str(ID), str(ID)),
                          overwrite=overwrite_mode,
                          quiet=quiet_mode)
                r.out_gdal(input='progressive_adim_%s' % str(ID),
                           output='%s/progressive_adim_%s.tif' % (tmp_output_dir, ID), createopt="COMPRESS=DEFLATE",
                           overwrite=overwrite_mode, quiet=quiet_mode)

                r.mapcalc("progressive_adim_no_inoutlet_%s=if(inoutlet_mask_%s==0,progressive_adim_%s,null())" % (
                    str(ID), str(ID), str(ID)), overwrite=overwrite_mode, quiet=quiet_mode)

                # Mask out zero values
                r.mapcalc("inlet_mask_%s=if(inlet_mask_%s==1,1,null())" % (str(ID), str(ID)),
                          overwrite=overwrite_mode,
                          quiet=quiet_mode)
                r.mapcalc("outlet_mask_%s=if(outlet_mask_%s==1,1,null())" % (str(ID), str(ID)),
                          overwrite=overwrite_mode, quiet=quiet_mode)

                # Create raster inlet and outlet position
                r.mapcalc("inlet_pos_%s=inlet_mask_%s*0" % (str(ID), str(ID)), overwrite=overwrite_mode,
                          quiet=quiet_mode)
                r.mapcalc("outlet_pos_%s=outlet_mask_%s*(%s-1)" % (str(ID), str(ID), progressive_i_count),
                          overwrite=overwrite_mode, quiet=quiet_mode)

                # Get outlet and inlet elevation
                r.mapcalc("inlet_elev_%s=inlet_mask_%s*raw_dem" % (str(ID), str(ID)), overwrite=overwrite_mode,
                          quiet=quiet_mode)
                r.out_gdal(input='inlet_elev_%s' % str(ID), output='%s/inlet_elev_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)
                r.mapcalc("outlet_elev_%s=outlet_mask_%s*raw_dem" % (str(ID), str(ID)), overwrite=overwrite_mode,
                          quiet=quiet_mode)
                r.out_gdal(input='outlet_elev_%s' % str(ID), output='%s/outlet_elev_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)

                # Convert raster to point vector
                r.to_vect(input='inlet_elev_%s' % str(ID), output='inlet_point_%s' % str(ID), type='point',
                          overwrite=overwrite_mode, quiet=quiet_mode)
                v.out_ogr(input='inlet_point_%s' % str(ID), type="point",
                          output='%s/inlet_elev_vect_%s' % (tmp_output_dir, str(ID)), format="ESRI_Shapefile",
                          overwrite=overwrite_mode, quiet=quiet_mode)
                inlet_point_shp = '%s/inlet_elev_vect_%s/inlet_point_%s.shp' % (tmp_output_dir, str(ID), str(ID))

                # Replace standard ID with the position one
                v.in_ogr(input=inlet_point_shp, output='inlet_point_%s_v2' % str(ID), quiet=quiet_mode)
                v.what_rast(map='inlet_point_%s_v2' % str(ID), raster='inlet_pos_%s' % str(ID), column='cat_',
                            quiet=quiet_mode)
                v.db_renamecolum(map='inlet_point_%s_v2' % str(ID), column=('cat_', 'POS'), quiet=quiet_mode)
                v.out_ogr(input='inlet_point_%s_v2' % str(ID), type="point",
                          output='%s/inlet_elev_vect_%s_v2' % (tmp_output_dir, str(ID)), format="ESRI_Shapefile",
                          overwrite=overwrite_mode, quiet=quiet_mode)
                inlet_point_shp_v2 = '%s/inlet_elev_vect_%s_v2/inlet_point_%s_v2.shp' % (
                    tmp_output_dir, str(ID), str(ID))

                # Convert raster to point vector
                r.to_vect(input='outlet_elev_%s' % str(ID), output='outlet_point_%s' % str(ID), type='point',
                          overwrite=overwrite_mode, quiet=quiet_mode)
                v.out_ogr(input='outlet_point_%s' % str(ID), type="point",
                          output='%s/outlet_elev_vect_%s' % (tmp_output_dir, str(ID)), format="ESRI_Shapefile",
                          overwrite=overwrite_mode, quiet=quiet_mode)
                outlet_point_shp = '%s/outlet_elev_vect_%s/outlet_point_%s.shp' % (tmp_output_dir, str(ID), str(ID))

                # Replace standard ID with the position one
                v.in_ogr(flags="o", input=outlet_point_shp, output='outlet_point_%s_v2' % str(ID), quiet=quiet_mode)
                v.what_rast(map='outlet_point_%s_v2' % str(ID), raster='outlet_pos_%s' % str(ID), column='cat_',
                            quiet=quiet_mode)
                v.db_renamecolum(map='outlet_point_%s_v2' % str(ID), column=('cat_', 'POS'), quiet=quiet_mode)
                v.out_ogr(input='outlet_point_%s_v2' % str(ID), type="point",
                          output='%s/outlet_elev_vect_%s_v2' % (tmp_output_dir, str(ID)), format="ESRI_Shapefile",
                          overwrite=overwrite_mode, quiet=quiet_mode)
                outlet_point_shp_v2 = '%s/outlet_elev_vect_%s_v2/outlet_point_%s_v2.shp' % (
                    tmp_output_dir, str(ID), str(ID))

                # Convert raster to point vector
                r.to_vect(input='progressive_adim_no_inoutlet_%s' % str(ID), output='extra_points_%s' % str(ID),
                          type='point', overwrite=overwrite_mode, quiet=quiet_mode)
                v.out_ogr(input='extra_points_%s' % str(ID), type="point",
                          output='%s/extra_vect_%s' % (tmp_output_dir, str(ID)), format="ESRI_Shapefile",
                          overwrite=overwrite_mode, quiet=quiet_mode)
                extra_points_shp = '%s/extra_vect_%s/extra_points_%s.shp' % (tmp_output_dir, str(ID), str(ID))

                # Sort values ascending and export sorted shapefile with new ID. ID is the cell position in downstream direction from 1 to non null cells (withouth inlet and outlet)
                extra_points_shp_ordered = '%s/extra_vect_orderd_%s.shp' % (tmp_output_dir, str(ID))
                shape = gpd.read_file(extra_points_shp)
                shape = shape.assign(POS=0)
                shape_sorted = shape.iloc[shape['value'].sort_values(ascending=[True]).index.values]
                all_values = shape_sorted['value'].values
                count = 0
                for value_i in all_values:
                    count += 1
                    shape_sorted.loc[shape_sorted['value'] == value_i, 'POS'] = count
                shape_sorted = shape_sorted[['cat', 'value', 'POS', 'geometry']]
                shape_sorted.to_file(driver='ESRI Shapefile', filename=extra_points_shp_ordered)
                v.in_ogr(flags="o", input=extra_points_shp_ordered, output='stream_points_ordered_%s' % str(ID),
                         quiet=quiet_mode)

                # Delete cat_ field
                v.db_dropcolumn(map='stream_points_ordered_%s' % str(ID), columns='cat_', quiet=quiet_mode)

                # Extract elevation from raster at point location (replace values into "value" Field)
                v.what_rast(map='stream_points_ordered_%s' % str(ID), raster='dem_i_%s' % str(ID), column='value',
                            quiet=quiet_mode)
                v.out_ogr(input='stream_points_ordered_%s' % str(ID), type="point",
                          output='%s/stream_elev_vect_%s' % (tmp_output_dir, str(ID)), format="ESRI_Shapefile",
                          overwrite=overwrite_mode, quiet=quiet_mode)
                stream_elev_points_shp = '%s/stream_elev_vect_%s/stream_points_ordered_%s.shp' % (
                    tmp_output_dir, str(ID), str(ID))

                # Extract profile array
                shape_str_elev = gpd.read_file(stream_elev_points_shp)
                shape_str_elev_sorted = shape_str_elev.iloc[
                    shape_str_elev['POS'].sort_values(ascending=[True]).index.values]
                str_elev_x = shape_str_elev_sorted['POS'].values
                str_elev_z = shape_str_elev_sorted['value'].values
                str_elev_geom = shape_str_elev_sorted['geometry'].values
                del shape_str_elev

                str_elev_x_lst = str_elev_x.tolist()
                str_elev_z_lst = str_elev_z.tolist()

                shape_inlet_elev = gpd.read_file(inlet_point_shp_v2)
                inlet_elev_x = shape_inlet_elev['POS'].values
                inlet_elev_z = shape_inlet_elev['value'].values
                inlet_elev_geom = shape_inlet_elev['geometry'].values
                del shape_inlet_elev
                inlet_elev_x_lst = inlet_elev_x.tolist()
                inlet_elev_z_lst = inlet_elev_z.tolist()

                shape_outlet_elev = gpd.read_file(outlet_point_shp_v2)
                outlet_elev_x = shape_outlet_elev['POS'].values
                outlet_elev_z = shape_outlet_elev['value'].values
                outlet_elev_geom = shape_outlet_elev['geometry'].values
                del shape_outlet_elev
                outlet_elev_x_lst = outlet_elev_x.tolist()
                outlet_elev_z_lst = outlet_elev_z.tolist()

                all_stream_elev_x = np.concatenate((inlet_elev_x, str_elev_x, outlet_elev_x))
                all_stream_elev_z = np.concatenate((inlet_elev_z, str_elev_z, outlet_elev_z))
                all_stream_elev_geom = np.concatenate((inlet_elev_geom, str_elev_geom, outlet_elev_geom))
                all_stream_elev_x_lst = all_stream_elev_x.tolist()
                all_stream_elev_z_lst = all_stream_elev_z.tolist()

                str_elev_x_m = all_stream_elev_x + 1
                str_elev_x_m = str_elev_x_m[:-2]
                str_elev_z_m = all_stream_elev_z[:-2]
                str_elev_x_m_lst = str_elev_x_m.tolist()
                str_elev_z_m_lst = str_elev_z_m.tolist()

                str_elev_x_p = all_stream_elev_x - 1
                str_elev_x_p = str_elev_x_p[2:]
                str_elev_z_p = all_stream_elev_z[2:]
                str_elev_x_p_lst = str_elev_x_p.tolist()
                str_elev_z_p_lst = str_elev_z_p.tolist()

                # Compute z-z_m and z_p-z
                delta_z_m = str_elev_z - str_elev_z_m
                delta_z_p = str_elev_z_p - str_elev_z
                delta_z_m_lst = delta_z_m.tolist()
                delta_z_p_lst = delta_z_p.tolist()

                # Compute |(dz_p-dz_m)*n_pixels|
                Dz_2 = abs((delta_z_p - delta_z_m) * 2)
                Dz_2_lst = Dz_2.tolist()

                # Compute average DZ
                avg_Dz = (inlet_elev_z[0] - outlet_elev_z[0]) / (progressive_i_count - 1)

                A = np.vstack((str_elev_x, str_elev_z, Dz_2)).T
                avg_Dz = avg_Dz / filtering_magnitude
                F = A[A[:, 2] < avg_Dz]

                X_filt = F[:, 0]
                Z_filt = F[:, 1]

                X_filt_lst = X_filt.tolist()
                Z_filt_lst = Z_filt.tolist()

                all_stream_filt_elev_x = np.concatenate((inlet_elev_x, X_filt, outlet_elev_x))
                all_stream_filt_elev_z = np.concatenate((inlet_elev_z, Z_filt, outlet_elev_z))
                all_stream_filt_elev_x_lst = all_stream_filt_elev_x.tolist()
                all_stream_filt_elev_z_lst = all_stream_filt_elev_z.tolist()

                window_length = progressive_i_count / 100.0
                odd_window_length_int = int(np.ceil(window_length / 2.) * 2 + 1)
                odd_window_length_int = 3 * odd_window_length_int
                # logging.info('window_length')
                # logging.info(window_length)
                # logging.info('progressive_i_count')
                # logging.info(progressive_i_count)
                # logging.info('odd_window_length_int')
                # logging.info(odd_window_length_int)

                f = interpolate.interp1d(all_stream_filt_elev_x_lst, all_stream_filt_elev_z_lst)
                str_elev_x_interp_lst = np.arange(0, progressive_i_count, 1)
                str_elev_z_interp_lst = f(str_elev_x_interp_lst)
                str_elev_x_interp = np.array(str_elev_x_interp_lst)
                str_elev_z_interp = np.array(str_elev_z_interp_lst)

                str_elev_z_interp_sm_lst = savgol_filter(str_elev_z_interp_lst, odd_window_length_int, 3)
                str_elev_z_interp_sm = np.array(str_elev_z_interp_sm_lst)

                str_elev_z_interp_adim = (
                        (all_stream_filt_elev_z - outlet_elev_z[0]) / (inlet_elev_z[0] - outlet_elev_z[0]))

                if flat_burning_mode_flag == False:
                    str_elev_z_burned = (all_stream_filt_elev_z - ((1 - str_elev_z_interp_adim) * burning_depth))
                    str_elev_z_burned_lst = str_elev_z_burned.tolist()
                else:
                    str_elev_z_burned = (all_stream_filt_elev_z - burning_depth)
                    str_elev_z_burned_lst = str_elev_z_burned.tolist()

                f = interpolate.interp1d(all_stream_filt_elev_x_lst, str_elev_z_burned_lst)
                str_elev_x_burn_interp_lst = np.arange(0, progressive_i_count, 1)
                str_elev_z_burn_interp_lst = f(str_elev_x_burn_interp_lst)
                str_elev_x_burn_interp = np.array(str_elev_x_burn_interp_lst)
                str_elev_z_burn_interp = np.array(str_elev_z_burn_interp_lst)

                str_elev_z_burn_interp_sm_lst = savgol_filter(str_elev_z_burn_interp, odd_window_length_int, 3)
                str_elev_z_burn_interp_sm = np.array(str_elev_z_burn_interp_sm_lst)

                out_csv_file = '%s/river_profile_%s_orig.csv' % (tmp_output_dir, str(ID))
                with open(out_csv_file, "w") as outfile:
                    entry = 'POS,Z'
                    outfile.write(entry)
                    outfile.write("\n")
                with open(out_csv_file, "a") as outfile:
                    for i in range(len(all_stream_elev_x_lst)):
                        entry = '%s,%s' % (str(all_stream_elev_x_lst[i]), str(all_stream_elev_z_lst[i]))
                        outfile.write(entry)
                        outfile.write("\n")

                out_csv_file_fil = '%s/river_profile_%s_filtered_%s_%s.csv' % (
                    tmp_output_dir, str(ID), str(buffer_distance_cells), str(filtering_magnitude))
                with open(out_csv_file_fil, "w") as outfile:
                    entry = 'POS_INT,Z_INT'
                    outfile.write(entry)
                    outfile.write("\n")
                with open(out_csv_file_fil, "a") as outfile:
                    for i in range(len(str_elev_x_interp_lst)):
                        entry = '%s,%s' % (str(str_elev_x_interp_lst[i]), str(str_elev_z_interp_lst[i]))
                        outfile.write(entry)
                        outfile.write("\n")

                out_csv_file_fil = '%s/river_profile_%s_filtered_sm_%s_%s.csv' % (
                    tmp_output_dir, str(ID), str(buffer_distance_cells), str(filtering_magnitude))
                with open(out_csv_file_fil, "w") as outfile:
                    entry = 'POS_INT_SM,Z_INT_SM'
                    outfile.write(entry)
                    outfile.write("\n")
                with open(out_csv_file_fil, "a") as outfile:
                    for i in range(len(str_elev_x_interp_lst)):
                        entry = '%s,%s' % (str(str_elev_x_interp_lst[i]), str(str_elev_z_interp_sm_lst[i]))
                        outfile.write(entry)
                        outfile.write("\n")

                out_csv_file_bur = '%s/river_profile_%s_burned_%s_%s_%s.csv' % (
                    main_output_dir, str(ID), str(buffer_distance_cells), str(filtering_magnitude), str(burning_depth))
                with open(out_csv_file_bur, "w") as outfile:
                    entry = 'POS_BUR,Z_BUR'
                    outfile.write(entry)
                    outfile.write("\n")
                with open(out_csv_file_bur, "a") as outfile:
                    for i in range(len(str_elev_x_burn_interp_lst)):
                        entry = '%s,%s' % (str(str_elev_x_burn_interp_lst[i]), str(str_elev_z_burn_interp_lst[i]))
                        outfile.write(entry)
                        outfile.write("\n")

                out_csv_file_bur = '%s/river_profile_%s_burned_sm_%s_%s_%s.csv' % (
                    tmp_output_dir, str(ID), str(buffer_distance_cells), str(filtering_magnitude), str(burning_depth))
                with open(out_csv_file_bur, "w") as outfile:
                    entry = 'POS_BUR_SM,Z_BUR_SM'
                    outfile.write(entry)
                    outfile.write("\n")
                with open(out_csv_file_bur, "a") as outfile:
                    for i in range(len(str_elev_x_burn_interp_lst)):
                        entry = '%s,%s' % (
                            str(str_elev_x_burn_interp_lst[i]), str(str_elev_z_burn_interp_sm_lst[i]))
                        outfile.write(entry)
                        outfile.write("\n")

                # Make a river 2D graph
                plt.plot(all_stream_elev_x_lst, all_stream_elev_z_lst, linewidth=0.3, label='Raw DEM',
                         color='black')
                plt.plot(str_elev_x_interp_lst, str_elev_z_interp_lst, linewidth=0.5, label='Filtered',
                         color='blue')
                plt.plot(str_elev_x_burn_interp_lst, str_elev_z_interp_sm_lst, linewidth=0.3, linestyle='dashed',
                         label='Filtered sm', color='blue')
                plt.plot(str_elev_x_burn_interp_lst, str_elev_z_burn_interp_lst, linewidth=0.5, label='Burned',
                         color='brown')
                plt.plot(str_elev_x_burn_interp_lst, str_elev_z_burn_interp_sm_lst, linewidth=0.3,
                         linestyle='dashed',
                         label='Burned sm', color='brown')
                plt.xlabel("Stream vertex no.")
                plt.ylabel("Elevation (m)")
                plt.title('Stream no. %s' % (str(ID)))
                plt.legend()
                plt.savefig('%s/river_profile_%s_%s_%s.png' % (
                    main_output_dir, str(ID), str(buffer_distance_cells), str(filtering_magnitude)))
                plt.close()

                # Import array and recreate shp
                final_points_shp = '%s/final_vect_%s.shp' % (tmp_output_dir, str(ID))
                str_elev_z_burn_interpm = str_elev_z_burn_interp[:-1]
                str_elev_z_burn_interpmp = str_elev_z_burn_interpm[1:]
                shape_str_elev = gpd.read_file(stream_elev_points_shp)
                shape_str_elev_sorted = shape_str_elev.iloc[
                    shape_str_elev['POS'].sort_values(ascending=[True]).index.values]
                all_pos = shape_str_elev_sorted['POS'].values
                count = 0
                for pos_i in all_pos:
                    shape_str_elev_sorted.loc[shape_str_elev_sorted['POS'] == pos_i, 'value'] = \
                        str_elev_z_burn_interpmp[count]
                    count += 1
                shape_str_elev_sorted = shape_str_elev_sorted[['cat', 'value', 'POS', 'geometry']]
                shape_str_elev_sorted.to_file(driver='ESRI Shapefile', filename=final_points_shp)
                v.in_ogr(flags="o", input=final_points_shp, output='final_points_vector_%s' % str(ID),
                         quiet=quiet_mode)

                r.mapcalc("outlet_elev_burned_%s=outlet_elev_%s-%s" % (str(ID), str(ID), str(burning_depth)),
                          overwrite=overwrite_mode, quiet=quiet_mode)
                v.to_rast(input='final_points_vector_%s' % str(ID), output='final_points_rast_%s' % str(ID),
                          type='point', use='attr', attribute_column='value', overwrite=overwrite_mode,
                          quiet=quiet_mode)
                r.patch(
                    input='inlet_elev_%s,final_points_rast_%s,outlet_elev_burned_%s' % (str(ID), str(ID), str(ID)),
                    output='dem_burned_i_%s' % str(ID), overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='dem_burned_i_%s' % str(ID),
                           output='%s/dem_burned_i_%s.tif' % (tmp_output_dir, ID),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)

                # Calculate depth to use for digging the dem

                r.mapcalc("dem_burned_ext_i_%s=if(!isnull(dem_burned_i_%s),1,null())" % (str(ID), str(ID)),
                          overwrite=overwrite_mode, quiet=quiet_mode)
                # r.out_gdal(input='dem_burned_ext_i_%s'%str(ID), output='%s/dem_burned_ext_i_%s.tif'%(tmp_output_dir,ID), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
                r.mapcalc("raw_dem_clip_burned_ext_i_%s=dem_burned_ext_i_%s*dem_i_%s" % (str(ID), str(ID), str(ID)),
                          overwrite=overwrite_mode, quiet=quiet_mode)
                r.mapcalc(
                    "burning_depth_i_%s=dem_burned_i_%s-raw_dem_clip_burned_ext_i_%s" % (str(ID), str(ID), str(ID)),
                    overwrite=overwrite_mode, quiet=quiet_mode)

                # Final expand buffered of 1 pixel
                # Method 1 linear degression from pixel distance (1 pixel)
                r.to_vect(input='dem_burned_ext_i_%s' % str(ID), output='dem_burning_line_i_%s' % str(ID), type='line',
                          overwrite=overwrite_mode, quiet=quiet_mode)
                v.generalize(input='dem_burning_line_i_%s' % str(ID), output='dem_burning_line_gen_%s' % str(ID),
                             method='reduction', threshold='1000000000', overwrite=overwrite_mode, quiet=quiet_mode)
                v.to_points_new(input='dem_burning_line_gen_%s' % str(ID),
                                output='dem_burning_line_gen_str_pnt_%s' % str(ID), use='start',
                                overwrite=overwrite_mode, quiet=quiet_mode)
                v.to_points_new(input='dem_burning_line_gen_%s' % str(ID),
                                output='dem_burning_line_gen_end_pnt_%s' % str(ID), use='end', overwrite=overwrite_mode,
                                quiet=quiet_mode)
                v.to_rast(input='dem_burning_line_gen_str_pnt_%s' % str(ID),
                          output='dem_burning_line_str_pnt_rst_%s' % str(ID), type='point', use='attr',
                          attribute_column='value', overwrite=overwrite_mode, quiet=quiet_mode)
                r.mapcalc("dem_burning_line_str_pnt_msk_%s=if(!isnull(dem_burning_line_str_pnt_rst_%s),1,0)" % (
                    str(ID), str(ID)), overwrite=overwrite_mode, quiet=quiet_mode)
                # r.out_gdal(input='dem_burning_line_str_pnt_msk_%s'%str(ID), output='%s/dem_burning_line_str_pnt_msk_%s.tif'%(tmp_output_dir,ID), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
                v.to_rast(input='dem_burning_line_gen_end_pnt_%s' % str(ID),
                          output='dem_burning_line_end_pnt_rst_%s' % str(ID), type='point', use='attr',
                          attribute_column='value', overwrite=overwrite_mode, quiet=quiet_mode)
                r.mapcalc("dem_burning_line_end_pnt_msk_%s=if(!isnull(dem_burning_line_end_pnt_rst_%s),1,0)" % (
                    str(ID), str(ID)), overwrite=overwrite_mode, quiet=quiet_mode)
                # r.out_gdal(input='dem_burning_line_end_pnt_msk_%s'%str(ID), output='%s/dem_burning_line_end_pnt_msk_%s.tif'%(tmp_output_dir,ID), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
                r.mapcalc(
                    "dem_half_burned_mask_i_%s=dem_burning_line_str_pnt_msk_%s+dem_burning_line_end_pnt_msk_%s" % (
                        str(ID), str(ID), str(ID)), overwrite=overwrite_mode, quiet=quiet_mode)
                # r.patch(input='dem_burning_line_str_pnt_msk_%s,dem_burning_line_end_pnt_msk_%s'%(str(ID),str(ID)), output='dem_half_burned_mask_i_%s'%str(ID), overwrite = overwrite_mode, quiet = quiet_mode)
                # r.out_gdal(input='dem_half_burned_mask_i_%s'%str(ID), output='%s/dem_half_burned_mask_i_%s.tif'%(tmp_output_dir,ID), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
                r.mapcalc("half_burning_depth_i_%s=burning_depth_i_%s/2.0" % (str(ID), str(ID)),
                          overwrite=overwrite_mode, quiet=quiet_mode)
                # r.out_gdal(input='half_burning_depth_i_%s'%str(ID), output='%s/half_burning_depth_i_%s.tif'%(tmp_output_dir,ID), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
                r.mapcalc(
                    "half_burning_depth_masked_i_%s=if(dem_half_burned_mask_i_%s==1,null(),half_burning_depth_i_%s)" % (
                        str(ID), str(ID), str(ID)), overwrite=overwrite_mode, quiet=quiet_mode)
                # r.out_gdal(input='half_burning_depth_masked_i_%s'%str(ID), output='%s/half_burning_depth_masked_i_%s.tif'%(tmp_output_dir,ID), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
                r.to_vect(input='half_burning_depth_masked_i_%s' % str(ID),
                          output='half_burning_depth_vect_i_%s' % str(ID), type='area', overwrite=overwrite_mode,
                          quiet=quiet_mode)

                if target_epsg != source_epsg:
                    v.buffer(flags="t", input='half_burning_depth_vect_i_%s' % str(ID),
                             output='half_burning_depth_vect_buff_i_%s' % str(ID), distance='%s' % (str(rrs_m)),
                             overwrite=overwrite_mode, quiet=quiet_mode)
                    v.out_ogr(input='half_burning_depth_vect_buff_i_%s' % str(ID), type="area",
                              output='%s/half_burning_depth_vect_buff_i_%s' % (tmp_output_dir, str(ID)),
                              format="ESRI_Shapefile", overwrite=overwrite_mode, quiet=quiet_mode)
                else:
                    v.buffer(flags="t", input='half_burning_depth_vect_i_%s' % str(ID),
                             output='half_burning_depth_vect_buff_i_%s' % str(ID), distance='%s' % (str(in_res_DD)),
                             overwrite=overwrite_mode, quiet=quiet_mode)
                    v.out_ogr(input='half_burning_depth_vect_buff_i_%s' % str(ID), type="area",
                              output='%s/half_burning_depth_vect_buff_i_%s' % (tmp_output_dir, str(ID)),
                              format="ESRI_Shapefile", overwrite=overwrite_mode, quiet=quiet_mode)

                v.to_rast(input='half_burning_depth_vect_buff_i_%s' % str(ID),
                          output='half_burning_depth_vect_buff_rst_i_%s' % str(ID), type='area',
                          overwrite=overwrite_mode, quiet=quiet_mode)
                # r.out_gdal(input='half_burning_depth_vect_buff_rst_i_%s'%str(ID), output='%s/half_burning_depth_vect_buff_rst_i_%s.tif'%(main_output_dir,ID), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
                r.patch(input='burning_depth_i_%s,half_burning_depth_vect_buff_rst_i_%s' % (str(ID), str(ID)),
                        output='burning_depth_extended_i_%s' % str(ID), overwrite=overwrite_mode, quiet=quiet_mode)
                r.out_gdal(input='burning_depth_extended_i_%s' % str(ID),
                           output='%s/burning_depth_extended_i_%s.tif' % (out_segment_dir, str(ID)),
                           createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)

                if os.path.isfile('%s/burning_depth_extended_i_%s.tif' % (
                        out_segment_dir, str(ID))):
                    segm_count += 1
                    if segm_count == 1:
                        segm_rst_list = str('burning_depth_extended_i_%s' % str(ID))
                    else:
                        segm_rst_list = '%s,%s' % (segm_rst_list, 'burning_depth_extended_i_%s' % str(ID))

            # Merge them all
            g.region(raster='raw_dem')
            g.region(zoom='raw_dem')
            r.mapcalc("dem_nan_mask=if(!isnull(raw_dem),1,0)", overwrite=overwrite_mode, quiet=quiet_mode)
            logging.info('-> dem segment conditioned rst list: %s' % segm_rst_list)
            r.mapcalc("burning_values=nmin(%s)" % (segm_rst_list), overwrite=overwrite_mode, quiet=quiet_mode)
            r.mapcalc("out_conditioned_dem_diff=raw_dem+burning_values", overwrite=overwrite_mode, quiet=quiet_mode)
            r.patch(input='out_conditioned_dem_diff,raw_dem', output='out_conditioned_dem_tmp',
                    overwrite=overwrite_mode,
                    quiet=quiet_mode)
            r.mapcalc("out_conditioned_dem=if(dem_nan_mask==1,out_conditioned_dem_tmp,null())",
                      overwrite=overwrite_mode,
                      quiet=quiet_mode)

            # Export and define output variable
            r.out_gdal(input='out_conditioned_dem', output='%s/%s.tif' % (main_output_dir, conditioned_dem),
                       overwrite=overwrite_mode, quiet=quiet_mode)
            # dem_input_cond_file='%s/%s.tif'%(main_output_dir,conditioned_dem)
            out_dem_cond_file = '%s/%s.tif' % (main_output_dir, conditioned_dem)

        else:
            logging.error(
                '-> ERROR! Choose DCOND_METHOD between A: simplified and B: complete!' + method_dem_burning + ' is not a good choice, verify config file!')
            raise NotImplementedError

    else:
        dem_input_cond_file = out_dem_filled_file
        g.region(raster='raw_dem')
        out_dem_cond_file = '%s/%s.tif' % (main_output_dir, conditioned_dem)
        r.external(input=dem_input_cond_file, output='conditioned_dem_rast', overwrite=overwrite_mode, quiet=quiet_mode)
        r.out_gdal(input='conditioned_dem_rast', output=out_dem_cond_file, overwrite=overwrite_mode, quiet=quiet_mode)

    # ------------------------------------------------------------------------------------------------------------------

    ################################################################################

    # Clean temporary files
    g.remove(flags="f", type='raster', pattern='tmp.*', quiet=quiet_mode)

    ################################################################################

    # Estimate total execution time
    tot_elapsed_time_sec = float(time.time() - code_start_time)
    tot_elapsed_time, time_units = Give_Elapsed_Time(tot_elapsed_time_sec)

    # Delete GRASS temporary database including mapset
    if os.path.isdir(grass_temp_db):
        shutil.rmtree(grass_temp_db)
    else:
        pass

    # # Delete output tmp sub-folder
    if os.path.isdir(tmp_output_dir):
        shutil.rmtree(tmp_output_dir)
    else:
        pass

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

# ------------------------------------------------------------------------------------------------------------------
# Project input raster to target EPSG
def Project_geotiff(input_dem_geotiff, source_epsg, target_epsg, output_dir):
    dst_crs = 'EPSG:%s' % (str(target_epsg))
    input_DTM_file_dirname = os.path.dirname(input_dem_geotiff)
    input_DTM_file_basename = os.path.basename(input_dem_geotiff)
    input_DTM_file_basename_we = os.path.splitext(input_DTM_file_basename)[0]
    output_DTM_file_prj = "%s/%s_%s.tif" % (str(output_dir), str(input_DTM_file_basename_we), str(target_epsg))
    input_dem = rasterio.open(input_dem_geotiff)
    source_epsg = input_dem.crs
    logging.info("Input raster %s..." % (str(source_epsg)))
    with rasterio.open(input_dem_geotiff) as src:
        transform, width, height = calculate_default_transform(src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
        with rasterio.open(output_DTM_file_prj, 'w', **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    compress='lzw',
                    nodata='-9999',
                    resampling=Resampling.bilinear)
                # output_dem = rasterio.open(output_DTM_file_prj)

    return output_DTM_file_prj

# ------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------
# Project input vector to target EPSG
def Project_shapefile(input_stream_shapefile, source_epsg, target_epsg, output_dir):
    input_stream_shapefile_dirname = os.path.dirname(input_stream_shapefile)
    input_stream_shapefile_basename = os.path.basename(input_stream_shapefile)
    input_stream_shapefile_basename_we = os.path.splitext(input_stream_shapefile_basename)[0]
    output_prj_stream_shpfile = '%s/%s_%s.shp' % (output_dir, input_stream_shapefile_basename_we, target_epsg)
    output_prj_stream_prjfile = '%s/%s_%s.prj' % (output_dir, input_stream_shapefile_basename_we, target_epsg)
    ogr2ogr_command = "ogr2ogr -s_srs EPSG:%s -t_srs EPSG:%s %s %s" % (
        str(source_epsg), str(target_epsg), output_prj_stream_shpfile, input_stream_shapefile)
    p = subprocess.Popen(ogr2ogr_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if err:
        logging.info('-> GDAL ogr2ogr Error: %s' % err)
    else:
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(target_epsg))
        srs.MorphToESRI()
        file_prj = open(output_prj_stream_prjfile, 'w')
        file_prj.write(srs.ExportToWkt())
        file_prj.close()

    return output_prj_stream_shpfile


# ------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-log_file', action="store", dest="alg_log")
    parser_handle.add_argument('-settings_file', action="store", dest="alg_settings")
    parser_handle.add_argument('-dem', action="store", dest="alg_dem")
    parser_handle.add_argument('-stream', action="store", dest="alg_stream")
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

    if parser_values.alg_dem:
        alg_dem = parser_values.alg_dem
    else:
        alg_dem = None

    if parser_values.alg_stream:
        alg_stream = parser_values.alg_stream
    else:
        alg_stream = None

    if parser_values.alg_output:
        alg_output = parser_values.alg_output
    else:
        alg_output = None

    return alg_settings, alg_log, alg_dem, alg_stream, alg_output

# -------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == "__main__":
    main()
# ------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------
