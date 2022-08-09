#!/usr/bin/env python
"""
REFlEx - Step1 - Hydrogical Derivatives
__date__ = '20220726'
__version__ = '2.0.1'
__author__ =
        'Mauro Arcorace' (mauro.arcorace@cimafoundation.org',
        'Alessandro Masoero (alessandro.masoero@cimafoundation.org',
        'Valerio Basso',
        'Alessia Matanò',
        'Andrea Libertino (andrea.libertino@cimafoundation.org)'
__library__ = 'REFlEx'
General command line:
### python reflex_step1_hydro_derivatives.py -log_file "/path/to/log.txt" -settings_file "settings.json" -base_path "/path/to/base_folder" (-dem "/path/to/dem.tif")
Version(s):
20190220 (1.0.0) --> Beta release
20220406 (2.0.0) --> Full revision - Integrated vector network and pfafstetter assignation for future integration of network editing
                                     Temporary removed MFD support
                                     Adopted pyshed dinf algorithms
20220726 (2.0.2) --> Fix basin delineation procedure
"""
# -------------------------------------------------------------------------------------
# Import python libraries
import os
import numpy as np
import shutil
import sys
import time
import subprocess
from argparse import ArgumentParser
import  logging
from osgeo import osr
from lib.reflex_tools_utils import Give_Elapsed_Time, Start_GRASS_py3, set_logging, read_file_json
from lib.reflex_tools_pfafstetter import assign_pfafstetter_code
from pysheds.grid import Grid
import geopandas as gpd

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'REFlEx - STEP 1 - Hydrological Derivatives'
alg_version = '2.0.1'
alg_release = '2022-07-26'
# Algorithm parameter(s)
time_format = '%Y%m%d%H%M'
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Script Main
def main():
    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    config_file, alg_log, input_dem_geotiff, base_path = get_args()

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
    step0_dir_name = data_settings["step_0"]["dir_name"].replace("{base_path}",base_path)
    step1_dir_name = data_settings["step_1"]["dir_name"].replace("{base_path}",base_path)

    drain_method_streams = data_settings["step_1"]["stream_definition_method"]
    drain_method_hand = data_settings["step_1"]["hand_definition_method"]
    facc_threshold = data_settings["step_1"]["flow_accumulation_threshold_cells"]
    iNum_col = 23               # Number of column of pffstetter in the shapefile
    ################################################################################

    ################################################################################
    
    # Define fullpath of input conditioned and filled DEM from step 0 or use provided dem
    grass_step0_db = step0_dir_name
    
    filled_dem='r_%s_%s_filled_dem'%(domain_name,rrs)
    conditioned_dem='r_%s_%s_conddem'%(domain_name,rrs)

    if input_dem_geotiff is None:
        dem_input_filled_file='%s/%s.tif'%(grass_step0_db,filled_dem)
        dem_input_conditioned_file='%s/%s.tif'%(grass_step0_db,conditioned_dem)
    else:
        dem_input_filled_file = input_dem_geotiff
        dem_input_conditioned_file = input_dem_geotiff
    ################################################################################
    
    
    ################################################################################
        
    # Record starting simulation time
    code_start_time = time.time()

    ################################################################################

    # INITIALIZE GRASS
    # ----------------

    #Start GRASS
    logging.info("--> Starting GRASS..")
    grass_temp_db_name='grass_temp'
    gisbase, grass_temp_db, location, mapset = Start_GRASS_py3(grass_bin, str(target_epsg), grass_temp_db_name, step1_dir_name)
    logging.info('--> GRASS temporary database path: %s'%grass_temp_db)
    logging.info('--> GRASS temporary database directory: %s'%location)
    logging.info('--> GRASS mapset: %s'%mapset)

    # Import GRASS Python bindings
    import grass.script as grass
    import grass.script.setup as gsetup

    # Launch session and do something
    gsetup.init(gisbase, grass_temp_db, location, mapset)
    logging.info('--> GRASS GIS environment: %s'%grass.gisenv())

    # Import GRASS modules
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.pygrass.modules.shortcuts import vector as v
    from grass.pygrass.modules.shortcuts import database as db
    from grass.pygrass.modules import Module

    # Get info about projection
    in_proj = grass.read_command('g.proj', flags = 'jf')
    in_proj = in_proj.strip()
    logging.info("--> GRASS projection parameters: '%s'"%in_proj)

    # Get info about region
    in_region = grass.region()
    logging.info("--> GRASS computational region: '%s'"%in_region)

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

    # ----------------

    ################################################################################

    ################################################################################

    # DEFINE VARIABLE NAMES
    # ---------------------
    
    input_dem='r_%s_%s_dem'%(domain_name,rrs)
    dem_nans='r_%s_%s_dem_nans'%(domain_name,rrs)
    raw_dem='r_%s_%s_raw_dem'%(domain_name,rrs)
    
    conditioned_dem='r_%s_%s_conddem'%(domain_name,rrs)
    filled_dem='r_%s_%s_filled_dem'%(domain_name,rrs)
  
    dem_coverage='r_%s_%s_dem_coverage'%(domain_name,rrs)
    elevation_shade='r_%s_%s_elevshade'%(domain_name,rrs)
    hillshade='r_%s_%s_hillshade'%(domain_name,rrs)
    
    if (drain_method_streams=='d8'):
        flow_dir='r_%s_%s_d8_flow_dir'%(domain_name,rrs)
        flow_dir4hand='r_%s_%s_d8_flow_dir4hand'%(domain_name,rrs)
        flow_acc='r_%s_%s_d8_flowacc'%(domain_name,rrs)
        flow_acc_skm='r_%s_%s_d8_flow_acc_skm'%(domain_name,rrs)
        flow_acc_sm='r_%s_%s_d8_flow_acc_sm'%(domain_name,rrs)
        streamrast='r_%s_%s_d8_streams'%(domain_name,rrs)
        streamrast_bin='r_%s_%s_d8_streamrast_bin'%(domain_name,rrs)
        whatersheds='r_%s_%s_d8_whatersheds'%(domain_name,rrs)
        subbasins='r_%s_%s_d8_subbasins'%(domain_name,rrs)
        mainbasins='r_%s_%s_d8_mainbasins'%(domain_name,rrs)
        sinks='r_%s_%s_d8_sinks'%(domain_name,rrs)
        stream_distance_rst='r_%s_%s_d8_stream_distance'%(domain_name,rrs)
    #elif (drain_method_streams=='MFD'):
    #    flow_dir='r_%s_%s_mfd_flow_dir'%(domain_name,rrs)
    #    flow_dir4hand='r_%s_%s_mfd_flow_dir4hand'%(domain_name,rrs)
    #    flow_acc='r_%s_%s_mfd_flowacc'%(domain_name,rrs)
    #    flow_acc_skm='r_%s_%s_mfd_flow_acc_skm'%(domain_name,rrs)
    #    flow_acc_sm='r_%s_%s_mfd_flow_acc_sm'%(domain_name,rrs)
    #    streamrast='r_%s_%s_mfd_streams'%(domain_name,rrs)
    #    streamrast_bin='r_%s_%s_mfd_streamrast_bin'%(domain_name,rrs)
    #    whatersheds='r_%s_%s_mfd_whatersheds'%(domain_name,rrs)
    #    subbasins='r_%s_%s_mfd_subbasins'%(domain_name,rrs)
    #    mainbasins='r_%s_%s_mfd_mainbasins'%(domain_name,rrs)
    #    sinks='r_%s_%s_mfd_sinks'%(domain_name,rrs)
    #    stream_distance_rst='r_%s_%s_mfd_stream_distance'%(domain_name,rrs)
    else:
        logging.error("Failed in using drainage method %s. Not Supported."%str(drain_method_hand))
        sys.exit(1)

    if (drain_method_hand=='dinf'):
        dinf_flow_dir='r_%s_%s_dinf_flow_dir'%(domain_name,rrs)
        dinf_stream_distance='r_%s_%s_dinf_stream_distance'%(domain_name,rrs)
    else:
        pass
    streamrast2='r_%s_%s_streams4dinf'%(domain_name,rrs)
    cell_area='r_%s_%s_areacell'%(domain_name,rrs)
    strah_streams = 'r_%s_%s_d8_strahler_streams' % (domain_name, rrs)
    stream_vector = 'v_%s_%s_d8_stream_vect' % (domain_name, rrs)
    stream_points = 'v_%s_%s_d8_stream_points' % (domain_name, rrs)
    stream_lines = 'v_%s_%s_d8_stream_lines' % (domain_name, rrs)
    ################################################################################

    ################################################################################
    
    # CREATE OUTPUT DIRECTORIES
    # -------------------------
    
    # Define output MAIN directory
    grass_output_db = step1_dir_name
    logging.info('--> Output main directory: %s'%grass_output_db)
    if not os.path.exists(grass_output_db):
        os.makedirs(grass_output_db)

    # Define output tmp sub-directory
    out_tmp_dirname='tmp'
    out_tmp_dir = os.path.join(grass_output_db, out_tmp_dirname)
    logging.info('--> Output tmp sub-directory: %s'%out_tmp_dir)
    if not os.path.exists(out_tmp_dir):
        os.makedirs(out_tmp_dir)

    # Define output MAIN directory
    logging.info('--> Output main directory: %s'%grass_output_db)
    if not os.path.exists(grass_output_db):
        os.makedirs(grass_output_db)

    # Define output raster sub-directory
    grass_output_db_rst_dirname='rst'
    grass_output_db_rst = os.path.join(grass_output_db, grass_output_db_rst_dirname)
    logging.info('--> Output raster sub-directory: %s'%grass_output_db_rst)
    if not os.path.exists(grass_output_db_rst):
        os.makedirs(grass_output_db_rst)

    # Define output vector sub-directory
    grass_output_db_vct_dirname='vct'
    grass_output_db_vct = os.path.join(grass_output_db, grass_output_db_vct_dirname)
    logging.info('--> Output vector sub-directory: %s'%grass_output_db_vct)
    if not os.path.exists(grass_output_db_vct):
        os.makedirs(grass_output_db_vct)
    ################################################################################

    ################################################################################
    # IMPORT DEM INTO GRASS
    # ---------------------

    # Import DEM file into GRASS
    r.external(input=dem_input_filled_file, output=filled_dem)

    ############################
    # Set computational region
    g.region(raster=filled_dem)
    ############################

    # Compute DEM coverage
    r.mapcalc("%s = isnull(%s)"%(dem_coverage,filled_dem), overwrite = overwrite_mode, quiet = quiet_mode)
    r.mapcalc("%s=not(if(%s))"%(dem_coverage,dem_coverage), overwrite = overwrite_mode, quiet = quiet_mode)
    r.out_gdal(input=dem_coverage, output='%s/%s.tif'%(grass_output_db_rst,dem_coverage), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
    r.out_gdal(input=filled_dem, output='%s/%s.tif' % (grass_output_db_rst, input_dem), type="Float32", createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)

    ################################################################################

    ################################################################################

    # IMPORT CONDITIONED DEM AND CONVERT TO CM
    # ----------------------------------------

    # Import into GRASS conditioned and filled DEMs or re-import raw dem if not found
    r.external(input=dem_input_conditioned_file, output=input_dem)

    dem_slope = 'r_%s_%s_dem_slope' % (domain_name, rrs)
    r.slope_aspect(elevation=filled_dem, slope=dem_slope, format='percent', overwrite=overwrite_mode)
    r.out_gdal(input=dem_slope, output=os.path.join(grass_output_db_rst,dem_slope + ".tif"), createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode,
               quiet=quiet_mode, type='Float32')
    
    # Convert from meter to cm
    input_dem_cm='r_%s_%s_dem_cm'%(domain_name,rrs)
    r.mapcalc("%s=if(%s>=-100,%s*100,null())"%(input_dem_cm,filled_dem,filled_dem), overwrite = overwrite_mode, quiet = quiet_mode)
    out_input_dem_cm_file='%s/%s.tif'%(grass_output_db_rst,input_dem_cm)
    r.out_gdal(input=input_dem_cm, output=out_input_dem_cm_file, createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Extract hydro derivatives

    # Derive from DEM flow direction, flow accumulation and streams information
    logging.info("--> Extracting from DEM flow direction, flow accumulation and streams information..")
    
    # Extract hydrological derivatives
    logging.info("--> Extracting hydrological derivatives from DEM ..")

    if (drain_method_streams=='d8'):
        logging.info("--> Extracting hydrological derivatives from DEM using the SFD (D8) method..")
        r.watershed(flags="sab", elevation=input_dem, threshold='%s'%facc_threshold, accumulation=flow_acc, drainage=flow_dir, basin=whatersheds, stream=streamrast, overwrite = overwrite_mode, quiet = quiet_mode)
        r.watershed(flags="sab", elevation=input_dem_cm, drainage=flow_dir4hand, overwrite = overwrite_mode, quiet = quiet_mode)
        r.out_gdal(input=flow_dir, output='%s/%s.tif'%(grass_output_db_rst,flow_dir), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
        r.out_gdal(input=streamrast, output='%s/%s.tif'%(grass_output_db_rst,streamrast), overwrite = overwrite_mode, quiet = quiet_mode)
        r.out_gdal(input=flow_dir4hand, output='%s/%s.tif'%(grass_output_db_rst,flow_dir4hand), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, nodata=-1, quiet = quiet_mode, flags='f')
        r.out_gdal(input=flow_acc, output='%s/%s.tif' % (grass_output_db_rst, flow_acc), createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)
    #elif (drain_method_streams=='MFD'):
    #    logging.info("--> Extracting hydrological derivatives from DEM using the MFD method....")
    #    r.watershed(elevation=input_dem, threshold='%s'%facc_threshold, accumulation=flow_acc, drainage=flow_dir, basin=whatersheds, stream=streamrast, overwrite = overwrite_mode, quiet = quiet_mode)
    #    r.watershed(elevation=input_dem_cm, drainage=flow_dir4hand, overwrite = overwrite_mode, quiet = quiet_mode)
    #    r.out_gdal(input=flow_dir, output='%s/%s.tif'%(grass_output_db_rst,flow_dir), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
    #    r.out_gdal(input=streamrast, output='%s/%s.tif'%(grass_output_db_rst,streamrast), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
    #    r.out_gdal(input=flow_dir4hand, output='%s/%s.tif'%(grass_output_db_rst,flow_dir4hand), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
    #    r.out_gdal(input=flow_acc, output='%s/%s.tif' % (grass_output_db_rst, flow_acc), createopt="COMPRESS=DEFLATE",overwrite=overwrite_mode, quiet=quiet_mode)
    else:
        logging.error("Failed in using drainage method %s. Not Supported."%str(drain_method_streams))
        sys.exit(1)

    r.mapcalc("slope_channel=if(isnull(" + streamrast + "),null()," + dem_slope + ")", overwrite=overwrite_mode, quiet=quiet_mode)
    r.mapcalc("dem_channel=if(isnull(" + streamrast + "),null()," + input_dem + ")", overwrite=overwrite_mode, quiet=quiet_mode)
    r.out_gdal(input="slope_channel", output=os.path.join(grass_output_db_rst, "r_{domain_name}_{rrs}_slope_channel.tif".format(domain_name=domain_name, rrs=rrs)),
               createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)
    r.out_gdal(input="dem_channel", output=os.path.join(grass_output_db_rst, "r_{domain_name}_{rrs}_dem_channel.tif".format(domain_name=domain_name, rrs=rrs)),
               createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)

    # Convert pixels count to areas using flow accumulation
    logging.info("--> Convert flow accumulation pixels count to areas..")
    r.mapcalc("%s=if(isnull(%s),null(),%s*%s)" % (cell_area, input_dem, rrs_km, rrs_km), overwrite=overwrite_mode,quiet=quiet_mode)
    r.mapcalc("%s=%s*%s" % (flow_acc_skm, flow_acc, cell_area), overwrite=overwrite_mode, quiet=quiet_mode)
    r.mapcalc("%s=%s*%s*10^6" % (flow_acc_sm, flow_acc, cell_area), overwrite=overwrite_mode, quiet=quiet_mode)
    r.out_gdal(input=flow_acc_skm, output='%s/%s.tif' % (grass_output_db_rst, flow_acc_skm),
               createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)
    r.out_gdal(input=flow_acc_sm, output='%s/%s.tif' % (grass_output_db_rst, flow_acc_sm), createopt="COMPRESS=DEFLATE",
               overwrite=overwrite_mode, quiet=quiet_mode)

    # Export whatersheds raster files
    r.out_gdal(input=whatersheds, output='%s/%s.tif'%(grass_output_db_rst,whatersheds), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode, type='Int32', flags='f')
    
    # Delineate basins and subbasins
    logging.info("--> Delineating basins and subbasins..")
    r.stream_basins(direction=flow_dir, stream_rast=streamrast, memory=300, basins=subbasins, overwrite = overwrite_mode, quiet = quiet_mode)
    r.stream_basins(flags="l", direction=flow_dir, stream_rast=streamrast, memory=300, basins=mainbasins, overwrite = overwrite_mode, quiet = quiet_mode)
    r.out_gdal(input=mainbasins, output='%s/%s.tif'%(grass_output_db_rst,mainbasins), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
    r.out_gdal(input=subbasins, output='%s/%s.tif'%(grass_output_db_rst,subbasins), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)


    # IDENTIFY SINKS FROM FLOW DIR
    # ----------------------------
    
    r.mapcalc("%s=if(%s<0,1,null())"%(sinks,flow_dir), overwrite = overwrite_mode)
    r.out_gdal(input=sinks, output='%s/%s.tif'%(grass_output_db_rst,sinks), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)   

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Calculate d-inf derivatives

    if (drain_method_hand=='dinf'):

        # COMPUTE DINF FLOW DIRECTION
        # -------------------------------------

        # Compute D-infinity flow direction:
        logging.info("--> Computing D-infinity Flow Direction and Slope..")

        grid = Grid.from_raster(dem_input_filled_file)
        dem = grid.read_raster(dem_input_filled_file, nodata=-9999)
        fdir = grid.flowdir(dem, routing='dinf')

        dinf_flow_dir_outfile = '%s/%s.tif' % (grass_output_db_rst, dinf_flow_dir)

        grid.to_raster(fdir, dinf_flow_dir_outfile, nodata=-1, dtype=np.float32)
    else:
        pass

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Calculate hand delineation for the whole basin
    # ------------------------------------------
    r.mapcalc("%s=if(%s==1,%s,null())"%(streamrast_bin,dem_coverage,streamrast), overwrite = overwrite_mode, quiet = quiet_mode)
    r.mapcalc("%s=if(isnull(%s),0,1)"%(streamrast_bin,streamrast_bin), overwrite = overwrite_mode, quiet = quiet_mode)
    out_streamrast_bin_file='%s/%s.tif'%(grass_output_db_rst,streamrast_bin)
    r.out_gdal(input=streamrast_bin, output=out_streamrast_bin_file, createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)

    if (drain_method_hand=='dinf'):
        
        # COMPUTE DINF HAND MAP FOR ENTIRE BASIN
        # --------------------------------------
        dinf_hand_map_outfile = '%s/%s.tif' % (grass_output_db_rst, dinf_stream_distance)

        stream = grid.read_raster(out_streamrast_bin_file, nodata=0)
        dem_filled = grid.read_raster(dem_input_filled_file, nodata=-9999)
        hand = grid.compute_hand(fdir=fdir, dem=dem_filled*100, mask=stream, routing='dinf')
        grid.to_raster(hand, dinf_hand_map_outfile, nodata=-1, dtype=np.float32)

        # --------------------------------------

    elif  (drain_method_hand=='d8'): # or (drain_method_hand=='MFD')):

        # COMPUTE D8 OR MFD HAND MAP FOR ENTIRE BASIN
        # -------------------------------------
        logging.info('-> Executing %s distance down average....'%str(drain_method_hand))
        
        r.stream_distance(direction=flow_dir4hand, stream_rast=streamrast_bin, elevation=input_dem_cm, difference=stream_distance_rst, overwrite = overwrite_mode, quiet = quiet_mode)
        r.out_gdal(input=stream_distance_rst, output='%s/%s.tif'%(grass_output_db_rst,stream_distance_rst), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)

        # -------------------------------------

    else:
        logging.error("Failed in using drainage method %s. Not Supported."%str(drain_method_hand))
        sys.exit(1)

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Calculate stream hierarchy

    # Calculates Strahler's streams hierarchy
    logging.info("--> Calculating streams hierarchy..")
    # Replace old flow dir with the one generated with r.accumulated
    r.accumulate(direction=flow_dir, accumulation=flow_acc + "_temp", overwrite=overwrite_mode, quiet=quiet_mode)  #####andrea

    r.stream_order(flags="a", stream_rast=streamrast, direction=flow_dir, elevation=input_dem,
                   accumulation = flow_acc + "_temp", stream_vect=stream_vector, strahler=strah_streams,
                   overwrite=overwrite_mode, quiet=quiet_mode)
    r.out_gdal(input=strah_streams, output='%s/%s.tif' % (grass_output_db_rst, strah_streams),
               overwrite=overwrite_mode, quiet=quiet_mode)

    # Convert stream lines to vector and edit table of contents
    logging.info("--> Converting stream lines to vector and edit table of contents..")
    v.extract(input=stream_vector, type='line', output=stream_lines, quiet=quiet_mode)
    stream_lines_shpfile = '%s/%s.gpkg' % (out_tmp_dir, stream_lines)
    v.out_ogr(input=stream_lines, output=stream_lines_shpfile, format="GPKG", type="line", overwrite=overwrite_mode,
              quiet=quiet_mode)

    # Polygonize macrobasin
    mbasin_gpkgfile = '%s/v_%s_%s_mbasin.gpkg' % (out_tmp_dir, domain_name, rrs)
    r.to_vect(input=mainbasins, output='v_mbas', type='area', column='MBID', overwrite=overwrite_mode,
              quiet=quiet_mode)
    v.out_ogr(input='v_mbas', output=mbasin_gpkgfile, format="GPKG", type="area", overwrite=overwrite_mode,
              quiet=quiet_mode)

    # Identify respective MBID for each STR and join them to the attribute table of the STR ORDER
    v.overlay(ainput=stream_lines, atype='line', binput='v_mbas', btype='area', operator='and',
              output='v_bulk_ord_line_mb', overwrite=overwrite_mode, quiet=quiet_mode)

    # Clean Attribute Table
    v.db_dropcolumn(map='v_bulk_ord_line_mb', columns='a_topo_dim,a_drwal_old,a_stright,a_sinosoid,b_cat,b_cat_',
                    quiet=quiet_mode)

    # Export shapefile/GeoPackage
    v.db_renamecolum(map='v_bulk_ord_line_mb', column=('b_MBID', 'MBID'))
    v.db_renamecolum(map='v_bulk_ord_line_mb', column=('a_next_stream', 'a_next_strea'))
    v.db_renamecolum(map='v_bulk_ord_line_mb', column=('a_scheidegger', 'a_scheidegge'))
    v.db_renamecolum(map='v_bulk_ord_line_mb', column=('a_source_elev', 'a_source_ele'))
    v.db_renamecolum(map='v_bulk_ord_line_mb', column=('a_outlet_elev', 'a_outlet_ele'))

    streams_mod_gpkg_filename = 'v_%s_%s_streams_mod.gpkg' % (domain_name, rrs)
    streams_mod_gpkgfile = '%s/%s' % (out_tmp_dir, streams_mod_gpkg_filename)
    v.out_ogr(input='v_bulk_ord_line_mb', output=streams_mod_gpkgfile, format="GPKG", type="line",
              overwrite=overwrite_mode, quiet=quiet_mode)

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Trace and dissolve sub-basins

    # Polygonize subbasins
    subbasins_tmp_shpfile_filename = 'v_%s_%s_subbasins_tmp' % (domain_name, rrs)
    subbasins_tmp_shpfile = '%s/%s.shp' % (out_tmp_dir, subbasins_tmp_shpfile_filename)
    r.to_vect(input=subbasins, output='v_bas', type='area', column='stream', overwrite=overwrite_mode, quiet=quiet_mode)

    subbasins_shpfile = '%s/%s.shp' % (grass_output_db_vct, 'v_%s_%s_subbasins' % (domain_name, rrs))
    subbasins_dissolved_shpfile = '%s/%s.shp' % (grass_output_db_vct, 'v_%s_%s_subbasins_dissolved' % (domain_name, rrs))

    v.out_ogr(input='v_bas', output=subbasins_shpfile, format="ESRI_Shapefile", type="area", overwrite=overwrite_mode, quiet=quiet_mode)

    temp = gpd.read_file(subbasins_shpfile)
    temp['geometry'] = temp.buffer(0.000000001)
    temp = temp[temp.geometry.type == 'Polygon']
    basin_dissolved = temp.dissolve(by='stream', as_index=False)
    basin_dissolved.to_file(subbasins_dissolved_shpfile)
    os.system('rm ' + subbasins_shpfile.replace(".shp",".*"))

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Pfafstetter labelling
    logging.info("--> Assign pfafstetter code ..")
    streams_pfaf_lab_shp = assign_pfafstetter_code(out_tmp_dir, acp_x, target_epsg, streams_mod_gpkg_filename,
                                                       iNum_col, grass_output_db_vct)
    logging.info("--> Pfafstetter_labelling part 1 output: %s" % str(streams_pfaf_lab_shp))

    # Clean temporary files
    g.remove(flags="f", type='raster', pattern='tmp.*', quiet = quiet_mode)

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Closing algorithm
    # Estimate total execution time
    tot_elapsed_time_sec = float(time.time() - code_start_time)
    tot_elapsed_time,time_units = Give_Elapsed_Time(tot_elapsed_time_sec)	
    logging.info("REFLEX STEP 1 execution time: %.2f %s"%(tot_elapsed_time,time_units))

    # Delete GRASS temporary database including mapset
    if os.path.isdir(grass_temp_db):
        shutil.rmtree(grass_temp_db)
    else:
        pass

    # Delete output tmp sub-folder
    if os.path.isdir(out_tmp_dir):
        shutil.rmtree(out_tmp_dir)
    else:
        pass

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
    parser_handle.add_argument('-dem', action="store", dest="alg_dem")
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

    if parser_values.alg_output:
        alg_output = parser_values.alg_output
    else:
        alg_output = None

    return alg_settings, alg_log, alg_dem, alg_output
# -------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == "__main__":
    main()
# ------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------
