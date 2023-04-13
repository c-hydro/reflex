#!/usr/bin/env python
"""
REFlEx - Step1 - Hydrogical Derivatives
__date__ = '20230330'
__version__ = '2.1.0'
__author__ =
        'Mauro Arcorace' (mauro.arcorace@cimafoundation.org',
        'Alessandro Masoero (alessandro.masoero@cimafoundation.org',
        'Valerio Basso',
        'Alessia Matanò',
        'Andrea Libertino (andrea.libertino@cimafoundation.org'
__library__ = 'REFlEx'
General command line:
### python reflex_step1_hydro_derivatives.py -log_file "/path/to/log.txt" -settings_file "settings.json" -base_path "/path/to/base_folder" (-dem "/path/to/dem.tif")
Version(s):
20190220 (1.0.0) --> Beta release
20220406 (2.0.0) --> Full revision - Integrated vector network and pfafstetter assignation for future integration of network editing
                                     Temporary removed MFD support
                                     Adopted pyshed dinf algorithms
20220726 (2.0.2) --> Fix basin delineation procedure
20230330 (2.1.0) --> Optimize size and number of produced hydroderivatives
                     Simplify integration with local grass and proj installations
                     Break backward compatibility with old static data
"""
# -------------------------------------------------------------------------------------
# Import python libraries
import os
import numpy as np
import shutil
import sys
import time
from argparse import ArgumentParser
import  logging
from lib.reflex_tools_utils import Give_Elapsed_Time, set_logging, read_file_json
from lib.reflex_tools_pfafstetter import assign_pfafstetter_code
from pysheds.grid import Grid
import geopandas as gpd

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_name = 'REFlEx - STEP 1 - Hydrological Derivatives'
alg_version = '2.1.0'
alg_release = '2023-03-30'
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
    # Get parameters from config file

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

    try:
        proj_lib = data_settings["algorithm"]["proj_db"]
    except:
        proj_lib = "/usr/share/proj"

    drain_method_streams = data_settings["step_1"]["stream_definition_method"]
    drain_method_hand = data_settings["step_1"]["hand_definition_method"]
    facc_threshold = data_settings["step_1"]["flow_accumulation_threshold_cells"]
    iNum_col = 23               # Number of column of pffstetter in the shapefile

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
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
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Define variables names

    cond_dem = 'r_%s_%s_cond_dem' % (domain_name, rrs)
    filled_dem = 'r_%s_%s_filled_dem' % (domain_name, rrs)

    if (drain_method_streams == 'd8'):
        flow_dir = 'r_%s_%s_d8_flow_dir' % (domain_name, rrs)
        flow_dir4hand = 'r_%s_%s_d8_flow_dir4hand' % (domain_name, rrs)
        flow_acc = 'r_%s_%s_d8_flowacc' % (domain_name, rrs)
        flow_acc_skm = 'r_%s_%s_d8_flow_acc_skm' % (domain_name, rrs)
        flow_acc_sm = 'r_%s_%s_d8_flow_acc_sm' % (domain_name, rrs)
        streamrast = 'r_%s_%s_d8_streams' % (domain_name, rrs)
        streamrast_bin = 'r_%s_%s_d8_streamrast_bin' % (domain_name, rrs)
        whatersheds = 'r_%s_%s_d8_whatersheds' % (domain_name, rrs)
        subbasins = 'r_%s_%s_d8_subbasins' % (domain_name, rrs)
        mainbasins = 'r_%s_%s_d8_mainbasins' % (domain_name, rrs)
        sinks = 'r_%s_%s_d8_sinks' % (domain_name, rrs)
        stream_distance_rst = 'r_%s_%s_d8_stream_distance' % (domain_name, rrs)
    else:
        logging.error("Failed in using drainage method %s. Not Supported." % str(drain_method_hand))
        sys.exit(1)

    if (drain_method_hand == 'dinf'):
        dinf_flow_dir = 'r_%s_%s_dinf_flow_dir' % (domain_name, rrs)
        dinf_stream_distance = 'r_%s_%s_dinf_stream_distance' % (domain_name, rrs)
    else:
        pass
    cell_area = 'r_%s_%s_areacell' % (domain_name, rrs)
    strah_streams = 'r_%s_%s_d8_strahler_streams' % (domain_name, rrs)
    stream_vector = 'v_%s_%s_d8_stream_vect' % (domain_name, rrs)
    stream_lines = 'v_%s_%s_d8_stream_lines' % (domain_name, rrs)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Create output directories

    # Define output MAIN directory
    grass_output_db = step1_dir_name
    logging.info('--> Output main directory: %s' % grass_output_db)
    if not os.path.exists(grass_output_db):
        os.makedirs(grass_output_db)

    # Define output tmp sub-directory
    out_tmp_dirname = 'tmp'
    out_tmp_dir = os.path.join(grass_output_db, out_tmp_dirname)
    logging.info('--> Output tmp sub-directory: %s' % out_tmp_dir)
    if not os.path.exists(out_tmp_dir):
        os.makedirs(out_tmp_dir)

    # Define output MAIN directory
    logging.info('--> Output main directory: %s' % grass_output_db)
    if not os.path.exists(grass_output_db):
        os.makedirs(grass_output_db)

    # Define output raster sub-directory
    grass_output_db_rst_dirname = 'rst'
    grass_output_db_rst = os.path.join(grass_output_db, grass_output_db_rst_dirname)
    logging.info('--> Output raster sub-directory: %s' % grass_output_db_rst)
    if not os.path.exists(grass_output_db_rst):
        os.makedirs(grass_output_db_rst)

    # Define output vector sub-directory
    grass_output_db_vct_dirname = 'vct'
    grass_output_db_vct = os.path.join(grass_output_db, grass_output_db_vct_dirname)
    logging.info('--> Output vector sub-directory: %s' % grass_output_db_vct)
    if not os.path.exists(grass_output_db_vct):
        os.makedirs(grass_output_db_vct)
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Initialise system

    # Record starting simulation time
    code_start_time = time.time()

    # Start GRASS
    logging.info("--> Starting GRASS..")

    os.environ['PROJ_LIB'] = proj_lib
    if not os.path.join(proj_lib,"proj.db"):
        logging.error(" ERROR! proj.db file not found at " + proj_lib + " Please add the 'proj_db' key under the 'algorithm' key of the setting file to specify the location!")
        raise FileNotFoundError
    from grass_session import TmpSession

    # Configure grass
    os.environ['GISDBASE'] = os.path.join(step1_dir_name, "grass_temp", "")
    os.environ['GRASSBIN'] = grass_bin
    os.environ['GISBASE'] = os.popen(grass_bin + ' --config path').read().strip('\n')
    gpydir = os.path.join(os.environ['GISBASE'], "etc", "python")
    sys.path.append(gpydir)
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.pygrass.modules.shortcuts import vector as v
    from grass.pygrass.modules.shortcuts import database as db
    from grass.pygrass.modules import Module
    from grass.exceptions import GrassError

    if os.path.isdir(os.path.join(os.environ['GISDBASE'])):
        shutil.rmtree(os.path.join(os.environ['GISDBASE']))

    with TmpSession(gisdb=os.environ['GISDBASE'], location="wgs84", create_opts='EPSG:4326'):

        # Define alias for GRASS modules
        r.out_gdal = Module('r.out.gdal')
        try:
            r.stream_basins = Module('r.stream.basins')
        except GrassError:
            logging.warning('r.stream.basins not found... INSTALL')
            g.extension(extension='r.stream.basins', operation="add")
            r.stream_basins = Module('r.stream.basins')
        try:
            r.stream_extract = Module('r.stream.extract')
        except GrassError:
            logging.warning('r.stream.extract not found... INSTALL')
            g.extension(extension='r.stream.extract', operation="add")
            r.stream_extract = Module('r.stream.extract')
        try:
            r.stream_distance = Module('r.stream.distance')
        except GrassError:
            logging.warning('r.stream.distance not found... INSTALL')
            g.extension(extension='r.stream.distance', operation="add")
            r.stream_distance = Module('r.stream.distance')
        try:
            r.stream_distance = Module('r.stream.order')
        except GrassError:
            logging.warning('r.stream.distance not found... INSTALL')
            g.extension(extension='r.stream.order', operation="add")
            r.stream_order = Module('r.stream.order')

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

        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Import dem and calculate morphology

        # Import dem
        logging.info("--> Import filled dem")
        r.external(input=dem_input_filled_file, output=filled_dem, quiet = quiet_mode)
        g.region(raster=filled_dem)

        # Import conditioned dem
        logging.info("--> Import conditioned dem")
        r.external(input=dem_input_conditioned_file, output=cond_dem, quiet = quiet_mode)

        # Calculate slope
        logging.info("--> Calculate slope")
        dem_slope = 'r_%s_%s_dem_slope_perthousend' % (domain_name, rrs)
        r.slope_aspect(elevation=filled_dem, slope=dem_slope, format='percent', overwrite=overwrite_mode)
        r.mapcalc("dem_slope_perthousend = " + dem_slope + " * 1000" , overwrite = overwrite_mode, quiet = quiet_mode)
        r.out_gdal(input="dem_slope_perthousend", output=os.path.join(grass_output_db_rst,dem_slope + ".tif"), createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode, type='Int32', flags="f")

        # Convert from meter to cm (for lighter computation)
        logging.info("--> Export dem in cm")
        filled_dem_cm='r_%s_%s_filled_dem_cm'%(domain_name,rrs)
        r.mapcalc("%s=int(%s*100)"%(filled_dem_cm,filled_dem), overwrite = overwrite_mode, quiet = quiet_mode)
        out_input_dem_cm_file='%s/%s.tif'%(grass_output_db_rst,filled_dem_cm)
        r.out_gdal(input=filled_dem_cm, output=out_input_dem_cm_file, createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode, type='Int32', flags="f")
        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Extract hydro derivatives
        # Extract hydrological derivatives
        logging.info("--> Extracting hydrological derivatives from DEM ..")
        if (drain_method_streams=='d8'):
            logging.info("--> Extracting hydrological derivatives from DEM using the SFD (D8) method..")
            r.watershed(flags="sab", elevation=cond_dem, threshold='%s'%facc_threshold, accumulation=flow_acc, drainage=flow_dir, basin=whatersheds, stream=streamrast, overwrite = overwrite_mode, quiet = quiet_mode)
            r.out_gdal(input=flow_dir, output='%s/%s.tif'%(grass_output_db_rst,flow_dir), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
            r.out_gdal(input=streamrast, output='%s/%s.tif'%(grass_output_db_rst,streamrast), overwrite = overwrite_mode, quiet = quiet_mode)
            r.out_gdal(input=flow_acc, output='%s/%s.tif' % (grass_output_db_rst, flow_acc), createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode)
            r.out_gdal(input=whatersheds, output='%s/%s.tif' % (grass_output_db_rst, whatersheds), createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode, type='Int32', flags='f')
        else:
            logging.error("Failed in using drainage method %s. Not Supported."%str(drain_method_streams))
            sys.exit(1)

        logging.info("--> Calculate slope channel ..")
        r.mapcalc("slope_channel=if(isnull(" + streamrast + "),null(),dem_slope_perthousend)", overwrite=overwrite_mode, quiet=quiet_mode)
        r.out_gdal(input="slope_channel", output=os.path.join(grass_output_db_rst, "r_{domain_name}_{rrs}_slope_channel_perthousend.tif".format(domain_name=domain_name, rrs=rrs)),
                   createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, quiet=quiet_mode, type='Int32', flags="f")

        # Convert pixels count to areas using flow accumulation
        logging.info("--> Convert flow accumulation pixels count to areas..")
        r.mapcalc(flow_acc_skm + "=" + flow_acc + "*" + str(rrs_km) + "*" + str(rrs_km), overwrite=overwrite_mode, quiet=quiet_mode)
        r.out_gdal(input=flow_acc_skm, output='%s/%s.tif' % (grass_output_db_rst, flow_acc_skm), createopt="COMPRESS=DEFLATE", overwrite=overwrite_mode, type='Float32', flags="f")

        # Delineate basins and subbasins
        logging.info("--> Delineating macro basins..")
        r.stream_basins(flags="l", direction=flow_dir, stream_rast=streamrast, memory=300, basins=mainbasins, overwrite=overwrite_mode, quiet=quiet_mode)
        r.out_gdal(input=mainbasins, output='%s/%s.tif'%(grass_output_db_rst,mainbasins), createopt="COMPRESS=DEFLATE", overwrite = overwrite_mode, quiet = quiet_mode)
        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Calculate d-inf derivatives
        if (drain_method_hand=='dinf'):

            # Compute D-infinity flow direction:
            logging.info("--> Computing D-infinity Flow Direction..")

            grid = Grid.from_raster(dem_input_filled_file)
            dem = grid.read_raster(dem_input_filled_file, nodata=-9999)
            fdir = grid.flowdir(dem, routing='dinf')

            dinf_flow_dir_outfile = '%s/%s.tif' % (grass_output_db_rst, dinf_flow_dir)

            grid.to_raster(fdir, dinf_flow_dir_outfile, nodata=-1, dtype=np.float32)
        else:
            pass

        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Calculate stream hierarchy

        # Calculates Strahler's streams hierarchy
        logging.info("--> Calculating streams hierarchy..")
        # Replace old flow dir with the one generated with r.accumulated
        r.accumulate(direction=flow_dir, accumulation=flow_acc + "_temp", overwrite=overwrite_mode, quiet=quiet_mode)

        r.stream_order(flags="a", stream_rast=streamrast, direction=flow_dir, elevation=cond_dem,
                       accumulation = flow_acc + "_temp", stream_vect=stream_vector, strahler=strah_streams,
                       overwrite=overwrite_mode, quiet=quiet_mode)
        r.out_gdal(input=strah_streams, output='%s/%s.tif' % (grass_output_db_rst, strah_streams),
                   overwrite=overwrite_mode, quiet=quiet_mode)

        # Convert stream lines to vector and edit table of contents
        logging.info("--> Converting stream lines to vector and edit table of contents..")
        v.extract(input=stream_vector, type='line', output=stream_lines, quiet=quiet_mode)
        stream_lines_shpfile = '%s/%s.gpkg' % (out_tmp_dir, stream_lines)
        v.out_ogr(input=stream_lines, output=stream_lines_shpfile, format="GPKG", type="line", overwrite=overwrite_mode, quiet=quiet_mode)

        # Polygonize macrobasin
        mbasin_gpkgfile = '%s/v_%s_%s_mbasin.gpkg' % (out_tmp_dir, domain_name, rrs)
        r.to_vect(input=mainbasins, output='v_mbas', type='area', column='MBID', overwrite=overwrite_mode, quiet=quiet_mode)
        v.out_ogr(input='v_mbas', output=mbasin_gpkgfile, format="GPKG", type="area", overwrite=overwrite_mode, quiet=quiet_mode)

        # Identify respective MBID for each STR and join them to the attribute table of the STR ORDER
        v.overlay(ainput=stream_lines, atype='line', binput='v_mbas', btype='area', operator='and', output='v_bulk_ord_line_mb', overwrite=overwrite_mode, quiet=quiet_mode)

        # Clean Attribute Table
        v.db_dropcolumn(map='v_bulk_ord_line_mb', columns='a_topo_dim,a_drwal_old,a_stright,a_sinosoid,b_cat,b_cat_', quiet=quiet_mode)

        # Export shapefile/GeoPackage
        v.db_renamecolum(map='v_bulk_ord_line_mb', column=('b_MBID', 'MBID'))
        v.db_renamecolum(map='v_bulk_ord_line_mb', column=('a_next_stream', 'a_next_strea'))
        v.db_renamecolum(map='v_bulk_ord_line_mb', column=('a_scheidegger', 'a_scheidegge'))
        v.db_renamecolum(map='v_bulk_ord_line_mb', column=('a_source_elev', 'a_source_ele'))
        v.db_renamecolum(map='v_bulk_ord_line_mb', column=('a_outlet_elev', 'a_outlet_ele'))

        streams_mod_gpkg_filename = 'v_%s_%s_streams_mod.gpkg' % (domain_name, rrs)
        streams_mod_gpkgfile = '%s/%s' % (out_tmp_dir, streams_mod_gpkg_filename)
        v.out_ogr(input='v_bulk_ord_line_mb', output=streams_mod_gpkgfile, format="GPKG", type="line",  overwrite=overwrite_mode, quiet=quiet_mode)

        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Trace and dissolve sub-basins
        # Polygonize subbasins
        r.to_vect(input=whatersheds, output='v_bas', type='area', column='stream', overwrite=overwrite_mode, quiet=quiet_mode)
        subbasins_shpfile = '%s/%s.shp' % (grass_output_db_vct, 'v_%s_%s_subbasins' % (domain_name, rrs))
        subbasins_dissolved_shpfile = '%s/%s.shp' % (grass_output_db_vct, 'v_%s_%s_subbasins' % (domain_name, rrs))
        v.out_ogr(input='v_bas', output=subbasins_shpfile, format="ESRI_Shapefile", type="area", overwrite=overwrite_mode, quiet=quiet_mode)

        temp = gpd.read_file(subbasins_shpfile)
        temp['geometry'] = temp.buffer(0.000000001)
        temp = temp[temp.geometry.type == 'Polygon']
        basin_dissolved = temp.dissolve(by='stream', as_index=False)
        basin_dissolved.to_file(subbasins_dissolved_shpfile)
        # -------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------
        # Pfafstetter labelling
        logging.info("--> Assign pfafstetter code ..")
        streams_pfaf_lab_shp = assign_pfafstetter_code(out_tmp_dir, acp_x, target_epsg, streams_mod_gpkg_filename, iNum_col, grass_output_db_vct)
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