# Eventually source to gdal-proj envs
# source /home/reflex/fp_system_libs_generic/fp_system_libs_generic_proj

REFLEX_PATH='/home/reflex/fp_apps_system/reflex/'
WORKING_PATH="/home/reflex/data/ISMEA/reflex/magra/"
LOG_PATH=$WORKING_PATH'/log/'
CONFIG_FILE=$WORKING_PATH"/magra.json"
DEM_FILE=$WORKING_PATH"/input/DEM_magra.tif"
STREAM_FILE=$WORKING_PATH"/input/network_magra.shp"
mkdir -p $LOG_PATH

#STEP 0
python3 $REFLEX_PATH/reflex_step0_dem_conditioning.py -log_file $LOG_PATH/log0.txt -settings_file $CONFIG_FILE -dem $DEM_FILE -stream $STREAM_FILE -base_path $WORKING_PATH

#STEP 1
python3 $REFLEX_PATH/reflex_step1_hydro_derivatives.py -log_file $LOG_PATH/log1.txt -settings_file $CONFIG_FILE -base_path $WORKING_PATH

#STEP 2
python3 $REFLEX_PATH/reflex_step2_static_preprocessing.py -log_file $LOG_PATH/log2.txt -settings_file $CONFIG_FILE -base_path $WORKING_PATH

#STEP3
python3 $REFLEX_PATH/reflex_step3_hand.py -log_file $LOG_PATH/log3.txt -settings_file $CONFIG_FILE -base_path $WORKING_PATH

#STEP4
python3 $REFLEX_PATH/reflex_step4_flood_mapping.py -log_file $LOG_PATH/log4.txt -settings_file $CONFIG_FILE -base_path $WORKING_PATH
