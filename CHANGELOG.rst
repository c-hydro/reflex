=========
Changelog
=========
Version 2.0.2 [2022-08-11]
**************************
**reflex_step2_static_preprocessing.py**
    - Fixed procedure for concentration time estimation
**reflex_step3_hand.py**
    - Review of the multiprocessing routine
    - Add check for missing hand maps
STEP 4: **reflex_step4_flood_mapping.py**
    - Review of the multiprocessing and of the merging routine
Version 2.0.0 [2022-04-06]
**************************
GLOBAL: **revision**
        - Adopted json settings file
        - Revised parsing scheme for inputs
        - Revised logging and exception management

STEP 0: **reflex_step0_dem_conditioning.py**
        - Added simplified carving 'A'
        - Reduced number of produced outputs, export in Float32 set-up as automatic
        - Removed taudem dependency
        - Removed proj automatic fix if missing in input file
        
STEP 1: **reflex_step1_hydro_derivatives.py**
        - Integrated vector network and pfafstetter assignation for future integration of network editing
        - Removed taudem dependency an pyshed dinf algorithm adopted
        - Temporary removed MFD support
        - External pre-processed dem can be provided
        
STEP 2: **reflex_step2_static_preprocessing.py**
        - Mask file produced as unique shapefile for each stream
        - Basin masks limited to max 2 basins upstrean and max 2 downstream
        - Revised stream management (volume extimation moved to STEP 4 for a pure "static" script)
        - Automatic selection of best epsg for proj
        - Parallel computation implementation
        
STEP 3: **reflex_step3_hand.py**
        - Removed taudem and grass dependencies, pure pyheds approach
        - Parallel computation implementation
        - Updated procedure for coastal expansion - extented hand maps merging
        - Head loss now provided in cm/m (previously it was cm)

STEP 4: **reflex_step4_flood_mapping.py**
        - Parallel computation implementation
        
OTHER: **libs**
        - Review of the program libraries system

Version 1.0.0 [2019-02-20]
**************************
        - Official first release
