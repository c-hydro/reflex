=========
Changelog
=========
Version 2.1.0 [2023-03-30]
**************************
**CAUTION**
Compatibility with previously created static data is broken when upgrading to versions above 2.1.x

**reflex_step0_dem_conditioning.py**
    - Simplify integration with local grass an proj installations

**reflex_step1_hydro_derivatives.py**
    - Simplify integration with local grass an proj installations
    - Rationalise the number of produced hydroderivatives
    - Dem is exported in cm and slope in per-thousend to use integer maps

**reflex_step2_static_preprocessing.py**
    - Optimize singleprocessing for big domains
    - Revise all the process, merged basin features scripts

**reflex_step3_hand.py**
    - Optimize singleprocessing for big domains

**reflex_step4_flood_mapping.py**
    - Add possibility to provide a roughness coefficient for each stream
    - Optimize singleprocessing for big domains

**reflex_tools_basins.py**
    - Merged basin masks and basin features functions for optimize files IO

**reflex_tools_flooding.py**
    - Pick the stream-specific roughness or rely on the provided one

Version 2.0.2 [2022-10-25]
**************************
**reflex_step2_static_preprocessing.py**
    - Bug fixes

**reflex_step4_flood_mapping.py**
    - Added the possibility to manually activate coastal expansion

Version 2.0.1 [2022-08-11]
**************************
**reflex_step2_static_preprocessing.py**
    - Fixed procedure for concentration time estimation

**reflex_step3_hand.py**
    - Review of the multiprocessing routine
    - Add check for missing hand maps

**reflex_step4_flood_mapping.py**
    - Review of the multiprocessing and of the merging routine
    - Manage grass bug that set to 0 gradients lower than 1E-06

Version 2.0.0 [2022-04-06]
**************************
GLOBAL: **revision**
        - Adopted json settings file
        - Revised parsing scheme for inputs
        - Revised logging and exception management

**reflex_step0_dem_conditioning.py**
        - Added simplified carving 'A'
        - Reduced number of produced outputs, export in Float32 set-up as automatic
        - Removed taudem dependency
        - Removed proj automatic fix if missing in input file
        
**reflex_step1_hydro_derivatives.py**
        - Integrated vector network and pfafstetter assignation for future integration of network editing
        - Removed taudem dependency an pyshed dinf algorithm adopted
        - Temporary removed MFD support
        - External pre-processed dem can be provided
        
**reflex_step2_static_preprocessing.py**
        - Mask file produced as unique shapefile for each stream
        - Basin masks limited to max 2 basins upstrean and max 2 downstream
        - Revised stream management (volume extimation moved to STEP 4 for a pure "static" script)
        - Automatic selection of best epsg for proj
        - Parallel computation implementation
        
**reflex_step3_hand.py**
        - Removed taudem and grass dependencies, pure pyheds approach
        - Parallel computation implementation
        - Updated procedure for coastal expansion - extented hand maps merging
        - Head loss now provided in cm/m (previously it was cm)

**reflex_step4_flood_mapping.py**
        - Parallel computation implementation
        
**libs**
        - Review of the program libraries system

Version 1.0.0 [2019-02-20]
**************************
        - Official first release
