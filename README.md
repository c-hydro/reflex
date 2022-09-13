# REFlEx
REFlEx (Rapid Estimation of FLood EXtent) is the geomorphologic inundation model developed at CIMA Foundation

### Intallation
Prerequisites:
- Grass Gis v.7.* with addons: r.stream.basin - r.stream.order - r.accumulate
- [Compatible gdal and proj libreries](https://github.com/c-hydro/fp-system-library)
- [Python3 reflex virtual environment](https://github.com/c-hydro/reflex/tree/main/install)

WARNING! The version of the gdal python package should match the one installed on the system, before to install replace in the "fp_python3_reflex_libraries.yaml" file the version of gdal (gdal=2.3.3) with the one obtained from the command "gdal-config --version"

------------------------
#### References
*Scientific references*

Arcorace M. - [Enhancing Existing Operational Forecast-Based Flood Detection Solutions through an Integrated Use of Satellite Earth Observations and Numerical Models.](https://iris.unige.it/retrieve/handle/11567/1047031/533661/phdunige_4181596.pdf) Doctorate Thesis in Computer Science
and Systems Engineering (Cycle XXXIII). 2021, Università di Genova.

*Technical references*

Bartos M. - [Pysheds: simple and fast watershed delineation in python.](https://github.com/mdbartos/pysheds) 2020. doi: 10.5281/zenodo.3822494

GRASS Development Team - [Geographic Resources Analysis Support System (GRASS GIS) Software, Version 7.8.](http://grass.osgeo.org) Open Source Geospatial Foundation, 2019.

