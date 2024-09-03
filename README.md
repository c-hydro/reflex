# REFlEx
REFlEx (Rapid Estimation of FLood EXtent) is the geomorphologic inundation model developed at CIMA Foundation.

![REFLEX_scheme_Arcorace](https://user-images.githubusercontent.com/57633516/234909358-24dd321d-ef52-4f53-ab2d-8074cc1f9fba.png)

The figure summarise the REFlEx modelling chain, outlining the interconnections between its diverse modules. The model's inputs include a digital elevation model (DEM), flowlines representing the river network, and information on the streamflow. Terrain elevation is represented by using a hydrologically conditioned DEM, which undergoes filtering, filling, and carving according to the input flowlines. The conditioned DEM is processed using the D8 and D-Infinity models to extract hydrological derivatives, such as drainage directions. To delineate floodplains, the HAND methodology is utilized, which derives local drainage potentials from normalized topography. Simultaneously, stream hierarchy is computed to develop HAND maps for each river order, starting from headwaters. Sequentially "filling" the HAND maps using the input surface runoff volumes enables the derivation of flood extent and depth information for each sub-basin. Balancing the volume underlying the HAND map with the observed one produces an optimal flood map. The output of REFLEX is a flood map that provides information on water extent and depth. 

### Intallation
Prerequisites:
- Grass Gis v.7.* with addons: r.stream.basin - r.stream.order - r.accumulate
- [Compatible gdal and proj libreries](https://github.com/c-hydro/fp-system-library)
- [Python3 reflex virtual environment](https://github.com/c-hydro/reflex/tree/main/install)

WARNING! The version of the gdal python package should match the one installed on the system, before to install replace in the "fp_python3_reflex_libraries.yaml" file the version of gdal (gdal=2.3.3) with the one obtained from the command "gdal-config --version"

------------------------
#### References
*Scientific references*

Arcorace, M., Libertino, A., Alfieri, L., Gabellani, S., Matanò, A., Masoero, A., Basso, V., & Boni, G. (2024). REFLEX—A novel method for the rapid estimation of flood extent. Journal of Flood Risk Management, e13034. DOI: [10.1111/jfr3.13034](https://doi.org/10.1111/jfr3.13034).

Arcorace M. - [Enhancing Existing Operational Forecast-Based Flood Detection Solutions through an Integrated Use of Satellite Earth Observations and Numerical Models.](https://iris.unige.it/retrieve/handle/11567/1047031/533661/phdunige_4181596.pdf) Doctorate Thesis in Computer Science
and Systems Engineering (Cycle XXXIII). 2021, Università di Genova.

*Technical references*

Bartos M. - [Pysheds: simple and fast watershed delineation in python.](https://github.com/mdbartos/pysheds) 2020. doi: 10.5281/zenodo.3822494

GRASS Development Team - [Geographic Resources Analysis Support System (GRASS GIS) Software, Version 7.8.](http://grass.osgeo.org) Open Source Geospatial Foundation, 2019.

