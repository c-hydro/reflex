**Installation**

Prerequisites:
- Grass Gis v.7.*
- Grass Gis addons: r.stream.basin - r.stream.order - r.accumulate
- Gdal and proj libreries (git@github.com:c-hydro/fp-system-library.git)
- Python3 reflex virtual environment (use the provided setup_fp_system_conda_python_reflex.sh)

WARNING! The version of the gdal python package should match the one installed on the system, before to install replace in the "fp_python3_reflex_libraries.yaml" file the version of gdal (gdal=2.3.3) with the one obtained from the command "gdal-config --version"
