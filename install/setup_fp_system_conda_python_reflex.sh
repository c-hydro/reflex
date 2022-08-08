#!/bin/bash -e

#-----------------------------------------------------------------------------------------
# Script information
script_name='FP ENVIRONMENT - PYTHON3 LIBRARIES FOR HMC PACKAGE - CONDA'
script_version="1.6.3"
script_date='2022/05/23'

# Define file reference path according with https link(s) --> https://repo.anaconda.com/miniconda/
fp_env_file_miniconda='https://repo.continuum.io/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh'

# Argument(s) default definition(s)
fp_env_folder_root_default=$HOME/fp_envs_python3
fp_env_file_reference_default='fp_system_conda_python3_reflex_settings'
fp_env_folder_libraries_default='fp_system_conda_python3_reflex_libraries'

fp_env_file_requirements_default='requirements_fp_system_conda_python3_reflex.yaml'

# Examples of generic command-line:
# conda create --yes --name $fp_env_folder_libraries numpy scipy pip python=3
# conda create --yes --nama $fp_env_folder_libraries -c conda-forge pyresample pygeobase
# conda install --yes -c conda-forge nbconvert
# conda env export | grep -v "^prefix: " > requirements_fp_env_python.yaml
# conda env create --file $fp_env_file_requirements # create virtual env using yaml file
# conda create -y --name $fp_env_folder_libraries --file $fp_env_file_requirements # create virtual env using ascii file
#-----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Info script start
echo " ==================================================================================="
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> START ..."

# Get arguments number and values
script_args_n=$#
script_args_values=$@

echo ""
echo " ==> Script arguments number: $script_args_n"
echo " ==> Script arguments values: $script_args_values"
echo ""
echo " ==> Script arguments 1 - Directory of libraries [string: path]-> $1"
echo " ==> Script arguments 2 - Filename of system environment [string: filename] -> $2"
echo " ==> Script arguments 3 - Name of virtual environment [string: name] -> $3"
echo " ==> Script arguments 4 - Filename of system requirements in yaml format [string: name] -> $4"
echo ""

# Get folder root path
if [ $# -eq 0 ]; then
    fp_env_folder_root=$fp_env_folder_root_default		
	fp_env_file_reference=$fp_env_file_reference_default
	fp_env_folder_libraries=$fp_env_folder_libraries_default
	fp_env_file_requirements=$fp_env_file_requirements_default
elif [ $# -eq 1 ]; then
	fp_env_folder_root=$1	
	fp_env_file_reference=$fp_env_file_reference_default
	fp_env_folder_libraries=$fp_env_folder_libraries_default
	fp_env_file_requirements=$fp_env_file_requirements_default
elif [ $# -eq 2 ]; then
	fp_env_folder_root=$1
	fp_env_file_reference=$2
	fp_env_folder_libraries=$fp_env_folder_libraries_default
	fp_env_file_requirements=$fp_env_file_requirements_default
elif [ $# -eq 3 ]; then
	fp_env_folder_root=$1
	fp_env_file_reference=$2
	fp_env_folder_libraries=$3
	fp_env_file_requirements=$fp_env_file_requirements_default
elif [ $# -eq 4 ]; then
	fp_env_folder_root=$1
	fp_env_file_reference=$2
	fp_env_folder_libraries=$3
	fp_env_file_requirements=$4
fi

# Create root folder
if [ ! -d "$fp_env_folder_root" ]; then
	mkdir -p $fp_env_folder_root
fi

# Define folder path(s)
fp_env_folder_root=$fp_env_folder_root

# Define environment filename
fp_env_path_reference=$fp_env_folder_root/$fp_env_file_reference

# multilines comment: if [ 1 -eq 0 ]; then ... fi
#-----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Install python environmente using miniconda
echo " ====> CHECK PYTHON ENVIRONMENT ... "
if [ -d "$fp_env_folder_root/bin/" ] > /dev/null 2>&1 ; then
	export PATH="$fp_env_folder_root/bin:$PATH"
	if source activate $fp_env_folder_root > /dev/null 2>&1 ; then
		echo " ====> CHECK PYTHON ENVIRONMENT ... FOUND."
    	fp_env_install=false
    else
   	 	echo " ====> CHECK PYTHON ENVIRONMENT ... FOUND ENVIRONMENT BUT LIBRARIES ARE NOT FOUND."
    	fp_env_install=false
    fi
else
	echo " ====> CHECK PYTHON ENVIRONMENT ... NOT FOUND."
	fp_env_install=true
fi

echo " ====> INSTALL PYTHON ENVIRONMENT ... "
if $fp_env_install; then
	# Download library source codes
	echo " =====> GET LIBRARY FILES ... "
	wget $fp_env_file_miniconda -O miniconda.sh
	echo " =====> GET LIBRARY FILES ... DONE!"

	if [ -d "$fp_env_folder_root" ]; then 
		rm -Rf $fp_env_folder_root; 
	fi

	bash miniconda.sh -b -p $fp_env_folder_root
	
	echo " ====> INSTALL PYTHON ENVIRONMENT ... DONE!"
	
else
	echo " ====> INSTALL PYTHON ENVIRONMENT ... DONE. PREVIOSLY INSTALLED"
fi
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Install python libraries
echo " ====> INSTALL PYTHON LIBRARIES ... "
export PATH="$fp_env_folder_root/bin:$PATH"

if [ -f $fp_env_file_requirements ] ; then
	echo " =====> USE OF CONDA REQUIREMENTS FILE YAML: $fp_env_file_requirements"
	# source activate $fp_env_folder_libraries
    # conda create -y --name $fp_env_folder_libraries --file $fp_env_file_requirements
    conda env create --file $fp_env_file_requirements
   
else

	echo " =====> USE OF CONDA GENERIC COMMAND-LINE"
	
	echo " =====> [1/2] CONDA-DEFAULT CHANNEL INSTALLATION ... "
	conda create --yes --name $fp_env_folder_libraries numpy scipy pandas matplotlib=3.1.3 rasterio geopandas cartopy=0.17 netCDF4 cython h5py proj4 xarray bottleneck dask pip python=3.7
	echo " =====> [1/2] CONDA-DEFAULT CHANNEL INSTALLATION ... DONE"
	
	echo " =====> [2/2] PYTHON-PIP INSTALLATION ... "
	source activate $fp_env_folder_libraries
	pip install h5netcdf
	echo " =====> [2/2] PYTHON-PIP INSTALLATION ... DONE"

fi

echo " ====> INSTALL PYTHON LIBRARIES ... DONE!"
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Create environmental file
echo " ====> CREATE ENVIRONMENTAL FILE ... "

# Delete old version of environmental file
cd $fp_env_folder_root

if [ -f $fp_env_file_reference ] ; then
    rm $fp_env_file_reference
fi

# Export BINARY PATH(S)
echo "PATH=$fp_env_folder_root/bin:"'$PATH'"" >> $fp_env_file_reference
echo "export PATH" >> $fp_env_file_reference

# Export VENV ACTIVATION
echo "source activate $fp_env_folder_libraries" >> $fp_env_file_reference

echo " ====> CREATE ENVIRONMENTAL FILE ... DONE!"
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Info script end
echo " ==> "$script_name" (Version: "$script_version" Release_Date: "$script_date")"
echo " ==> ... END"
echo " ==> Bye, Bye"
echo " ==================================================================================="
# ----------------------------------------------------------------------------------------

