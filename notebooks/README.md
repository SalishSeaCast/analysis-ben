The Jupyter Notebooks in this directory are made by Ben for
quick sharing of results.

The links below are to static renderings of the notebooks via
[nbviewer.ipython.org](http://nbviewer.ipython.org/).
Descriptions below the links are from the first cell of the notebooks
(if that cell contains Markdown or raw text).

* ##[timeseries_tools_dev.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/timeseries_tools_dev.ipynb)  
    
    **`timeseries_tools` Development Notebook**  
    This notebook is for developing a prototype Nowcast timeseries analysis package. The primary goal of this package is to be memory-efficient. As such, the results arrays are flattened to 2-D (time, space) so that land points can be removed. The 2-D dimensions are also ideal for some analyses like PCA. The basic workflow proceeds as follows:  
       * Flatten the model grid and mask to 2-D (time, space) and remove land indices  
       * Load, process, and flatten hourly Nowcast Results to 2-D (time, space) and remove land indices  
       * Concatenate consecutive 24 hour periods  
       * Reshape the concatenated `Numpy ndarray` to (time, depth, y, x)  

* ##[Hakai_data.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/Hakai_data.ipynb)  
    
* ##[HRDPS_loading_numpy.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/HRDPS_loading_numpy.ipynb)  
    
    **HRDPS loading in Numpy**  

* ##[timeseries_tools_usage.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/timeseries_tools_usage.ipynb)  
    
    **Timeseries Tools Usage**  
    This notebook demonstrates loading in a Nowcast timeseries using `timeseries_tools`. `Timeseries_tools` extracts only the required fields from the Nowcast netCDF files and concatenates the flattened (time, space) numpy arrays together. This uses a minimum amount of memory relative to bulkier routines such as `mfdataset`. The flattened arrays are reshaped before they are returned, but there is an option to keep them flattened, which is necessary for certain data analysis methods like principal component analysis.  

* ##[maps_cartopy.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/maps_cartopy.ipynb)  
    
    **Maps in Cartopy**  

* ##[WindAnalysis.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/WindAnalysis.ipynb)  
    
* ##[nowcast_analysis.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/nowcast_analysis.ipynb)  
    
    **Nowcast Analysis Notebook**  
    This notebook demonstrates the process of loading and analyzing long (6+ month) timeseries of SalishSeaCast Nowcast (or Hindcast) results.  
      
    The analysis is organized into the following sections:  
      
       1. [Load Results](#Load-Results)  
       2. [Postprocessing](#Postprocessing)  
       3. [Wind Averaging](#Wind-Averaging)  
       4. Spectral Coherence  
       5. [Principal Component Analysis](#Principal-Component-Analysis)  
       6. Canonical Correlation Analysis  

* ##[OceanTimeseries.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/OceanTimeseries.ipynb)  
    
    **WOA data**  
      
    Temperature and nitrate climatologies from  
      
    [https://www.nodc.noaa.gov/OC5/woa18/woa18data.html](https://www.nodc.noaa.gov/OC5/woa18/woa18data.html)  
      
    ***  

* ##[Tug_spill_preliminary.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/Tug_spill_preliminary.ipynb)  
    
* ##[GRL2016_data.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/GRL2016_data.ipynb)  
    
* ##[master_hindcast_extractor.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/master_hindcast_extractor.ipynb)  
    
    **Master hindcast extractor**  
      
    My latest code for optimized extraction of long hindcast fields  
      
    The local extraction method uses [Dask](https://dask.org) and [xarray.open_mfdataset](http://xarray.pydata.org/en/stable/io.html#reading-multi-file-datasets) based on Doug's optimized workflow in [dask-expts](https://nbviewer.jupyter.org/urls/bitbucket.org/salishsea/analysis-doug/raw/default/notebooks/dask-expts/dask_expts.ipynb). \  
    *I get about 1 h extraction time per year of hourly surface fields across 12 Salish workers using this method.*  
      
    The ERDDAP method is much slower but publically accessible and requires no user-side CPU resources. \  
    *Extraction time is about 12 h per year of hourly surface fields.*  
      
    Both methods extract in monthy chunks and save to yearly files. The code is designed to combine variables on the *same* model grid (temperature and nitrate for example) but not across different grids (so velocity separate from temperate separate from wind fields, run multiple times to extract these separately).  
      
    ***  

* ##[MPI_decomposition_calculator.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/MPI_decomposition_calculator.ipynb)  
    
    **MPI decomposition calculator**  
    For calculating decompositions equal to whole nodes  

* ##[matrix_interpolation.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/matrix_interpolation.ipynb)  
    
    **Interpolation matrices for OpenDrift**  
      
    ***  

* ##[wind_decomposition.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/wind_decomposition.ipynb)  
    
    **SoG wind decomposition recipes**  
      
    ***  

* ##[NEMOConfigSetup.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/NEMOConfigSetup.ipynb)  
    
* ##[maps_for_Genna.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/maps_for_Genna.ipynb)  
    
    **Maps for Genna**  

* ##[CPCA_varimax.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/CPCA_varimax.ipynb)  
    
    **Principal Component Analysis of Nowcast Velocities**  
    This notebook describes the process of performing complex principal component analysis (C-PCA) on the SalishSeaCast Nowcast velocity record and explores the results.  
      
    This notebook is organized into the following sections:  
      
       1. [Load Nowcast and GEM Record](#Load-Nowcast-and-GEM-Record)  
       2. [Principal Component Analysis](#Principal-Component-Analysis)  
       3. [Varimax Rotation](#Varimax-Rotation)  
       4. [Explore Results](#2016-HRDPS-Wind-EOFs)  

* ##[ONC_API_usage.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/ONC_API_usage.ipynb)  
    
    **ONC API 2.0 Usage**  
    Python API usage  
      
    [`https://wiki.oceannetworks.ca/display/O2A/Python+Client+Library`](https://wiki.oceannetworks.ca/display/O2A/Python+Client+Library)  
      
    Timeseries parameters  
      
    [`https://wiki.oceannetworks.ca/display/DP/1`](https://wiki.oceannetworks.ca/display/DP/1)  
      
    ADCP parameters  
      
    [`https://wiki.oceannetworks.ca/display/DP/5`](https://wiki.oceannetworks.ca/display/DP/5)  
      
    ***  
      
    **Installation**  
      
    `$ pip install onc`  
      
    ***  
      
    **Import libraries**  

* ##[make_nowcast_timeseries.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/make_nowcast_timeseries.ipynb)  
    
    **Make Nowcast Timeseries**  
    This notebook describes the process of extracting a timeseries of SalishSeaCast Nowcast and HRDPS results for analysis.  

* ##[Maps.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/Maps.ipynb)  
    
    **Salish Sea Maps**  
    This notebook plots various maps relevant to the SalishSeaCast project. The first map is a general study area map. The second map shows the model domain and bathymetry. The third map shows the Strait of Georgia with relevant observation platform locations overplotted.  

* ##[gemlam_weirdness.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/gemlam_weirdness.ipynb)  
    
    **Gemlam weirdness**  
      
    I was messing around with some of the early HRDPS record (gemlam) and noticed some inconsistencies in the gridding. I've included a few example dates below.  
      
    ***  

* ##[NEMO_grid_rotation.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/NEMO_grid_rotation.ipynb)  
    
* ##[visualization_workflows_xarray.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/visualization_workflows_xarray.ipynb)  
    
    **Visualisation of NEMO/GEM/Observations using XArray**  
    This notebook demonstrates the use of several tools for easy loading and visualization of model results and drifter observations. Most of these tools require `xarray`, and rather than make them flexible to other packages, I thought I would just demonstrate their use with `xarray` here.  
      
    First we'll import the necessary libraries and set our preferred notebook formatting.  

* ##[SalishSeaCast_timeseries_from_Skookum.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/SalishSeaCast_timeseries_from_Skookum.ipynb)  
    
    **Working with daily SalishSeaCast results from Skookum**  

* ##[OS2018_plots.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/OS2018_plots.ipynb)  
    
    **Hindcast timeseries analysis**  

* ##[working_with_HRDPS.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/working_with_HRDPS.ipynb)  
    
    **Working with HRDPS**  

* ##[HRDPS_correction.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/HRDPS_correction.ipynb)  
    
    **HRDPS Correction Check**  
    This notebook checks the HRDPS correction obtained from Luc Fillion's group  

* ##[SedimentModel.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/SedimentModel.ipynb)  
    

##License

These notebooks and files are copyright 2013-2020
by the Salish Sea MEOPAR Project Contributors
and The University of British Columbia.

They are licensed under the Apache License, Version 2.0.
http://www.apache.org/licenses/LICENSE-2.0
Please see the LICENSE file for details of the license.
