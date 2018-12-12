The Jupyter Notebooks in this directory are made by Ben for
quick sharing of results.

The links below are to static renderings of the notebooks via
[nbviewer.ipython.org](http://nbviewer.ipython.org/).
Descriptions below the links are from the first cell of the notebooks
(if that cell contains Markdown or raw text).

* ##[Hakai_data.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Hakai_data.ipynb)  
    
* ##[visualization_workflows_xarray.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/visualization_workflows_xarray.ipynb)  
    
    **Visualisation of NEMO/GEM/Observations using XArray**  
    This notebook demonstrates the use of several tools for easy loading and visualization of model results and drifter observations. Most of these tools require `xarray`, and rather than make them flexible to other packages, I thought I would just demonstrate their use with `xarray` here.  
      
    First we'll import the necessary libraries and set our preferred notebook formatting.  

* ##[HRDPS_correction.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/HRDPS_correction.ipynb)  
    
    **HRDPS Correction Check**  
    This notebook checks the HRDPS correction obtained from Luc Fillion's group  

* ##[timeseries_tools_usage.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/timeseries_tools_usage.ipynb)  
    
    **Timeseries Tools Usage**  
    This notebook demonstrates loading in a Nowcast timeseries using `timeseries_tools`. `Timeseries_tools` extracts only the required fields from the Nowcast netCDF files and concatenates the flattened (time, space) numpy arrays together. This uses a minimum amount of memory relative to bulkier routines such as `mfdataset`. The flattened arrays are reshaped before they are returned, but there is an option to keep them flattened, which is necessary for certain data analysis methods like principal component analysis.  

* ##[make_nowcast_timeseries.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/make_nowcast_timeseries.ipynb)  
    
    **Make Nowcast Timeseries**  
    This notebook describes the process of extracting a timeseries of SalishSeaCast Nowcast and HRDPS results for analysis.  

* ##[SeicheScaling.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/SeicheScaling.ipynb)  
    
    **Scaling Upwelling and Seicheing in the Strait of Georgia**  

* ##[nowcast_analysis.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/nowcast_analysis.ipynb)  
    
    **Nowcast Analysis Notebook**  
    This notebook demonstrates the process of loading and analyzing long (6+ month) timeseries of SalishSeaCast Nowcast (or Hindcast) results.  
      
    The analysis is organized into the following sections:  
      
       1. [Load Results](#Load-Results)  
       2. [Postprocessing](#Postprocessing)  
       3. [Wind Averaging](#Wind-Averaging)  
       4. Spectral Coherence  
       5. [Principal Component Analysis](#Principal-Component-Analysis)  
       6. Canonical Correlation Analysis  

* ##[CPCA_varimax.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/CPCA_varimax.ipynb)  
    
    **Principal Component Analysis of Nowcast Velocities**  
    This notebook describes the process of performing complex principal component analysis (C-PCA) on the SalishSeaCast Nowcast velocity record and explores the results.  
      
    This notebook is organized into the following sections:  
      
       1. [Load Nowcast and GEM Record](#Load-Nowcast-and-GEM-Record)  
       2. [Principal Component Analysis](#Principal-Component-Analysis)  
       3. [Varimax Rotation](#Varimax-Rotation)  
       4. [Explore Results](#2016-HRDPS-Wind-EOFs)  

* ##[upwelling_geostrophic_velocities.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/upwelling_geostrophic_velocities.ipynb)  
    
    **Upwelling velocities notebook**  

* ##[timeseries_tools_dev.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/timeseries_tools_dev.ipynb)  
    
    **`timeseries_tools` Development Notebook**  
    This notebook is for developing a prototype Nowcast timeseries analysis package. The primary goal of this package is to be memory-efficient. As such, the results arrays are flattened to 2-D (time, space) so that land points can be removed. The 2-D dimensions are also ideal for some analyses like PCA. The basic workflow proceeds as follows:  
       * Flatten the model grid and mask to 2-D (time, space) and remove land indices  
       * Load, process, and flatten hourly Nowcast Results to 2-D (time, space) and remove land indices  
       * Concatenate consecutive 24 hour periods  
       * Reshape the concatenated `Numpy ndarray` to (time, depth, y, x)  

* ##[GRL2016_data.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/GRL2016_data.ipynb)  
    
* ##[NEMOConfigSetup.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/NEMOConfigSetup.ipynb)  
    
* ##[ONC_API_usage.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/ONC_API_usage.ipynb)  
    
* ##[DensitySections.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/DensitySections.ipynb)  
    
    **Nowcast Upwelling Analysis**  

* ##[Tug_spill_preliminary.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Tug_spill_preliminary.ipynb)  
    
* ##[OS2018_plots.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/OS2018_plots.ipynb)  
    
    **Hindcast timeseries analysis**  

* ##[Maps.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Maps.ipynb)  
    
    **Salish Sea Maps**  
    This notebook plots various maps relevant to the SalishSeaCast project. The first map is a general study area map. The second map shows the model domain and bathymetry. The third map shows the Strait of Georgia with relevant observation platform locations overplotted.  

* ##[ShallowWaterModelSummary.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/ShallowWaterModelSummary.ipynb)  
    
    **Shallow water models of coastal upwelling**  
      
    This notebook summarizes several attempts to model coastal upwelling in the Strait of Georgia using the shallow water equations. These modelling attempts are an effort to explain the dynamics of enhanced upwelling around capes that we observe in the SalishSeaCast model (see [**Motivation**](#Motivation)). The following cases are explored:  
      
       1. [**Linear time dependence**](#Linear-time-dependence)  
           * [Surface height sink near coast](#Surface-height-sink-near-coast)  
           * [Linear bottom friction](#Linear-bottom-friction)  
           * [Longshore windstress](#Longshore-windstress)  
       2. [**Steady state**](#Steady-state)  
           * [Linear bottom friction, non-zero depth gradient](#Linear-bottom-friction,-non-zero-depth-gradient)  
           * [Longshore bottom friction, cross-shelf bottom slope](#Longshore-bottom-friction,-cross-shelf-bottom-slope)  

* ##[SedimentModel.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/SedimentModel.ipynb)  
    

##License

These notebooks and files are copyright 2013-2018
by the Salish Sea MEOPAR Project Contributors
and The University of British Columbia.

They are licensed under the Apache License, Version 2.0.
http://www.apache.org/licenses/LICENSE-2.0
Please see the LICENSE file for details of the license.
