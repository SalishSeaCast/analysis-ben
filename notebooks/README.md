The Jupyter Notebooks in this directory are made by Ben for
quick sharing of results.

The links below are to static renderings of the notebooks via
[nbviewer.ipython.org](http://nbviewer.ipython.org/).
Descriptions below the links are from the first cell of the notebooks
(if that cell contains Markdown or raw text).

* ##[Wind_comparisons_jan2016.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Wind_comparisons_jan2016.ipynb)  
    
    **HRDPS Evaluation January 2016**  
    This notebook compares the GEM 2.5 km HRDPS wind product to local windstation records for January 2016  

* ##[Hakai_data.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Hakai_data.ipynb)  
    
* ##[visualization_workflows_xarray.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/visualization_workflows_xarray.ipynb)  
    
    **Visualisation of NEMO/GEM/Observations using XArray**  
    This notebook demonstrates the use of several tools for easy loading and visualization of model results and drifter observations. Most of these tools require `xarray`, and rather than make them flexible to other packages, I thought I would just demonstrate their use with `xarray` here.  
      
    First we'll import the necessary libraries and set our preferred notebook formatting.  

* ##[HRDPS_correction.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/HRDPS_correction.ipynb)  
    
    **HRDPS Correction Check**  
    This notebook checks the HRDPS correction obtained from Luc Fillion's group  

* ##[VelocityComparisons.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/VelocityComparisons.ipynb)  
    
* ##[DICTA_world_rivers.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/DICTA_world_rivers.ipynb)  
    
* ##[timeseries_tools_usage.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/timeseries_tools_usage.ipynb)  
    
    **Timeseries Tools Usage**  
    This notebook demonstrates loading in a Nowcast timeseries using `timeseries_tools`. `Timeseries_tools` extracts only the required fields from the Nowcast netCDF files and concatenates the flattened (time, space) numpy arrays together. This uses a minimum amount of memory relative to bulkier routines such as `mfdataset`. The flattened arrays are reshaped before they are returned, but there is an option to keep them flattened, which is necessary for certain data analysis methods like principal component analysis.  

* ##[AnalyticalModel.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/AnalyticalModel.ipynb)  
    
    **Analytical Upwelling Model**  

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

* ##[SOG_chem_data.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/SOG_chem_data.ipynb)  
    
* ##[ForceBalanceModelSummary.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/ForceBalanceModelSummary.ipynb)  
    
* ##[CPCA_varimax.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/CPCA_varimax.ipynb)  
    
    **Principal Component Analysis of Nowcast Velocities**  
    This notebook describes the process of performing complex principal component analysis (C-PCA) on the SalishSeaCast Nowcast velocity record and explores the results.  
      
    This notebook is organized into the following sections:  
      
       1. [Load Nowcast and GEM Record](#Load-Nowcast-and-GEM-Record)  
       2. [Principal Component Analysis](#Principal-Component-Analysis)  
       3. [Varimax Rotation](#Varimax-Rotation)  
       4. [Explore Results](#2016-HRDPS-Wind-EOFs)  

* ##[Wind_comparisons_mar2016.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Wind_comparisons_mar2016.ipynb)  
    
    **HRDPS Evaluation March 2016**  
    This notebook compares the GEM 2.5 km HRDPS wind product to local windstation records for March 2016  

* ##[ATW_model.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/ATW_model.ipynb)  
    
* ##[upwelling_geostrophic_velocities.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/upwelling_geostrophic_velocities.ipynb)  
    
    **Upwelling velocities notebook**  

* ##[timeseries_tools_dev.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/timeseries_tools_dev.ipynb)  
    
    **`timeseries_tools` Development Notebook**  
    This notebook is for developing a prototype Nowcast timeseries analysis package. The primary goal of this package is to be memory-efficient. As such, the results arrays are flattened to 2-D (time, space) so that land points can be removed. The 2-D dimensions are also ideal for some analyses like PCA. The basic workflow proceeds as follows:  
       * Flatten the model grid and mask to 2-D (time, space) and remove land indices  
       * Load, process, and flatten hourly Nowcast Results to 2-D (time, space) and remove land indices  
       * Concatenate consecutive 24 hour periods  
       * Reshape the concatenated `Numpy ndarray` to (time, depth, y, x)  

* ##[SpectralAnalysis.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/SpectralAnalysis.ipynb)  
    
* ##[TA.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/TA.ipynb)  
    
    **Fraser/Strait of Georgia Freshwater Chemistry Analysis**  
    This is the data notebook  
    **Sections**  
       1. [Local Functions](#Local-Functions)  
       2. [Load Data](#Load-Data)  
       3. [TA regressions](#TA-regressions)  
       4. [Fraser River Buoy pH Data](#Fraser-River-Buoy-pH-Data)  
       5. [Freshwater DIC:TA scenarios](#Freshwater-DIC:TA-scenarios)  

* ##[NEMOConfigSetup.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/NEMOConfigSetup.ipynb)  
    
* ##[ONC_API_usage.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/ONC_API_usage.ipynb)  
    
* ##[DensitySections.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/DensitySections.ipynb)  
    
    **Nowcast Upwelling Analysis**  

* ##[KelvinWaves.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/KelvinWaves.ipynb)  
    
* ##[SOG_river_figures.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/SOG_river_figures.ipynb)  
    
    **SOG Freshwater Chemistry Sensitivity Notebook**  
    This is the model sensitivity notebook  
    **Sections**  
       1. [Define Local Functions](#Define-Local-Functions)  
       2. [Load and Process Data](#Load-and-Process-Data)  
       3. [**Make Figures**](#Make-Figures)  
           * [Flow Dependence Figures](#Flow-Dependence-Figures)  
           * [SOG Timeseries Figures](#SOG-Timeseries-Figures)  
           * [SOG Salinity Average Figures](#SOG-Salinity-Average-Figures)  
           * [Salinity Space Figures](#Salinity-Space-Figures)  

* ##[AnalyticalModel_friction.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/AnalyticalModel_friction.ipynb)  
    
* ##[Tug_spill_preliminary.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Tug_spill_preliminary.ipynb)  
    
* ##[OS2018_plots.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/OS2018_plots.ipynb)  
    
    **Hindcast timeseries analysis**  

* ##[Maps.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Maps.ipynb)  
    
    **Salish Sea Maps**  
    This notebook plots various maps relevant to the SalishSeaCast project. The first map is a general study area map. The second map shows the model domain and bathymetry. The third map shows the Strait of Georgia with relevant observation platform locations overplotted.  

* ##[Upwelling_analysis.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Upwelling_analysis.ipynb)  
    
* ##[ShallowWaterModelSummary.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/ShallowWaterModelSummary.ipynb)  
    
* ##[WindAnalysis.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/WindAnalysis.ipynb)  
    
    **Wind Analysis**  

* ##[SedimentModel.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/SedimentModel.ipynb)  
    

##License

These notebooks and files are copyright 2013-2018
by the Salish Sea MEOPAR Project Contributors
and The University of British Columbia.

They are licensed under the Apache License, Version 2.0.
http://www.apache.org/licenses/LICENSE-2.0
Please see the LICENSE file for details of the license.
