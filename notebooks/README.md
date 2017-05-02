The Jupyter Notebooks in this directory are made by Ben for
quick sharing of results.

The links below are to static renderings of the notebooks via
[nbviewer.ipython.org](http://nbviewer.ipython.org/).
Descriptions below the links are from the first cell of the notebooks
(if that cell contains Markdown or raw text).

* ##[Wind_comparisons_jan2016.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Wind_comparisons_jan2016.ipynb)  
    
    **HRDPS Evaluation January 2016**  
    This notebook compares the GEM 2.5 km HRDPS wind product to local windstation records for January 2016  

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

* ##[SeicheScaling.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/SeicheScaling.ipynb)  
    
    **Scaling Upwelling and Seicheing in the Strait of Georgia**  

* ##[Wind_comparisons_mar2016.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Wind_comparisons_mar2016.ipynb)  
    
    **HRDPS Evaluation March 2016**  
    This notebook compares the GEM 2.5 km HRDPS wind product to local windstation records for March 2016  

* ##[NowcastTimeseries.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/NowcastTimeseries.ipynb)  
    
    **SalishSeaCast Nowcast Timeseries Analysis**  
    This notebook loads hourly Nowcast results 1 day at a time, and extracts selected slices into numpy arrays. The notebook itself is a development environment for the script `analysis-ben/scripts/prepare_nowcast_timeseries.py`.  

* ##[ONC_API_usage.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/ONC_API_usage.ipynb)  
    
* ##[Maps.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/Maps.ipynb)  
    
    **Salish Sea Maps**  
    This notebook plots various maps relevant to the SalishSeaCast project. The first map is a general study area map. The second map shows the model domain and bathymetry. The third map shows the Strait of Georgia with relevant observation platform locations overplotted.  

* ##[WindAnalysis.ipynb](http://nbviewer.ipython.org/urls/bitbucket.org/salishsea/analysis-ben/raw/tip/notebooks/WindAnalysis.ipynb)  
    
    **Wind Analysis**  


##License

These notebooks and files are copyright 2013-2017
by the Salish Sea MEOPAR Project Contributors
and The University of British Columbia.

They are licensed under the Apache License, Version 2.0.
http://www.apache.org/licenses/LICENSE-2.0
Please see the LICENSE file for details of the license.
