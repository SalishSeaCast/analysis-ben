## OpenDrift

The Jupyter Notebooks in this directory are made by Ben
Moore-Maley for quick sharing of results.

The links below are to static renderings of the notebooks via
[nbviewer.ipython.org](http://nbviewer.ipython.org/).
Descriptions below the links are from the first cell of the notebooks
(if that cell contains Markdown or raw text).

***
* ### [DrifterSimulations_parcels.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/OpenDrift/DrifterSimulations_parcels.ipynb)  
    
    **Drifter Simulations**  
      
    ***  

***
* ### [OceanParcels_workflow.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/OpenDrift/OceanParcels_workflow.ipynb)  
    
    **OceanParcels workflow**  
      
    Useful resources:  
      
    [OceanParcels Homepage](http://oceanparcels.org/)  
      
    [GMD discussion paper](https://doi.org/10.5194/gmd-2018-339)  
      
    [Coastline interaction scripts](https://github.com/OceanParcels/Parcelsv2.0PaperNorthSeaScripts) (from GMD paper)  
      
    Source code modification (to handle WW3 longitude):  
      
    `~/anaconda3/envs/py3_parcels/lib/python3.7/site-packages/parcels/field.py`  
      
    `@ field.py: line 1385`  
      
    `if lon.max().values > 180: lon = lon - 360`  
      
    ***  

***
* ### [drifter_comparisons.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/OpenDrift/drifter_comparisons.ipynb)  
    
    **Drifter comparisons**  

***
* ### [drifter_comparisons2.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/OpenDrift/drifter_comparisons2.ipynb)  
    
    **Drifter evaluation example**  

***
* ### [drifter_evaluation_master.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/OpenDrift/drifter_evaluation_master.ipynb)  
    
    **Master Drifter Evaluation Notebook**  
      
    ***  

***
* ### [opendrift_debug.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/OpenDrift/opendrift_debug.ipynb)  
    
***
* ### [opendrift_forcing.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/OpenDrift/opendrift_forcing.ipynb)  
    
    **OpenDrift**  
      
    **Documentation**  
      
    [https://github.com/opendrift/opendrift/wiki](https://github.com/opendrift/opendrift/wiki)  
      
    ***  
      
    **Installation**  
      
    `git clone https://github.com/OpenDrift/opendrift.git`  
      
    `python setup.py develop --user`  
      
    ***  

***
* ### [sample_opendrift_simulation.ipynb](http://nbviewer.ipython.org/urls/github.com/SalishSeaCast/analysis-ben/blob/master/notebooks/OpenDrift/sample_opendrift_simulation.ipynb)  
    
    **OpenDrift simulations forced by a regional ocean model (NEMO)**  
      
    ***  


## License

These notebooks and files are copyright 2013-2020
by the Salish Sea MEOPAR Project Contributors
and The University of British Columbia.

They are licensed under the Apache License, Version 2.0.
http://www.apache.org/licenses/LICENSE-2.0
Please see the LICENSE file for details of the license.
