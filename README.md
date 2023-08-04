# scRepo
A simple to use webtool for visualising your single cell data.

Given a data table scRepo will generate a webserver using dash.

The only requirement is that the data should be in h5ad format and it should have at least only dimensionality reduction plots in adata.obsm.

# Why use this?
Two important aspect of scRepo is 1) ease of deployment and 2) ease of adding new data to the dataset.
With least investment of time you can customise the deployment to your needs.

* An ideal scenario is to use the tool to maintain all the single cell data available in a group.

* Another scenario is for running a remote job for visualising all the data that you or your collaborators have in the HPC.


# Installation 

`sudo singularity build screpo.sif Singularity.v1.txt`

or if you do not have the privileges to generate a sif file, then use requirements.txt and pip to generate a virtual enviroment.

`pip -r requirements.txt`

# Deployment

## Prepare the data table.
The data table should be in csv format, each row for a h5ad file.

Only column that is mandatory is **"File"** which should correspond to where to find the h5ad file. The path can be relative or absolute and should be accessible to the user running the server.

An example datatable is provided - `data.csv`.

## Basic running of webserver.

```
cd $HOME
singualarity exec /path/to/screpo.sif screpo local data.csv 8221
```
Visit 127.0.0.1:8221 in our remote browser to access the webserver
