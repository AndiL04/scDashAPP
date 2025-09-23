# Building Dash app for brain single cell vascular project
This code is for building a Dash app to visualize gene expression data from a single-cell dataset.
It allows users to select a gene and view its average expression across different cell types in a bar plot.
The app uses Plotly for visualization and Dash for the web interface.  

# Data requied
Toy data: "/path/to/the/scDashAPP/merged_normalized_sub.h5ad"

## Install prerequisites
        conda env create -f environment.yml

## initlialize the conda environment
        conda activate dash_app 

## nevigate to the DashAPP directory
        cd /path/to/the/scDashAPP

## bilding the cell browser based on the anndata object
## Based on the tutorial: https://cellbrowser.readthedocs.io/en/master/scanpy.html
## This step only needs to be done once
        cbImportScanpy -i ./merged_normalized_sub.h5ad -o vas_atlas -n vas_atlas

## after building the cell browser, need to modify the cellbrowser.conf file to point to the correct color files
        name='vas_atlas'
        shortLabel='vas_atlas'
        exprMatrix='matrix.mtx.gz'
        #tags = ["10x", 'smartseq2']
        meta='meta.tsv'
        geneIdType='auto'
        defColorField='cell.class'
        labelField='cell.class'
        enumFields=['cell.class']
        coords=[
        {
                "file": "umap.harmony_coords.tsv",
                "shortLabel": "umap.harmony"
        }
        ]
        alpha=0.8
        #body_parts=["embryo", "heart", "brain"]
        radius=2
        markers = [{"file": "markers.tsv", "shortLabel":"Cluster Markers"}]
        quickGenesFile="quickGenes.tsv"
        colors={
        "region_abb": "region_abb_colors.tsv",
        "cell.class": "cell.class_colors.tsv",
        "cell.type": "cell.type_colors.tsv",
        "cell.type.middle": "cell.type.middle_colors.tsv"
        }

## start the cell browser server
        cd ./vas_atlas

## Start the Dash APP
        python Dash_app.py

## Open the add at http://127.0.0.1:8050/


