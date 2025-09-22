# Building Dash app for brain single cell vascular project
This code is for building a Dash app to visualize gene expression data from a single-cell dataset.
It allows users to select a gene and view its average expression across different cell types in a bar plot.
The app uses Plotly for visualization and Dash for the web interface.  

# Prerequest
1.Corresponding singel cell object & UCSC single cell brwoser (https://cellbrowser.readthedocs.io/en/master/installation.html)
2. Conda environment:
  import scanpy as sc
  import pandas as pd
  import numpy as np
  import plotly.express as px
  from dash import Dash, dcc, html, Input, Output
  import subprocess
  import atexit
  import webbrowser
  import socket
