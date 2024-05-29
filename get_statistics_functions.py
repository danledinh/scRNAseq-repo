"""
This script stores all functions needed to perform the get_statistics.py script

It includes all functions required for loading (h5 files), processing, and summarizing the results
of a scRNA-seq analysis in the context of the scRNA Multimethod comparison project

These functions will be imported into get_statistics.py
"""

##Import necessary packages

import numpy as np
import scanpy as sc
import sys
import os
import functools
import warnings
from anndata import AnnData

#set environment variables so that we can connect to R
# Get the current PATH
current_path = os.getenv('PATH')

# Define the path to the R binary
conda_path = sys.prefix ##sys.prefix actually gives python path, but we are using anaconda for python, so this should work

# Add the R bin path to the current PATH and adding env variables for R
os.environ['PATH'] = current_path

##import custom functions

##define function location
# Get the directory of this function file (located within whole function directory)
functions_dir = os.path.dirname(os.path.abspath(__file__))

##if the directory is not in path, add
if functions_dir not in sys.path:
    sys.path.append(functions_dir)
    print(f"Appended '{functions_dir}' to sys.path.")
else:
    print(f"'{functions_dir}' already exists in sys.path.")
  
##import functions
import highly_deviant_genes as hdg ##import all relevant highly deviant gene detection functions

##Functions

##define decorator function that will allow us to suppress print calls if verbose=True
def suppress_print(func):
    
    """
    Decorator function to suppress print calls if verbose=False.

    Inputs:
        func (function): The function to be decorated.

    Outputs:
        function: The decorated function.
    """
    
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Extract the 'verbose' argument from the keyword arguments or set it to True by default
        verbose = kwargs.pop('verbose', True)
        
        if not verbose:
            # If verbose is False, redirect the standard output to /dev/null (Linux/macOS) or 'nul' (Windows)
            original_stdout = sys.stdout
            sys.stdout = open('/dev/null', 'w')  # For Linux/macOS
            # sys.stdout = open('nul', 'w')  # For Windows

        try:
            # Call the original function with the provided arguments and keyword arguments
            return func(*args, **kwargs)
        finally:
            if not verbose:
                # If verbose is False, restore the original standard output and close the temporary output stream
                sys.stdout.close()
                sys.stdout = original_stdout

    return wrapper

    
##function to load AnnData
@suppress_print
def load_AnnData(filename, **kwargs):
    
    """
    Inputs:
    
    filename: full path to file to read (h5 or h5ad file)
    **kwargs: will accept any key word arguments that can be passed to scanpy's read_10x_h5 function (not relevant if .h5ad file)
    
    Output:
    Anndata object with gene expression data. The object will also have mitochondrial genes annotated and some summary
    statistics calculated via sc.pp.calculate_qc_metrics()
    
    QC Metrics Include:
        1. log1p(total counts)
        2. log1p(gene counts)
        3. Percentage of top 20 genes in counts
        4. mitochondrial_percentage
        5. log1p(mitochondrial_percentage)
        
    
    An AnnData object (https://anndata.readthedocs.io/en/latest/anndata.AnnData.html) holds many slots for annotations 
    and different representations of the data. It also comes with its own HDF5-based file format: `.h5ad`.
    
    Note that we set the default value of the keyword argument 'gex_only' to True so we only read in expression data
    
    Also note that we assume any .h5 formatted file is a 10X h5 file. We currently do not process non-10X .h5 files
    
    """
    
    ##Check type of input file (.h5: read_10X_h5 and .h5ad)
    if filename.endswith('h5'): ##Assuming 10X file
        
        ##check gex_only in keyword arguments
        gex_only = kwargs.get('gex_only', True)  # Get the value of 'gex_only' or set it to True

        # Modify kwargs if 'gex_only' is not provided (defaults to True if not provided)
        if 'gex_only' not in kwargs:
            kwargs['gex_only'] = gex_only

        ##suppress var names warning since we call var_names_make_unique
        warnings.filterwarnings('ignore', message="Variable names are not unique. To make them unique, call `.var_names_make_unique`.")

        ##read h5 file
        adata = sc.read_10x_h5(filename = filename, ##file path for .h5 file
                               **kwargs) ##uses any additional keyword arguments provided (gex_only = True by default)
        
        ##check if variable names are unique; if not, make them unique
        if (len(adata.var.index) - len(np.unique(adata.var.index))) > 0:

            ##make variables have unique names
            adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

        # Restore default warning behavior (optional)
        warnings.resetwarnings()
        
    if filename.endswith('h5ad'): ##Any h5ad file
        
        ##suppress var names warning since we call var_names_make_unique
        warnings.filterwarnings('ignore', message="Variable names are not unique. To make them unique, call `.var_names_make_unique`.")

        ##read in h5ad file
        adata = sc.read_h5ad(filename)
        
        ##check if variable names are unique; if not, make them unique
        if (len(adata.var.index) - len(np.unique(adata.var.index))) > 0:

            ##make variables have unique names
            adata.var_names_make_unique() 

        # Restore default warning behavior (optional)
        warnings.resetwarnings()
        
    
    ##add a stats attribute to adata for storing summary statistics
    adata.uns['stats'] = dict()
    
    return adata

##define function to read in raw data
def load_raw_AnnData(filename):
    
    """
    Read in raw data (ONLY h5 file for 10X data provides this) as an AnnData object.

    Input:
    filename (str): Filename of the raw data (h5 file)

    Output:
    AnnData: Raw AnnData object.
    """
    
    ##suppress var names warning since we call var_names_make_unique
    warnings.filterwarnings('ignore', message="Variable names are not unique. To make them unique, call `.var_names_make_unique`.")
    
    adata_raw = sc.read_10x_h5(filename=filename)
    adata_raw.var_names_make_unique()
    
    # Restore default warning behavior (optional)
    warnings.resetwarnings()
    
    return adata_raw

##define function to normalize data
@suppress_print
def log1p_normalize(adata):
    
    """
    Input: adata object with counts data stored in adata.X (should be unnormalized)
    
    Output: adata object with log1p(size corrected counts) stored in x and unnormalized counts 
            stored in adata.layers['unnormalized']
            
            
    Theory:
    
    The counts we have in our data are generated after cell capture, reverse transcription, amplification, and sequencing. 
    These steps inherently vary per cell, so the counts we see represent not only the biological variation per cell, but
    the technical variation as well. 
    
    Normalizing data is a valuable preprocessing step to adjust these counts in the dataset for technical variance by 
    scaling the observed variance into a specific range. There are many techniques (log shifted transformation, pearson r
    esiduals, etc.) used to make subsequent analysis and statistics applicable, but their usage depends on the situation at 
    hand.
    
    According to the Theis et. al's Single Cell Best Practices book, the shifted logarithm approach is useful for stabilizing
    variance and identifying differentially expressed genes while the pearson residual approach is useful for identifying biologically 
    relevant genes and rare cell types. 
    
    A recent benchmark by [Ahlmann-Eltze & Huber (2023)](https://doi.org/10.1038/s41592-023-01814-1)  
    revealed that the shifted logarithm approach,  demonstrates superior performance compared to other methods in
    uncovering underlying latent structures and stabilizing variance for identifying differentially expressed genes, 
    especially when followed by PCA. In this approach, the log-shifted counts are defined as:
    
   We will use the log-shifted normalization method because it performs well in the benchmark study and works to 
   identify differentially expressed genes. However, we will use the default size factor scaling
   provided in scanpy's `sc.pp.normalize_total()` function instead of the average counts as in the Ahlmann-Eltze & Huber
   paper. scanpy's default size factor scale is the median count across all cells, which should produce similar results 
   to the Ahlmann-Eltze & Huber paper (the only difference is that we use median instead of mean counts). This is performed
   below. 
   
   Note that we will not set a target_sum in `sc.pp.normalize_total` because we don't want to define a set scaling factor 
   for the data (ie $L=10^6$ would give us counts per million).
   """
    
    ## Normalizing Data

    ##recalculate metrics to get total_counts before normalizing
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], log1p = True, percent_top=[20], inplace=True)

    ##check if data has already been normalized (will have an 'unnormalized' layer)
    if adata.layers.get('unnormalized') is None:

        warnings.warn("Normalizing (and log1p) total counts from adata.X. Make sure adata.X is not already normalized, or it may cause problems.")

        #save unnormalized count data to adata.layers['unnormalized'] for feature selection using deviance
        adata.layers['unnormalized'] = adata.X.copy()

    else: ##ie it has already been normalized with log1p_normalize

        warnings.warn(f"Data already normalized. Using adata.layers['unnormalized'] to normalize and log1p data. Storing at adata.X")
        adata.X = adata.layers['unnormalized'].copy()

    ##Either way, scale and log-shift transform data
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)##does the y/s_c calculation
    adata.X = sc.pp.log1p(scales_counts["X"], copy=True) ##does log-shift ie log(y_scaled + 1)
    
    return adata