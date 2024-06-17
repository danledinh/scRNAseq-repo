##The following functions serve to identify highly deviant genes in scRNA data
##This method was suggested by Theis in the single-cell bast practices book (https://www.sc-best-practices.org/preprocessing_visualization/feature_selection.html)
##The specific method was identified and tested in Townes et. al (https://doi.org/10.1186/s13059-019-1861-6)
##The theory is based off the source code from the R scry package (https://rdrr.io/bioc/scry/src/R/featureSelection.R)
##This is the python implementation of this code. It has been tested against the scry results for synthetic matrices
##The sparse and dense implementations have also been tested against each other and with real 10X data

import numpy as np
from scipy.sparse import csr_matrix
# from scipy.sparse import coo_matrix
from scipy.sparse import issparse
from sklearn.impute import SimpleImputer

def validate_inputs(X, sz):
    '''
    Validate the inputs for the deviance computation functions.
    X should be a scipy sparse matrix or a numpy array.
    sz should be a numpy array or a single numeric value.
    '''
    
    # Validate X
    if not (issparse(X) or isinstance(X, np.ndarray)):
        raise ValueError("X should be a scipy sparse matrix or a numpy array")
    
    # Validate sz
    if not (isinstance(sz, np.ndarray) or np.isscalar(sz)):
        raise ValueError("sz should be a numpy array or a single numeric value")


def sparseBinomialDeviance(X, sz):
    '''
    This function calculates binomial deviance per gene in a cell x gene matrix of single-cell RNA-Seq data, 
    which can be used for feature selection.

    Parameters:
    X: A scipy sparse matrix (features in columns and cells in rows) representing gene counts.
    sz: A scale factor, representing total reads per cell (individual cell library size).

    Returns:
    Deviance: A measure of the goodness of fit of a binomial model to the data. 

    Theory:
    In general, deviance serves to calculate the difference between two models' fit to the data.
    In this case, it is defined as the difference between log-likelihoods of a saturated model 
    (a perfect model that fits the data perfectly) and a null model 
    (a model that assumes constant gene counts across all cells for a given gene).
    
    The deviance is defined as D = 2 * (log_lik(saturated) - log_like(null)).

    Assumption:
    In this function, we assume the count data follows a binomial distribution. In other words,
    it assumes the counts for each gene arise from a series of Bernoulli trials, 
    where each trial represents the chance of a particular read being assigned to a particular gene.
    This type of model is a valid approach when attempting feature selection: 
    See (Townes et. al: https://doi.org/10.1186/s13059-019-1861-6)
    
    Application:
    In this case, a high deviance indicates a case where the null model poorly describes the data 
    (i.e., the gene is highly variable across cells). Conversely, a low deviance indicates a case where the null
    model describes the data well (i.e., gene counts are not highly variable across cells).

    Notes:
    The log-likelihood for a Binomial model (for one cell) is:

    l(p, X) = X*log(p) + (n-X)*log(1-p)

    where p is the probability of a success (assigning a read to a gene for a cell), 
    X is the number of successes (assigned reads for a given cell), and n is the total number
    of trials (total number of reads in a cell across all genes).

    To calculate the total log-likelihood for a gene, you sum up the individual log-likelihoods for each cell.
    
    Reference:
    Townes, F. W. et. al (2019). Feature selection and dimension reduction for single-cell RNA-Seq 
    based on a multinomial model. Genome Biology, 20, 295.
    https://doi.org/10.1186/s13059-019-1861-6
    '''
    # Validate Inputs
    validate_inputs(X, sz) 

    # Convert X to CSR format if not already
    X = csr_matrix(X)  

    # Check to make sure X has eliminated zeros, for some reason it doesn't always
    X.eliminate_zeros()

    ##Calculate stats for saturated likelihood calculations
    logp =  X.multiply(1/sz)  # Compute log(p) for each cell/feature combo
    log1mp = logp.copy() # Make a copy for log(1-p) calculation
    logp.data = np.log(logp.data)  # Log transform non-zero probabilities
    log1mp.data = np.log1p(-log1mp.data)  # Log transform (1-p) with log1p(-x) which is more stable

    ##Calculate stats for null likelihood
    sz_sum = np.sum(sz) # Total number of reads for all cells
    feature_sums = np.sum(X, axis=0) # Total reads for each gene across all cells
    p_null = feature_sums / sz_sum  # Average proportion of reads for a given gene across all cells
    log1mp_null = np.log1p(-p_null)  # Log transform (1-p) with log1p(-x) which is more stable

    ##Calculate likelihoods
    ##saturated likelihood uses true data for X and p
    ##Null likelihood uses constant gene expression 
    ##and average gene proportion across all cells for X and p 
    ll_sat = X.multiply(logp - log1mp) + log1mp.multiply(sz)   # Compute saturated log-likelihood for each gene (element-wise)
    ll_sat = np.sum(ll_sat, axis=0)  # Add up log-likelihoods for each gene across all cells (cols)
    ll_null = np.asarray(feature_sums) * np.asarray(np.log(p_null) - log1mp_null) + sz_sum * log1mp_null  # Compute null log-likelihood (convert to array for braodcasting)

    ##Calculate deviance and return
    return (2 * np.asarray(ll_sat - ll_null)).squeeze()  # Compute and return the deviance as an array (with row dim dropped)


def denseBinomialDeviance(X, sz):
    '''
    Function to compute Binomial Deviance for a dense matrix X.
    X should be a numpy array with features in columns and observations in rows.
    sz is a scale factor.

    This function calculates binomial deviance per gene in a cell x gene matrix.
    Deviance is a generalized calculation of residual sum of squares that finds the difference in
    log-likelihoods between two models. It can be calculated as:
    
    D = 2(log_lik(saturated) - log_like(null/proposed))

    In this function, we assume the count data follows a binomial distribution. In other words,
    it assumes the counts for each gene arise from a series of Bernoulli trials, 
    where each trial represents the chance of a particular read being assigned to a particular gene.
    This type of model is a valid approach when attempting feature selection: 
    See (Townes et. al: https://doi.org/10.1186/s13059-019-1861-6)
    
    In this case, It calculates the deviance between the saturated model (perfect model) and the 
    null model (model that assumes gene counts will be constant across all cells for a given gene).
    As such, a high deviance describes a case where the null model poorly describes the data 
    (ie the gene is highly variable across cells) and a low deviance describes the case where the null
    model describes the data well (ie genes are not highly variable).

    We can then use the deviance values to select features (genes) that best-describe the variability
    in the data by selecting those with high deviance (ie those whose genes are highly variable).

    Note that the log-likelihood for a Binomial model (for one cell) is:

    l(p, X) = X*log(p) + (n-X)*log(1-p)

    where p is the probability of a success (assigning a read to a gene for a cell), 
    X is the number of successes (assigned reads for a given cell), and n is the total number
    of trials (total number of reads in a cell across all genes).

    To calculate the total log-likelihood for a gene, you would sum up the individual
    log-likelihoods for each cell.
    '''
    # Validate Inputs
    validate_inputs(X, sz) 

    ##ingore log error with zeros
    with np.errstate(divide='ignore', invalid='ignore'): 

        ##Calculate stats for saturated likelihood calculations
        logp =  X / sz  # For each cell/feature combo, compute prob of assigning a read to a gene
        log1mp = logp.copy() # Make a copy of prob for log(1-p) calc
        logp = np.log(logp)  # Log transform non-zero probabilities
        log1mp = np.log1p(-log1mp)  # Log transform (1-p) (log1p(-x) = log(1-x) but is more stable)

        ##Calculate stats for null likelihood
        sz_sum = np.sum(sz) # sz is library size per cell; so sz_sum is library for whole sample
        feature_sums = np.sum(X, axis=0) # Column sum to get total gene reads across all cells (ignore nan)
        p_null = feature_sums / sz_sum  # Compute average proportion of reads for a given gene across all cells
        log1mp_null = np.log1p(-p_null)  # Log transform (1-p_null) (log1p(-x) = log(1-x) but is more stable)

        ##Calculate likelihoods
        ##saturated likelihood uses true data for X and p
        ##Null likelihood uses constant gene expression 
        ##and average gene proportion across all cells for X and p 
        ll_sat = np.asarray(X) * np.asarray(logp - log1mp) + np.asarray(sz) * np.asarray(log1mp) # Compute saturated log-likelihood for each gene (element-wise)           
        ll_sat = np.nansum(ll_sat, axis=0)  # Column sum to add up log-likelihoods for each gene across all cells (ignore nan)
        ll_null = np.asarray(feature_sums) * np.asarray((np.log(p_null) - log1mp_null)) + np.asarray(sz_sum) * np.asarray(log1mp_null)  # Compute null log-likelihood (convert to array for braodcasting)

    ##Calculate deviance and return
    return (2 * (ll_sat - ll_null)).squeeze()  # Compute and return the deviance (make sure 1 dimensional)



def sparsePoissonDeviance(X, sz):
    """
    This function calculates Poisson deviance per gene in a cell x gene matrix of single-cell RNA-Seq data, 
    which can be used for feature selection.

    Parameters:
    X: A scipy debse matrix (features in columns and cells in rows) representing gene counts.
    sz: A scale factor, representing total reads per cell (individual cell library size).

    Returns:
    Deviance: A measure of the goodness of fit of a Poisson model to the data.

    Theory:
    In general, deviance serves to calculate the difference between two models' fit to the data.
    In this case, it is defined as the difference between log-likelihoods of a saturated model 
    (a perfect model that fits the data perfectly) and a null model 
    (a model that assumes constant gene expression rate across all cells for a given gene).
    The deviance is defined as D = 2 * (log_lik(saturated) - log_lik(null)).

    Assumption:
    In this function, we assume the count data follows a Poisson distribution. 
    This is a common assumption in RNA-Seq analysis, where each gene's read count is modeled as a Poisson 
    variable with the rate parameter lambda representing the true expression level of the gene. 
    Each gene's expression level can vary across cells, thus each cell could be said to have its own lambda. 

    This type of model is a valid approach when attempting feature selection: 
    See (Townes et. al: https://doi.org/10.1186/s13059-019-1861-6)

    Application:
    In this case, a high deviance indicates a case where the null model poorly describes the data 
    (i.e., the gene expression rate is highly variable across cells). Conversely, a low deviance 
    indicates a case where the null model describes the data well (i.e., gene expression rate is 
    not highly variable across cells).

    Notes:
    The log-likelihood for a Poisson model (for one cell) is:

    l(lambda, X) = X*log(lambda) - lambda - log(X!)

    where lambda is the rate parameter (average gene expression level for a given gene), 
    and X is the observed read count for a given cell.  

    To calculate the total log-likelihood for a gene, you sum up the individual log-likelihoods for each cell.
    In our case, the lambdas will be unique for each cell in the saturated model 
    and constant across cells in the null model (representing no variation between cells for each gene)

    In this function, for computational efficiency the factorial term in the Poisson likelihood is ignored.
    It is common practice in statistics to ignore constant terms unrelated to the parameter of interest (lambda here). 
    As such, this simplification is justified as the relative deviance between genes, used for feature selection, remains unaffected.

    In addition, we also ignore the `- lambda` term as is done in the Townes paper because including it
    it could potentially bias the deviance calculation toward cells with greater sequencing depth
    (This is particularly because sz for the Poisson is scaled such that is has a geometric mean of 1,
    which has a disproportionate affect on the `- lambda` term compared to the log(lambda) term)

    The method here (also used in Townes) should provide a sufficiently accurate representation
    of deviance that will allow us to rank the genes (Higher approximate deviance refering to more diverse expression)

    Reference:
    Townes, F. W. et. al (2019). Feature selection and dimension reduction for single-cell RNA-Seq 
    based on a multinomial model. Genome Biology, 20, 295.
    https://doi.org/10.1186/s13059-019-1861-6
    """

    # Validate Inputs
    validate_inputs(X, sz) 

    # Convert X to CSR format if not already
    X = csr_matrix(X)  

    # Check to make sure X has eliminated zeros, for some reason it doesn't always
    X.eliminate_zeros()

    ##Calculate stats for saturated likelihood calculations
    logp =  X.multiply(1/sz)  # Compute log(p) for each cell/feature combo
    lambda_sat = logp.copy() # Individual lambdas for the saturated model
    logp.data = np.log(logp.data)  # Log transform non-zero probabilities

    ##Calculate stats for null likelihood
    sz_sum = np.sum(sz) # Total number of reads for all cells
    feature_sums = np.sum(X, axis=0) # Total reads for each gene across all cells
    n = X.shape[0] # Total number of cells
    lambda_null = feature_sums / sz_sum  # Null lambda assuming constant gene expression across all cells

    ##Calculate likelihoods

    #saturated likelihood uses true data for X and lambda
    ##Null likelihood uses constant gene expression (lambda) across all cells for each gene 
    ll_sat = X.multiply(logp) # Compute saturated log-likelihood for each gene (element-wise)
    ll_sat = np.sum(ll_sat, axis = 0) # Add up log-likelihoods for each gene across all cells (cols; unique lambda per cell)
    ll_null = np.asarray(feature_sums) * np.asarray(np.log(lambda_null)) # Compute null log-likelihood (convert to array for braodcasting)

    ##Calculate deviance and return
    return (2 * np.asarray(ll_sat - ll_null)).squeeze()  # Compute and return the deviance as an array (with row dim dropped)


def densePoissonDeviance(X, sz):
    """
    This function calculates Poisson deviance per gene in a cell x gene matrix of single-cell RNA-Seq data, 
    which can be used for feature selection.

    Parameters:
    X: A scipy debse matrix (features in columns and cells in rows) representing gene counts.
    sz: A scale factor, representing total reads per cell (individual cell library size).

    Returns:
    Deviance: A measure of the goodness of fit of a Poisson model to the data.

    Theory:
    In general, deviance serves to calculate the difference between two models' fit to the data.
    In this case, it is defined as the difference between log-likelihoods of a saturated model 
    (a perfect model that fits the data perfectly) and a null model 
    (a model that assumes constant gene expression rate across all cells for a given gene).
    The deviance is defined as D = 2 * (log_lik(saturated) - log_lik(null)).

    Assumption:
    In this function, we assume the count data follows a Poisson distribution. 
    This is a common assumption in RNA-Seq analysis, where each gene's read count is modeled as a Poisson 
    variable with the rate parameter lambda representing the true expression level of the gene. 
    Each gene's expression level can vary across cells, thus each cell could be said to have its own lambda. 

    This type of model is a valid approach when attempting feature selection: 
    See (Townes et. al: https://doi.org/10.1186/s13059-019-1861-6)

    Application:
    In this case, a high deviance indicates a case where the null model poorly describes the data 
    (i.e., the gene expression rate is highly variable across cells). Conversely, a low deviance 
    indicates a case where the null model describes the data well (i.e., gene expression rate is 
    not highly variable across cells).

    Notes:
    The log-likelihood for a Poisson model (for one cell) is:

    l(lambda, X) = X*log(lambda) - lambda - log(X!)

    where lambda is the rate parameter (average gene expression level for a given gene), 
    and X is the observed read count for a given cell.  

    To calculate the total log-likelihood for a gene, you sum up the individual log-likelihoods for each cell.
    In our case, the lambdas will be unique for each cell in the saturated model 
    and constant across cells in the null model (representing no variation between cells for each gene)

    In this function, for computational efficiency the factorial term in the Poisson likelihood is ignored.
    It is common practice in statistics to ignore constant terms unrelated to the parameter of interest (lambda here). 
    As such, this simplification is justified as the relative deviance between genes, used for feature selection, remains unaffected.

    In addition, we also ignore the `- lambda` term as is done in the Townes paper because including it
    it could potentially bias the deviance calculation toward cells with greater sequencing depth
    (This is particularly because sz for the Poisson is scaled such that is has a geometric mean of 1,
    which has a disproportionate affect on the `- lambda` term compared to the log(lambda) term)

    The method here (also used in Townes) should provide a sufficiently accurate representation
    of deviance that will allow us to rank the genes (Higher approximate deviance refering to more diverse expression)

    Reference:
    Townes, F. W. et. al (2019). Feature selection and dimension reduction for single-cell RNA-Seq 
    based on a multinomial model. Genome Biology, 20, 295.
    https://doi.org/10.1186/s13059-019-1861-6
    """

    # Validate Inputs
    validate_inputs(X, sz) 

    ##ingore log error with zeros
    with np.errstate(divide='ignore', invalid='ignore'): 

        ##Calculate stats for saturated likelihood calculations
        logp =  X / sz  # For each cell/feature combo, compute prob of assigning a read to a gene
        lambda_sat = logp.copy() # Individual lambdas for the saturated model
        logp = np.log(logp)  # Log transform non-zero probabilities

        ##Calculate stats for null likelihood
        sz_sum = np.sum(sz) # Total number of reads for all cells
        feature_sums = np.sum(X, axis=0) # Total reads for each gene across all cells
        n = X.shape[0] # Total number of cells
        lambda_null = feature_sums / sz_sum  # Null lambda assuming constant gene expression across all cells

        ##Calculate likelihoods

        #saturated likelihood uses true data for X and lambda
        ##Null likelihood uses constant gene expression (lambda) across all cells for each gene 
        ll_sat = np.asarray(X)*np.asarray(logp) # Compute saturated log-likelihood for each gene (element-wise)
        ll_sat = np.nansum(ll_sat, axis = 0) # Add up log-likelihoods for each gene across all cells (cols; unique lambda per cell)
        ll_null = np.asarray(feature_sums) * np.asarray(np.log(lambda_null)) # Compute null log-likelihood 


    ##Calculate deviance and return
    return (2 * (ll_sat - ll_null)).squeeze()  # Compute and return the deviance as an array (make sure 1 dimensional)


def check_replace_nan(matrix):
    ##for both sparse and dense matrices, check if there are nans
    ##if there are, output a warning for how many there are and that we are replacing with "0"
    ##Then replace with "0"

    if issparse(matrix):
        if np.isnan(matrix.data).any():
            nan_count = np.isnan(matrix.data).sum()
            imp = SimpleImputer(missing_values=np.nan, strategy='constant', fill_value=0)
            matrix = imp.fit_transform(matrix)
            matrix.eliminate_zeros() ##remove zeros created by imputation to retain sparsity
            print(f" \n WARNING: Replaced {nan_count} NaN value(s) with 0 in sparse matrix. \n")
    else:
        if np.isnan(matrix).any():
            nan_count = np.isnan(matrix).sum()
            imp = SimpleImputer(missing_values=np.nan, strategy='constant', fill_value=0)
            matrix = imp.fit_transform(matrix)
            print(f" \n WARNING: Replaced {nan_count} NaN value(s) with 0 in dense matrix. \n")
            
    return matrix

def compute_deviance(m, fam='binomial'):
    '''
    m is either a dense numpy array or a scipy sparse matrix
    m is a data matrix with genes as columns and cells as rows
    '''

    # validate family parameter
    fam_options = ['binomial', 'poisson']
    if fam not in fam_options:
        raise ValueError(f"fam must be one of {fam_options}")

    ##For m, replace nan with 0 values. This is not ideal for traditional imputation,
    ##but considering count data is mostly zeros, filling in np.nan with 0 is a valid approach
    ##and should not change deviance much, especiallly because we don't expect many nans
    ##however, we will supply a warning if there are any nans and the strategy used to fill them
    m = check_replace_nan(m)

    # compute total gene counts per cell (row sums so axis = 1; ie add all columns)
    sz = np.sum(m, axis=1) ##nans won't be counted due to imputation above
    sz = sz.reshape((-1,1)) ##make sure it has a second dimension so it can broadcast

    if fam == "poisson":
        # In the case of a Poisson model, we want to scale the size factors (sz) so that they 
        # have a geometric mean of 1. This is done to ensure similar sequencing depth across cells,
        # as variations in sz can cause variations in counts not due to biological variations.
        lsz = np.log(sz)
        sz = np.exp(lsz - np.mean(lsz))

    # choose computation method based on matrix sparsity and family
    if issparse(m):
        if fam == "binomial":
            return sparseBinomialDeviance(m, sz)
        else:  # fam == "poisson"
            return sparsePoissonDeviance(m, sz)
    else:
        #m is either 1) an ordinary dense array or matrix
        # 2) a non-sparse Matrix
        # 3) a dense object like HDF5Array (on disk)
        if fam == "binomial":
            return denseBinomialDeviance(m, sz)
        else:  # fam == "poisson"
            return densePoissonDeviance(m, sz)


def highly_deviant_genes(adata=None, X=None, top_pct=None, top_n=None, family='binomial', gene_names=None):
    '''
    Input:

    adata: anndata.AnnData object. The data matrix is assumed to be stored in adata.X if adata is provided
    X: np.array or None, optional data matrix (cols by genes format).
       If specified, this matrix is used for calculation and the output is a 
       dictionary of deviance scores, deviance masks, and gene names (if provided)
    top_pct: float, top percentage of genes to identify (0.0 to 1.0); default is 0.15
    top_n: integer or None, number of top most deviant genes to identify; default is None
    family: string, either 'binomial' or 'poisson'. Specifies the method to calculate deviance
    gene_names: list or None, optional list of gene names corresponding to the columns of X 
                (only to be used) if X is specified

    Output:
    If adata provided, adata ouput with deviance scores, deviance mask, and deviance calculation method
    If X matrix provided, matrix with deviance scores, deviance mask, (and gene names if provided)
    '''

    # Validate family parameter
    fam_options = ['binomial', 'poisson']
    if family not in fam_options:
        raise ValueError(f"family must be one of {fam_options}")

    # Check if both adata and X are provided
    if adata is not None and X is not None:
        raise ValueError("Only of of adata or X should be provided, not both.")

    # If X is not specified, use adata.X as the data matrix (if adata has an X attribute)
    if X is None and hasattr(adata, 'X'):
        X = adata.X
    elif X is None and not hasattr(adata, 'X'):
        raise ValueError("\nadata object does not have a default adata.X attribute. \nProvide an adata object with an X attribute (adata.X) or use argument 'X' to provide a matrix.")

    # Check if both top_n and top_pct are specified
    if top_n is not None and top_pct is not None:
        raise ValueError("Only one of top_n or top_pct should be specified, not both.")

    # Calculate deviance
    deviance = compute_deviance(X, family)

    # Determine top genes based on top_n or top_pct
    if top_n is not None:
        top_n = int(top_n)
        if top_n <= 0 or top_n > X.shape[1]:
            raise ValueError("top_n must be a positive integer within the range of available genes")
    else:
        if top_pct is not None:
            top_pct = float(top_pct)
        else:
            top_pct = 0.15 ##default is 15%
        if top_pct <= 0.0 or top_pct > 1.0:
            raise ValueError("top_pct must be a float between 0.0 and 1.0")
        top_n = int(X.shape[1] * top_pct)
    
    idx = deviance.argsort()[-top_n:] # index of top_n

    # Create boolean mask for top genes
    mask = np.zeros(X.shape[1], dtype=bool)
    mask[idx] = True # set those that are "top" as true

    # If input was adata, add results back, if not output a dictionary
    if adata is not None:
        # Add computed results back to adata object
        adata.var["highly_deviant"] = mask
        adata.var["deviance"] = deviance
        adata.var["deviance_calculation_method"] = family

        # return adata object
        return adata

    else:
        # Add computed results to a dictionary
        results = {"deviance": deviance, "highly_deviant": mask}

        # Add gene names if provided
        if gene_names is not None:
            if len(gene_names) != X.shape[1]:
                raise ValueError("The number of gene names provided does not match the number of genes in the data matrix.")
            results["gene_names"] = gene_names

        return results
