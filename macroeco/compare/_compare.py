from __future__ import division

import numpy as np
import scipy as sp
import scipy.stats as stats
import pandas as pd

from ..misc import doc_sub

_data_doc = \
    """data : iterable
        Data for analysis"""

_model_doc = \
    """model : obj
        Scipy frozen distribution object. When freezing, keyword args ``loc``
        and ``scale`` should only be included if they represent a parameter.
    """

_obs_pred_doc = \
    """obs, pred : array-like objects
        Observed and predicted data
    """


@doc_sub(_data_doc, _model_doc)
def nll(data, model):
    """
    Negative log likelihood given data and a model

    Parameters
    ----------
    {0}
    {1}

    Returns
    -------
    float
        Negative log likelihood

    Examples
    ---------

    >>> import macroeco.models as md
    >>> import macroeco.compare as comp

    >>> # Generate random data
    >>> rand_samp = md.logser.rvs(p=0.9, size=100)

    >>> # Get nll for p = 0.9
    >>> comp.nll(rand_samp, md.logser(p=0.9))
    237.6871819262054

    >>> # Get the nll for the MLE for p
    >>> mle_p = md.logser.fit_mle(rand_samp)
    >>> comp.nll(rand_samp, md.logser(*mle_p))
    235.2841347820297

    """

    try:
        log_lik_vals = model.logpmf(data)
    except:
        log_lik_vals = model.logpdf(data)
    return -np.sum(log_lik_vals)


@doc_sub(_data_doc)
def lrt(data, model_full, model_reduced, df=None):
    """
    Compare two nested models using a likelihood ratio test

    Parameters
    ----------
    {0}
    model_full : obj
        A frozen scipy distribution object representing the full model
        (more complex model).
    model_reduced : scipy distribution object
        A frozen scipy distribution object representing the reduced model
        (simpler model).
    df : int
        The degrees of freedom for the lrt (optional). If none, df is
        calculated as the difference between the number of parameters in the
        full and reduced models.

    Returns
    -------
    tuple
        G^2 statistic, p-value

    Notes
    -----
    Parameters of distribution objects must be given as keyword arguments. Ex.
    ``norm = stats.norm(loc=0, scale=1)``

    A p-value < alpha suggests significant evidence for the full (more complex)
    model. In other words, the null hypothesis is that the reduced model is
    correct

    The LRT only applies to nested models. The G^2 statistic and G-test rely on
    the assumption that the test statistic is approximately chi-squared
    distributed. This assumption breaks down for small samples sizes.

    Examples
    --------

    >>> import macroeco.models as md
    >>> import macroeco.compare as comp

    >>> # Generate random data
    >>> rand_samp = md.nbinom_ztrunc.rvs(20, 0.5, size=100)

    >>> # Fit Zero-truncated NBD (Full model)
    >>> mle_nbd = md.nbinom_ztrunc.fit_mle(rand_samp)

    >>> # Fit a logseries (limiting case of Zero-truncated NBD, reduced model)
    >>> mle_logser = md.logser.fit_mle(rand_samp)

    >>> # Compare models with LRT
    >>> comp.lrt(rand_samp, md.nbinom_ztrunc(mu=mle_nbd[0], k_agg=mle_nbd[1]), md.logser(p=mle_logser[0]))
    (15.33429080890221, 9.0066719644695982e-05)

    >>> # Reject the null hypothesis that the logseries is a better model and
    >>> # choose the Zero-truncated NBD.

    """

    # Calculate G^2 statistic
    ll_full = nll(data, model_full) * -1
    ll_reduced = nll(data, model_reduced) * -1
    test_stat = 2 * (ll_full - ll_reduced)

    # Set df if necessary
    if not df:
        df =  ( len(model_full.args) + len(model_full.kwds)
              - len(model_reduced.args) - len(model_reduced.kwds) )

    return test_stat, stats.chisqprob(test_stat, df)


@doc_sub(_data_doc, _model_doc)
def AIC(data, model, params=None, corrected=True):
    """
    Akaike Information Criteria given data and a model

    Parameters
    ----------
    {0}
    {1}
    params : int
        Number of parameters in the model. If None, calculated from model
        object.
    corrected : bool
        If True, calculates the small-sample size correct AICC. Default True.

    Returns
    -------
    float
        AIC(C) value

    Notes
    -----
    AICC should be used when the number of observations is < 40.

    Examples
    --------

    >>> import macroeco.models as md
    >>> import macroeco.compare as comp

    >>> # Generate random data
    >>> rand_samp = md.nbinom_ztrunc.rvs(20, 0.5, size=100)

    >>> # Fit Zero-truncated NBD (Full model)
    >>> mle_nbd = md.nbinom_ztrunc.fit_mle(rand_samp)

    >>> # Fit a logseries (limiting case of Zero-truncated NBD, reduced model)
    >>> mle_logser = md.logser.fit_mle(rand_samp)

    >>> # Get AIC for ztrunc_nbinom
    >>> comp.AIC(rand_samp, md.nbinom_ztrunc(*mle_nbd))
    765.51518598676421

    >>> # Get AIC for logser
    >>> comp.AIC(rand_samp, md.logser(*mle_logser))
    777.05165086534805

    >>> # Support for for zero-truncated NBD over logseries because AIC is
    >>> # smaller

    >>> # Call AIC with params given as 2 (should be the same as above)
    >>> comp.AIC(rand_samp, md.nbinom_ztrunc(*mle_nbd), params=2)
    765.51518598676421

    >>> # Call AIC without sample size correction
    >>> comp.AIC(rand_samp, md.nbinom_ztrunc(*mle_nbd), params=2, corrected=False)
    765.39147464655798

    References
    ----------
    .. [#]
       Burnham, K and Anderson, D. (2002) Model Selection and Multimodel
       Inference: A Practical and Information-Theoretic Approach (p. 66). New
       York City, USA: Springer.

    """
    n = len(data)  # Number of observations
    L = nll(data, model)

    if not params:
        k = len(model.kwds) + len(model.args)
    else:
        k = params

    if corrected:
        aic_value = 2 * k + 2 * L + (2 * k * (k + 1)) / (n - k - 1)
    else:
        aic_value = 2 * k + 2 * L

    return aic_value


def AIC_compare(aic_list):
    """
    Calculates delta AIC and AIC weights from a list of AIC values

    Parameters
    -----------------
    aic_list : iterable
        AIC values from set of candidat models

    Returns
    -------------
    tuple
        First element contains the delta AIC values, second element contains
        the relative AIC weights.

    Notes
    -----
    AIC weights can be interpreted as the probability that a given model is the
    best model in the set.

    Examples
    --------

    >>> # Generate random data
    >>> rand_samp = md.nbinom_ztrunc.rvs(20, 0.5, size=100)

    >>> # Fit Zero-truncated NBD (Full model)
    >>> mle_nbd = md.nbinom_ztrunc.fit_mle(rand_samp)

    >>> # Fit a logseries (limiting case of Zero-truncated NBD, reduced model)
    >>> mle_logser = md.logser.fit_mle(rand_samp)

    >>> # Get AIC for ztrunc_nbinom
    >>> nbd_aic = comp.AIC(rand_samp, md.nbinom_ztrunc(*mle_nbd))

    >>> # Get AIC for logser
    >>> logser_aic = comp.AIC(rand_samp, md.logser(*mle_logser))

    >>> # Make AIC list and get weights
    >>> aic_list = [nbd_aic, logser_aic]
    >>> comp.AIC_compare(aic_list)
    (array([  0.        ,  19.11806518]),
    array([  9.99929444e-01,   7.05560486e-05]))

    >>> # Zero-truncated NBD is a far superior model based on AIC weights

    """

    aic_values = np.array(aic_list)
    minimum = np.min(aic_values)
    delta = aic_values - minimum
    values = np.exp(-delta / 2)
    weights = values / np.sum(values)

    return delta, weights


def sum_of_squares(obs, pred):
    """
    Sum of squares between observed and predicted data

    Parameters
    ----------
    obs : iterable
        Observed data
    pred : iterable
        Predicted data

    Returns
    -------
    float
        Sum of squares

    Notes
    -----
    The length of observed and predicted data must match.

    """

    return np.sum((np.array(obs) - np.array(pred)) ** 2)


def r_squared(obs, pred, one_to_one=False, log_trans=False):
    """
    R^2 value for a regression of observed and predicted data

    Parameters
    ----------
    obs : iterable
        Observed data
    pred : iterable
        Predicted data
    one_to_one : bool
        If True, calculates the R^2 based on the one-to-one line (see [#]_),
        and if False, calculates the standard R^2 based on a linear regression.
        Default False.
    log_trans : bool
        If True, log transforms obs and pred before R^2 calculation.

    Returns
    -------
    float
        R^2 value

    Notes
    -----
    Using the traditional R^2 to compare the fit of observed and predicted
    values may be misleading as the relationship may not be one-to-one but the
    R^2 value may be quite high. The one-to-one option alleviates this problem.
    Note that with the one-to-one option R^2 can be negative.

    Examples
    --------

    >>> import numpy as np
    >>> import macroeco.compare as comp

    >>> # Generate some data
    >>> x_vals = np.linspace(1, 20, num=100)
    >>> y_vals = np.random.normal(4 + x_vals*2, 1)

    >>> # Standard R^2
    >>> comp.r_squared(x_vals, y_vals)
    0.99336568326291697

    >>> # R^2 about the 1:1 line, will be a poor fit (possibly negative)
    >>> comp.r_squared(x_vals, y_vals, one_to_one=True)
    -6.8621799432144988

    >>> # Generate some other data
    >>> y_vals = np.random.normal(x_vals, 1)

    >>> # Normal R^2
    >>> comp.r_squared(x_vals, y_vals)
    0.97651897660174425

    >>> # R^2 on to the one to one line
    >>> comp.r_squared(x_vals, y_vals, one_to_one=True)
    0.97591430200514639

    References
    ----------
    .. [#]
       White, E., Thibault, K., & Xiao, X. (2012). Characterizing the species
       abundance distributions across taxa and ecosystems using a simple
       maximum entropy model. Ecology, 93(8), 1772-8

    """

    if log_trans:
        obs = np.log(obs)
        pred = np.log(pred)

    if one_to_one:
        r_sq = 1 - (sum_of_squares(obs, pred) /
                    sum_of_squares(obs, np.mean(obs)))
    else:
        b0, b1, r, p_value, se = stats.linregress(obs, pred)
        r_sq = r ** 2

    return r_sq

def preston_bin(data, max_num):
    """
    Bins data on base 2 using Preston's method

    Parameters
    ----------
    data : array-like
        Data to be binned
    max_num :  float
        The maximum upper value of the data

    Returns
    -------
    tuple
        (binned_data, bin_edges)

    Notes
    -----
    Uses Preston's method of binning, which has exclusive lower boundaries and
    inclusive upper boundaries. Densities are not split between bins.

    Examples
    --------

    >>> import macroeco.compare as comp
    >>> import numpy as np

    >>> # Load some data and get Preston bins
    >>> data = np.array([1, 1, 1, 1, 4, 5, 6, 7, 12, 34, 56])
    >>> comp.preston_bin(data, np.max(data))
    (array([4, 0, 1, 3, 1, 0, 2]),
    array([  1.,   2.,   3.,   5.,   9.,  17.,  33.,  65.]))

    References
    ----------
    .. [#]
       Preston, F. (1962). The canonical distribution of commonness and rarity.
       Ecology, 43, 185-215

    """

    log_ub = np.ceil(np.log2(max_num))

    # Make an exclusive lower bound in keeping with Preston
    if log_ub == 0:
        boundaries = np.array([0, 1])
    elif log_ub == 1:
        boundaries = np.arange(1, 4)
    else:
        boundaries = 2 ** np.arange(0, log_ub + 1)
        boundaries = np.insert(boundaries, 2, 3)
        boundaries[3:] = boundaries[3:] + 1

    hist_data = np.histogram(data, bins=boundaries)
    return hist_data


def pueyo_bins(data):
    """
    Binning method based on Pueyo (2006)

    Parameters
    ----------
    data : array-like data
        Data to be binned

    Returns
    -------
    : tuple of arrays
        binned data, empirical probability density

    Notes
    -----
    Bins the data in into bins of length 2**i, i=0, 1, 2 ...
    The empirical probability densities will sum to 1 if multiplied by the
    respective 2**i.

    """
    log_ub = np.ceil(np.log2(np.max(data)))
    bins = 2**np.arange(log_ub + 1)
    binned_data = np.histogram(data, bins=bins)[0]
    epdf = (1 / bins[:-1]) * binned_data / len(data)
    return binned_data, epdf


