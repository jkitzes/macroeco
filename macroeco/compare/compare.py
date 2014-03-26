from __future__ import division

import numpy as np
import scipy as sp
import scipy.stats as stats
import pandas as pd

from ..misc import doc_sub

_data_doc = \
    """data : array-like
        data from which to caculate the the likelihood
    """

_model_doc = \
    """model : scipy distribution object
        A frozen scipy model object.  Needs to have the attribute *.shape
    """

_obs_pred_doc = \
    """obs, pred : array-like objects
        Observed and predicted data
    """


@doc_sub(_data_doc, _model_doc)
def nll(data, model):
    """
    Calculate the neagtive log likelihood given data and a model

    Parameters
    ----------
    {0}
    {1}

    Returns
    -------
    float
        Negative log likelihood

    """

    try:
        log_lik_vals = model.logpmf(data)
    except:
        log_lik_vals = model.logpdf(data)
    return -np.sum(log_lik_vals)


@doc_sub(_data_doc)
def lrt(data, model_null, model_alt, df=None):
    """
    This functions compares two nested models using the likelihood ratio
    test.

    Parameters
    ----------
    {0}
    model_null : scipy distribution object
        The null model as a frozen scipy distribution object. Parameters of
        distribution must be given as keyword arguments.
        Ex. ``norm = stats.norm(loc=0, scale=1)``

    model_alt : scipy distribution object
        The alternative model as a a frozen scipy distribution object.

    df : int
        Optional. Specify the degrees of freedom for the lrt.  Calculated
        as the number of parameters in model_alt - number of parameters in
        model_null.  If None, the df is calculated from the model
        objects.

    Returns
    -------
    tuple
        (G^2 statistic, p-value)

    Notes
    -----

    Interpretation: p-value < alpha suggests signficant evidence for your
    alternative model

    The LRT only applies to nested models. The variable test_stat is known as
    the G^2 statistic.  The G-test uses the fact that -2log(Likelihood_null /
    Likelihood_alt) is approximately chi-squared.  This assumption breaks down
    for small samples sizes.

    """

    # Calculate G^2 statistic
    ll_null = nll(data, model_null) * -1
    ll_alt = nll(data, model_alt) * -1
    test_stat = -2 * (ll_null - ll_alt)

    # Set df if necessary
    if not df:
        df = len(model_alt.kwds) - len(model_null.kwds)

    return (test_stat, stats.chisqprob(test_stat, df))


@doc_sub(_data_doc, _model_doc)
def AIC(data, model, params=None, corrected=True):
    """
    Calculate AIC given values of a model given data and model parameters

    Parameters
    ----------
    {0}
    {1}
    params : int
        The number of parameters in the model. If None, calculates the number
        of parameters from the distribution object

    corrected : bool
        If True, calculates the corrected AICC, if False calculates the
        uncorrected AIC.

    Returns
    -------
    float
        AIC(C) value

    Notes
    -----
    AICC should be used when the number of observations is < 40.

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
        k = len(model.kwds)
    else:
        k = params

    if corrected:
        aic_value = 2 * k + 2 * L + (2 * k * (k + 1)) / (n - k - 1)
    else:
        aic_value = 2 * k + 2 * L

    return aic_value


def AIC_weights(aic_list):
    """
    Calculates the AIC weights for a given set of models.

    Parameters
    -----------------
    aic_list : array-like object
        Array-like object containing AIC values from different models

    Returns
    -------------
    tuple
        First element contains the relative AIC weights, second element
        contains the delta AIC values.

    Notes
    -----
    AIC weights can be interpreted as the probability that a given model is the
    best model in comparison to the other models
    """

    aic_values = np.array(aic_list)
    minimum = np.min(aic_values)
    delta = aic_values - minimum
    values = np.exp(-delta / 2)
    weights = values / np.sum(values)

    return weights, delta


@doc_sub(_obs_pred_doc)
def sum_of_squares(obs, pred):
    """
    Calculates the sum of squares between observed (X) and predicted (Y) data.
    Attempts to braodcast arrays if lengths don't match.

    Parameters
    ----------
    {0}
    Returns
    -------
    float
        Sum of squares
    """
    obs, pred = tuple(np.broadcast_arrays(obs, pred))
    return np.sum((np.array(obs) - np.array(pred)) ** 2)


@doc_sub(_obs_pred_doc)
def r_squared(obs, pred, one_to_one=False, log_trans=True):
    """
    Get's the R^2 value for a regression of observed (X) and predicted (Y)
    data

    Parameters
    ----------
    {0}
    one_to_one : bool
        If True, calculates the R^2 based on the one-to-one line as done in
        [#]_.  If False, calculates the standard R^2 from a regression fit.

    log_trans : bool
        If True, log transforms obs and pred.

    Returns
    -------
    float
        R^2 value

    Notes
    -----
    Using just R^2 to compare the fit of observed and predicted values can be
    misleading because the relationship may not be one-to-one but the R^2
    value may be quite high. The one-to-one option alleviates this problem.

    References
    ----------
    .. [#]
        White, E., Thibault, K., & Xiao, X. (2012). Characterizing the species
        abundance distributions across taxa and ecosystems using a simple
        maximum entropy model. Ecology, 93(8), 1772-8

    """

    # Sort obs and pred
    obs = np.sort(obs)
    pred = np.sort(pred)

    if log_trans:
        obs = np.log(obs)
        pred = np.log(pred)

    if one_to_one:
        # Equation from White et al 2012
        r_sq = 1 - sum_of_squares(obs, pred) / \
                        sum_of_squares(obs, np.mean(obs))
    else:
        b0, b1, r, p_value, se = stats.linregress(obs, pred)
        r_sq = r ** 2

    return r_sq


def bin_data(data, max_num):
    """
    Bins the data on base 2.  Uses Preston's method of binning which has
    exclusive lower boundaries and inclusive upper boundaries. Densities are
    not split between bins.

    Parameters
    ----------
    data : array-like
        Data to be binned

    max_num :  float
        The maximum upper most boundary of the data

    Returns
    -------
    tuple
        (binned_data, bin_edges)

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
