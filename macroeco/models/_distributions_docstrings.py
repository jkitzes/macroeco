# From scipy 0.14 - docstring internals keep changing with scipy updates
# So for now, include it here for import into _distributions

# These are the docstring parts used for substitution in specific
# distribution docstrings

docheaders = {'methods': """\nMethods\n-------\n""",
              'parameters': """\nParameters\n---------\n""",
              'notes': """\nNotes\n-----\n""",
              'examples': """\nExamples\n--------\n"""}

_doc_rvs = """\
rvs(%(shapes)s, loc=0, scale=1, size=1)
    Random variates.
"""
_doc_pdf = """\
pdf(x, %(shapes)s, loc=0, scale=1)
    Probability density function.
"""
_doc_logpdf = """\
logpdf(x, %(shapes)s, loc=0, scale=1)
    Log of the probability density function.
"""
_doc_pmf = """\
pmf(x, %(shapes)s, loc=0, scale=1)
    Probability mass function.
"""
_doc_logpmf = """\
logpmf(x, %(shapes)s, loc=0, scale=1)
    Log of the probability mass function.
"""
_doc_cdf = """\
cdf(x, %(shapes)s, loc=0, scale=1)
    Cumulative density function.
"""
_doc_logcdf = """\
logcdf(x, %(shapes)s, loc=0, scale=1)
    Log of the cumulative density function.
"""
_doc_sf = """\
sf(x, %(shapes)s, loc=0, scale=1)
    Survival function (1-cdf --- sometimes more accurate).
"""
_doc_logsf = """\
logsf(x, %(shapes)s, loc=0, scale=1)
    Log of the survival function.
"""
_doc_ppf = """\
ppf(q, %(shapes)s, loc=0, scale=1)
    Percent point function (inverse of cdf --- percentiles).
"""
_doc_isf = """\
isf(q, %(shapes)s, loc=0, scale=1)
    Inverse survival function (inverse of sf).
"""
_doc_moment = """\
moment(n, %(shapes)s, loc=0, scale=1)
    Non-central moment of order n
"""
_doc_stats = """\
stats(%(shapes)s, loc=0, scale=1, moments='mv')
    Mean('m'), variance('v'), skew('s'), and/or kurtosis('k').
"""
_doc_entropy = """\
entropy(%(shapes)s, loc=0, scale=1)
    (Differential) entropy of the RV.
"""
_doc_fit = """\
fit(data, %(shapes)s, loc=0, scale=1)
    Parameter estimates for generic data.
"""
_doc_expect = """\
expect(func, %(shapes)s, loc=0, scale=1, lb=None, ub=None, conditional=False, **kwds)
    Expected value of a function (of one argument) with respect to the distribution.
"""
_doc_expect_discrete = """\
expect(func, %(shapes)s, loc=0, lb=None, ub=None, conditional=False)
    Expected value of a function (of one argument) with respect to the distribution.
"""
_doc_median = """\
median(%(shapes)s, loc=0, scale=1)
    Median of the distribution.
"""
_doc_mean = """\
mean(%(shapes)s, loc=0, scale=1)
    Mean of the distribution.
"""
_doc_var = """\
var(%(shapes)s, loc=0, scale=1)
    Variance of the distribution.
"""
_doc_std = """\
std(%(shapes)s, loc=0, scale=1)
    Standard deviation of the distribution.
"""
_doc_interval = """\
interval(alpha, %(shapes)s, loc=0, scale=1)
    Endpoints of the range that contains alpha percent of the distribution
"""
_doc_allmethods = ''.join([docheaders['methods'], _doc_rvs, _doc_pdf,
                           _doc_logpdf, _doc_cdf, _doc_logcdf, _doc_sf,
                           _doc_logsf, _doc_ppf, _doc_isf, _doc_moment,
                           _doc_stats, _doc_entropy, _doc_fit,
                           _doc_expect, _doc_median,
                           _doc_mean, _doc_var, _doc_std, _doc_interval])

# Note that the two lines for %(shapes) are searched for and replaced in
# rv_continuous and rv_discrete - update there if the exact string changes
_doc_default_callparams = """
Parameters
----------
x : array_like
    quantiles
q : array_like
    lower or upper tail probability
%(shapes)s : array_like
    shape parameters
loc : array_like, optional
    location parameter (default=0)
scale : array_like, optional
    scale parameter (default=1)
size : int or tuple of ints, optional
    shape of random variates (default computed from input arguments )
moments : str, optional
    composed of letters ['mvsk'] specifying which moments to compute where
    'm' = mean, 'v' = variance, 's' = (Fisher's) skew and
    'k' = (Fisher's) kurtosis. (default='mv')
"""
_doc_default_longsummary = """\
Continuous random variables are defined from a standard form and may
require some shape parameters to complete its specification.  Any
optional keyword parameters can be passed to the methods of the RV
object as given below:
"""
_doc_default_frozen_note = """
Alternatively, the object may be called (as a function) to fix the shape,
location, and scale parameters returning a "frozen" continuous RV object:

rv = %(name)s(%(shapes)s, loc=0, scale=1)
    - Frozen RV object with the same methods but holding the given shape,
      location, and scale fixed.
"""
_doc_default_example = """\
Examples
--------
>>> from scipy.stats import %(name)s
>>> import matplotlib.pyplot as plt
>>> fig, ax = plt.subplots(1, 1)

Calculate a few first moments:

%(set_vals_stmt)s
>>> mean, var, skew, kurt = %(name)s.stats(%(shapes)s, moments='mvsk')

Display the probability density function (``pdf``):

>>> x = np.linspace(%(name)s.ppf(0.01, %(shapes)s),
...               %(name)s.ppf(0.99, %(shapes)s), 100)
>>> ax.plot(x, %(name)s.pdf(x, %(shapes)s),
...          'r-', lw=5, alpha=0.6, label='%(name)s pdf')

Alternatively, freeze the distribution and display the frozen pdf:

>>> rv = %(name)s(%(shapes)s)
>>> ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')

Check accuracy of ``cdf`` and ``ppf``:

>>> vals = %(name)s.ppf([0.001, 0.5, 0.999], %(shapes)s)
>>> np.allclose([0.001, 0.5, 0.999], %(name)s.cdf(vals, %(shapes)s))
True

Generate random numbers:

>>> r = %(name)s.rvs(%(shapes)s, size=1000)

And compare the histogram:

>>> ax.hist(r, normed=True, histtype='stepfilled', alpha=0.2)
>>> ax.legend(loc='best', frameon=False)
>>> plt.show()
"""

_doc_default = ''.join([_doc_default_longsummary,
                        _doc_allmethods,
                        _doc_default_callparams,
                        _doc_default_frozen_note,
                        _doc_default_example])

_doc_default_before_notes = ''.join([_doc_default_longsummary,
                                     _doc_allmethods,
                                     _doc_default_callparams,
                                     _doc_default_frozen_note])

docdict = {
    'rvs': _doc_rvs,
    'pdf': _doc_pdf,
    'logpdf': _doc_logpdf,
    'cdf': _doc_cdf,
    'logcdf': _doc_logcdf,
    'sf': _doc_sf,
    'logsf': _doc_logsf,
    'ppf': _doc_ppf,
    'isf': _doc_isf,
    'stats': _doc_stats,
    'entropy': _doc_entropy,
    'fit': _doc_fit,
    'moment': _doc_moment,
    'expect': _doc_expect,
    'interval': _doc_interval,
    'mean': _doc_mean,
    'std': _doc_std,
    'var': _doc_var,
    'median': _doc_median,
    'allmethods': _doc_allmethods,
    'callparams': _doc_default_callparams,
    'longsummary': _doc_default_longsummary,
    'frozennote': _doc_default_frozen_note,
    'example': _doc_default_example,
    'default': _doc_default,
    'before_notes': _doc_default_before_notes
}

# Reuse common content between continuous and discrete docs, change some
# minor bits.
docdict_discrete = docdict.copy()

docdict_discrete['pmf'] = _doc_pmf
docdict_discrete['logpmf'] = _doc_logpmf
docdict_discrete['expect'] = _doc_expect_discrete
_doc_disc_methods = ['rvs', 'pmf', 'logpmf', 'cdf', 'logcdf', 'sf', 'logsf',
                     'ppf', 'isf', 'stats', 'entropy', 'expect', 'median',
                     'mean', 'var', 'std', 'interval']
for obj in _doc_disc_methods:
    docdict_discrete[obj] = docdict_discrete[obj].replace(', scale=1', '')
docdict_discrete.pop('pdf')
docdict_discrete.pop('logpdf')

_doc_allmethods = ''.join([docdict_discrete[obj] for obj in _doc_disc_methods])
docdict_discrete['allmethods'] = docheaders['methods'] + _doc_allmethods

docdict_discrete['longsummary'] = _doc_default_longsummary.replace(
    'Continuous', 'Discrete')
_doc_default_frozen_note = """
Alternatively, the object may be called (as a function) to fix the shape and
location parameters returning a "frozen" discrete RV object:

rv = %(name)s(%(shapes)s, loc=0)
    - Frozen RV object with the same methods but holding the given shape and
      location fixed.
"""
docdict_discrete['frozennote'] = _doc_default_frozen_note

_doc_default_discrete_example = """\
Examples
--------
>>> from scipy.stats import %(name)s
>>> import matplotlib.pyplot as plt
>>> fig, ax = plt.subplots(1, 1)

Calculate a few first moments:

%(set_vals_stmt)s
>>> mean, var, skew, kurt = %(name)s.stats(%(shapes)s, moments='mvsk')

Display the probability mass function (``pmf``):

>>> x = np.arange(%(name)s.ppf(0.01, %(shapes)s),
...               %(name)s.ppf(0.99, %(shapes)s))
>>> ax.plot(x, %(name)s.pmf(x, %(shapes)s), 'bo', ms=8, label='%(name)s pmf')
>>> ax.vlines(x, 0, %(name)s.pmf(x, %(shapes)s), colors='b', lw=5, alpha=0.5)

Alternatively, freeze the distribution and display the frozen ``pmf``:

>>> rv = %(name)s(%(shapes)s)
>>> ax.vlines(x, 0, rv.pmf(x), colors='k', linestyles='-', lw=1, 
...         label='frozen pmf')
>>> ax.legend(loc='best', frameon=False)
>>> plt.show()

Check accuracy of ``cdf`` and ``ppf``:

>>> prob = %(name)s.cdf(x, %(shapes)s)
>>> np.allclose(x, %(name)s.ppf(prob, %(shapes)s))
True

Generate random numbers:

>>> r = %(name)s.rvs(%(shapes)s, size=1000)
"""
docdict_discrete['example'] = _doc_default_discrete_example

_doc_default_before_notes = ''.join([docdict_discrete['longsummary'],
                                     docdict_discrete['allmethods'],
                                     docdict_discrete['callparams'],
                                     docdict_discrete['frozennote']])
docdict_discrete['before_notes'] = _doc_default_before_notes

_doc_default_disc = ''.join([docdict_discrete['longsummary'],
                             docdict_discrete['allmethods'],
                             docdict_discrete['frozennote'],
                             docdict_discrete['example']])
docdict_discrete['default'] = _doc_default_disc
