"""
===========================
Main (:mod:`macroeco.main`)
===========================

This module contains functions that execute macroecological analyses specified 
by user-generated `parameters.txt` configuration files. Instructions for 
creating parameter files can be found here.

.. autosummary::
   :toctree: generated/

   main

"""

from __future__ import division
import os
import shutil
import inspect
import configparser

from pandas import DataFrame
import matplotlib.pyplot as plt
from matplotlib.mlab import rec2csv, rec_append_fields
from mpltools import style
style.use('ggplot')

from twiggy_setup import get_log
from empirical import *
from distributions2 import *
from compare import *


# Dictionary with keys for allowable metrics and func type
metric_types = {
    'sad': 'dist',
    'sar': 'curve',
    'ear': 'curve',
    'ssad': 'dist',
}


def main(param_path='parameters.txt'):
    """
    Entry point function for analysis based on parameter files.

    Parameters
    ----------
    param_dir : str
        Path to directory containing user-generated parameter file

    """
  
    # Confirm file is present and extract dir name
    # TODO: Because of log catch in twiggy_setup, this doesn't print anything
    if not os.path.isfile(param_path):
        raise IOError, "Parameter file not found at %s" % param_path
    param_dir = os.path.dirname(param_path)
    
    # Get logger and announce start
    log = get_log(param_dir, clear=True)
    log.info('Starting analysis')

    # Read parameter file into params object
    params = configparser.ConfigParser()
    params.read(param_path)

    # Do analysis for each run with options dict (params + addl options)
    run_names = params.sections()
    for run_name in run_names:
        log.info('Starting run %s' % run_name)
        options = dict(params[run_name])
        options['param_dir'] = os.path.abspath(param_dir)
        options['run_dir'] = os.path.join(param_dir, run_name)
        options['metric_type'] = _check_metric(options)
        _do_analysis(options)
    log.info('Finished analysis successfully')


def _check_metric(options):
    """
    Checks if metric is in options list and returns string for metric type.

    Parameters
    ----------
    options : dict
        Option names and values for analysis

    Returns
    -------
    str
        'dist' for distribution, 'curve' for curve type, None if no metric is 
        specified.
    """
    if not 'metric' in options:
        return None
    try:
        return metric_types[options['metric']]
    except Exception:
        raise NotImplementedError, ("No analysis for metric %s is currently "
                                    "possible." % options['metric'])

def _do_analysis(options):
    """
    Do analysis for a single run, as specified by options.

    Parameters
    ----------
    options : dict
        Option names and values for analysis

    """

    if 'metric' in options:
        emp_results = _analyze_empirical(options)
    else:
        emp_results = None

    if 'models' in options:
        mod_results = _analyze_models(options, emp_results)
    else:
        mod_results = None

    _save_results(options, emp_results, mod_results)


def _analyze_empirical(options):
    """
    Perform empirical analysis of metric on data set

    Parameters
    ----------
    options : dict
        Option names and values for analysis

    Returns
    -------
    list of tuples
        Each tuple corresponds to a combination (see XXX), with first element 
        of the tuple giving a dictionary describing the combination and the 
        second element giving the result of the analysis. Any additional 
        elements are not used.

    """
    # TODO: (In empirical) Create result objects rather than strange lists of
    # nested tuples.
    
    # If no metadata path is given or data path invalid, raise error
    metadata_path = os.path.normpath(os.path.join(options['param_dir'], 
                                                  options['metadata']))
    if not os.path.isfile(metadata_path):
        raise IOError, "Path to metadata file %s is invalid." % metadata_path

    # Create Patch object for this data
    patch = Patch(metadata_path)

    # Get cols and splits variable (req by all metrics) and add to options
    options['cols'], options['splits'] = _get_cols_splits(options)

    # Get names of args and kwargs to method specified by metric option
    exec ("arg_and_kwd_names, _, _, kw_defaults = "
          "inspect.getargspec(%s)" % options['metric'])
    if kw_defaults:
        arg_names = arg_and_kwd_names[1:-len(kw_defaults)]  # Ignore patch
    kw_names = arg_and_kwd_names[-len(kw_defaults):]
    else:
        arg_names = arg_and_kwd_names[1:]
        kw_names = []

    # Create list with vals for all args and dict with vals for all kwargs
    # All required args must be in options
    args = []  # Patch is always first argument
    for arg_name in arg_names:
        try:
            exec 'args.append(eval("%s"))' % options[arg_name]
        except SyntaxError: # eval failing because option is a string
            exec 'args.append("%s")' % options[arg_name]
        except:
            raise ValueError, ("Value for required argument %s not provided"
                               % arg_name)

    kwargs = {}
    for kw_name in kw_names:
        if kw_name in options.keys():
        try:
                exec 'kwargs[kw_name]=eval("%s")' % options[kw_name]
            except SyntaxError:  # eval failing because option is a string
                exec 'kwargs[kw_name]="%s"' % options[kw_name]
            except:
                raise ValueError, ("Value for optional argument %s is invalid" 
                                   % kw_name)

    # Call Patch method with appropriate args and return result
    return eval("%s(patch, *args, **kwargs)" % options['metric'])


def _get_cols_splits(options):
    """
    Notes
    -----
    Always returns strings, even if dictionary or list is constructed here, to 
    ensure consistency with provided options.

    """

    # Splits may be given as option, else is set to empty
    if 'splits' in options.keys():
        splits = options['splits']
    else:
        splits = None

    # Cols may be given as option or individual col options may be options
    if 'cols' in options.keys():
        cols = options['cols']
    else:
        cols = {}
        for col in ['spp_col', 'count_col', 'energy_col', 'mass_col']:
            cols[col] = options.get(col, None)

    return str(cols), str(splits)


def _analyze_models(options, emp_results):
    """
    Perform theoretical analysis based on empirical data or options

    Parameters
    ----------
    options : dict
        Option names and values for analysis
    emp_results : list of tuples
        Output of method of `empirical.Patch`, or None if no data given

    Returns
    -------
    list of tuples
        Each tuple corresponds to a combination in emp_result, with one element 
        in each tuple for the result of each model comparison. The result 
        object is another tuple of fitted parameters (tuple), values (array), 
        comparison statistic names (list), and comparison statistic values 
        (list).

    """
    
    if emp_results:
        mod_results = _analyze_models_from_data(options, emp_results)
    else:
        mod_results = _analyze_models_from_options(options)

    return mod_results


def _analyze_models_from_data(options, emp_results):
    """
    Perform model analysis based on empirical data

    Parameters
    ----------
    options : dict
        Option names and values for analysis
    emp_results : list of tuples
        Output of method of `empirical.Patch`

    Returns
    -------
    list of dicts
        Each dict in the list corresponds to the similarly indexed combination 
        in emp_result. Dicts have a key for each given model name, with values 
        that are a four element list of fitted parameters (tuple), values 
        (array), comparison statistic names (tuple), and comparison statistic 
        values (tuple).

    """

    # Get list of model names
    models = options['models'].replace(' ', '').split(';')

    # Fit theories to all emp_results
    # TODO: Make work for 2D results, i.e., curves, comm_sep, o_ring
    # TODO: Make work for curves in general
    output_all = []
    for emp_result in emp_results:
        output_emp_result = {}
        for model in models:
            data = emp_result[1]['y'].values
            fits = _get_fits(data, model)
            values = _get_values(data, model, fits)
            stat_names, stats = _get_comparison_statistic(values, fits)
            output_emp_result[model] = [fits, values, stat_names, stats]
        output_all.append(output_emp_result)

    return output_all


def _analyze_models_from_options(options):
    """
    Perform model analysis based on options

    Parameters
    ----------
    options : dict
        Option names and values for analysis

    Returns
    -------
    list of tuples
        List of length 1 containing 1 tuple of length 1 (parallel structure to 
        _analyze_models_with_data). Content of that tuple is fitted parameters 
        (tuple).

    """
    raise NotImplementedError, "Models cannot be analyzed without data"

    #_get_fits_from_options should call model.translate_args (if exists)


def _get_fits(data, model):
    return eval("%s.fit2(data)" % model)


def _get_values(data, model, fits):

    try:
        values = eval("%s.pdf(data, *fits)" % model)
    except AttributeError:
        values = eval("%s.pmf(data, *fits)" % model)
    except:
        pass
    
    return values

def _get_comparison_statistic(data, fits):
    return ['AIC'], [get_AIC(data, fits)]


def _save_results(options, emp_results, mod_results):
    """
    Save results of analysis as tables and figures

    Parameters
    ----------
    options : dict
        Option names and values for analysis
    emp_results : list
        Results of empirical metric analysis from _analyze_empirical
    mod_results : list
        Results of theoretical metric analysis from _analyze_theoretical

    """

    # Ensure that output dir for this run exists and is empty
    shutil.rmtree(options['run_dir'], ignore_errors=True)
    os.makedirs(options['run_dir'])

    # Write outputs depending on pres/abs of emp and mod and dist/curve metric
    _write_combination_index_file(options, emp_results)
    _write_output(options, emp_results, mod_results)


def _write_combination_index_file(options, emp_results):
    """
    Write index of combinations table, giving number and combination
    """
    
    if not emp_results:
        return None

    f_path = os.path.join(options['run_dir'], '_combination_index.csv')
    with open(f_path, 'a') as f:
        for i,emp_result in enumerate(emp_results):
            f.write("%i,%s\n" % (i+1, str(emp_result[0])))


def _write_output(options, emp_results, mod_results):
    """
    Three groups of output
    - Fitted params (always if there is a model)
    - Data and pred (always if there is data, although no pred if no models)
    - Test statistis (only if both data and model)
    """

    # Get combinations from either emp or mod - if both exist must be same
    try:
        n_combs = len(emp_results)
    except:
        n_combs = len(mod_results)

    # Get list of names of models
    try:
        models = options['models'].replace(' ','').split(",")
    except:
        models = None        

    # Loop through all combinations
    for cidx in range(n_combs):
        if mod_results:
            _write_fitted_params(cidx, models, options, mod_results)
        if emp_results:
            _write_and_plot_data_pred(cidx, models, options, emp_results, 
                                      mod_results)
        if mod_results and emp_results:
            _write_test_statistics(cidx, models, options, mod_results)


def _write_fitted_params(cidx, models, options, mod_results):

    f = open(_get_file_path(cidx, options, "fitted_params.csv"), 'w')
    f.write("Model, Fit Parameters\n")

    for model in models:
        mod_result = mod_results[cidx][model]
        mod_fits = str(mod_result[0])[1:-1]  # Drop parens around tuple
        f.write("%s,%s\n" % (model, mod_fits))
    f.close()


def _write_and_plot_data_pred(cidx, models, options, emp_results, mod_results):
    """
    For distributions, will write and plot three kinds of comparisons
    - pdf/pmf vs histogram
    - cdf vs emp cdf
    - rad vs rad

    For curves, we'll only do data vs pred (note will have x and y values)
    """

    if options['metric_type'] == 'dist':
        _data_pred_dist(cidx, models, options, emp_results, mod_results)
    elif options['metric_type'] == 'curve':
        _data_pred_curve(cidx, models, options, emp_results, mod_results)


def _data_pred_dist(cidx, models, options, emp_results, mod_results):
    """
    These tables have column for data and each model.
    - pdf/pmf vs histogram
    - cdf vs emp cdf
    - rad vs rad
    Also make plots for all three
    """

    emp_result = emp_results[cidx][1]['y'].values
    n_vals = len(emp_result)

    # CDF
    # TODO: This goes up by integers to max value, can be too large
    x, emp_cdf = get_empirical_cdf(emp_result)

    def calc_func(model, x, shapes):
        return eval("%s.cdf(x, *shapes)" % model)

    plot_exec_str = "ax.step(x, emp, color='k')"

    _save_table_and_plot(cidx, models, options, mod_results, 'data_pred_cdf', 
                         x, emp_cdf, calc_func, plot_exec_str)

    # RAD
    x = np.arange(n_vals)/float(n_vals) + 0.5/float(n_vals)
    emp_rad = np.sort(emp_result)[::-1]

    def calc_func(model, x, shapes):
        return eval("%s.ppf(x, *shapes)" % model)[::-1]

    plot_exec_str = "ax.step(x * x_plot_mult, emp, color='k')"

    _save_table_and_plot(cidx, models, options, mod_results, 'data_pred_rad', 
                         x, emp_rad, calc_func, plot_exec_str, 
                         x_plot_mult=n_vals)

    # PDF/PMF
    hist_bins = 11
    emp_hist, edges = np.histogram(emp_result, hist_bins, normed=True)
    x = (np.array(edges[:-1]) + np.array(edges[1:])) / 2

    def calc_func(model, x, shapes):
        try:
            return eval("%s.pmf(np.floor(x), *shapes)" % model)
        except:
            return eval("%s.pdf(x, *shapes)" % model)

    plot_exec_str = "ax.bar(x-width/2, emp, width=width, color='gray')"

    _save_table_and_plot(cidx, models, options, mod_results, 'data_pred_pdf', 
                         x, emp_hist, calc_func, plot_exec_str)


def _save_table_and_plot(cidx, models, options, mod_results, name, x, emp, 
                         calc_func, plot_exec_str):

    f_path = _get_file_path(cidx, options, '%s.csv' % name)
    p_path = _get_file_path(cidx, options, '%s.png' % name)


    df = DataFrame({'x': x})
    df['empirical'] = emp
    for model in models:
        mod_result = mod_results[cidx][model]
        shapes = mod_result[0]
        result = calc_func(model, x, shapes)
        df[model] = result

    df.to_csv(f_path, index=False, float_format='%.4f')  # Table

    df_plt = df.set_index('x')  # Figure
    emp = df_plt['empirical']
    df_plt = df_plt.drop('empirical',1)

    width = x[1] - x[0]
    ax = df_plt.plot()
    exec plot_exec_str
    ax = _pad_plot_frame(ax)
    fig = ax.get_figure()
    fig.savefig(p_path)

    plt.close('all')

def _pad_plot_frame(ax, pad=0.01):
    """
    Provides padding on sides of frame equal to pad fraction of plot
    """

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xrange = xmax - xmin
    yrange = ymax - ymin

    ax.set_xlim(xmin - xrange*pad, xmax + xrange*pad)
    ax.set_ylim(ymin - yrange*pad, ymax + yrange*pad)

    return ax


def _data_pred_curve(cidx, models, options, emp_results, mod_results):
    raise NotImplementedError, "Data and curve comparison not implemented"


def _write_test_statistics(cidx, models, options, mod_results):
    # TODO: Add delta test statistics columns

    f = open(_get_file_path(cidx, options, "test_statistics.csv"), 'w')

    # Gets stat name list from any element of result dict - same for all models
    stat_names_list = next(mod_results[cidx].itervalues())[2]
    stat_names_str = str(stat_names_list)[1:-1].strip("'")

    f.write("Theory, %s\n" % stat_names_str)

    for model in models:
        mod_result = mod_results[cidx][model]
        mod_stats = str(mod_result[3])[1:-1]
        f.write("%s,%s\n" % (model, mod_stats))
    f.close()


def _get_file_path(cidx, options, file_name):
    return os.path.join(options['run_dir'],
                        '%i_%s' % (cidx+1, file_name))


if __name__ == '__main__':
    main(sys.argv[1])
