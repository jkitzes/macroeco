from __future__ import division
import sys
import os
import shutil
import warnings
import inspect
import configparser
import threading as thread
import traceback
import copy
import logging

import numpy as np
import pandas as pd

import matplotlib as mpl
mpl.use('agg')  # Prevents crash when GUI runs matplotlib in thread on Linux
import matplotlib.pyplot as plt

from .. import empirical as emp
from .. import models as mod
from .. import compare as comp
from .. import misc

import time
def _better_time(gmtime=None):
    return 

def main(param_path='parameters.txt'):
    """
    Entry point function for analysis based on parameter files.

    Parameters
    ----------
    param_path : str
        Path to user-generated parameter file

    """

    # Confirm parameters file is present
    if not os.path.isfile(param_path):
        raise IOError, "Parameter file not found at %s" % param_path

    # Get raw params and base options (non-run-dependent options)
    params, base_options = _get_params_base_options(param_path)

    # Configure and start logging
    # Done here instead of in function so will affect all subsequent calls
    log_path = os.path.join(base_options['results_dir'], '_log.txt')
    if os.path.isfile(log_path):
        os.remove(log_path)

    logging.basicConfig(level=logging.INFO, format='%(message)s')

    fileh = logging.FileHandler(log_path)
    fileh.setLevel(logging.DEBUG)
    filefmt = logging.Formatter(
            time.strftime("%Y/%m/%d %H:%M:%S %p", time.localtime()) + 
            ' - %(name)s - %(levelname)s - %(message)s')
    fileh.setFormatter(filefmt)
    logging.getLogger('').addHandler(fileh)

    def log_uncaught(type1, value1, traceback1):
        tb_list = traceback.format_exception(type1, value1, traceback1)
        tb_str = ''.join(tb_list)
        logging.critical('\n\n'+tb_str)
    sys.excepthook = log_uncaught

    logging.info('Running macroeco') # v%s' % __version__)
    logging.info('Parameters file at %s' % os.path.abspath(param_path))

    # Preliminary check for errors in parameters file
    bad_params = misc.check_parameter_file(param_path)
    if len(bad_params[0]) > 0:
        logging.warning("Possible formatting error(s) in" +
                    " %s: parameters %s on lines %s"
                      % (param_path, bad_params[0], bad_params[1]))

    logging.info('Starting analysis')

    # Do analysis for each run
    for run_name in base_options['run_names']:
        logging.info('Starting run %s' % run_name)
        options = dict(params[run_name])  # All parameters from this run
        options.update(base_options)  # Add base parameters
        options['run_dir'] = os.path.join(base_options['results_dir'],run_name)
        if 'format' in options['analysis']:
            _do_format(options)
        else:
            _do_analysis(options)
        logging.info('Finished run %s' % run_name)
    logging.info('Finished analysis successfully')
    logging.info('Results available at %s' % options['param_dir'])

    # Close logging - releases log file lock in Windows GUI
    logging.shutdown()
    

def _get_params_base_options(param_path):

    # Read parameter file into params object
    params = configparser.ConfigParser()
    try:
        params.read(param_path)
    except:
        raise ValueError, "Parameter file is invalid"

    # Setup param_dir and results_dir, get run_names
    param_dir = os.path.abspath(os.path.dirname(param_path))
    results_dir = os.path.join(param_dir, 'results')

    if os.path.isdir(results_dir):
        shutil.rmtree(results_dir)
    os.makedirs(results_dir)

    run_names = params.sections()

    # Check there's at least one run
    if not run_names:
        raise NameError, "Parameters file must contain at least one run"

    # Create options dict
    base_options = {}
    base_options['param_dir'] = param_dir
    base_options['results_dir'] = results_dir
    base_options['run_names'] = run_names

    return params, base_options


def _do_format(options):

    datapath = os.path.normpath(os.path.join(options['param_dir'],
                                                  options['data']))
    out_path = os.path.splitext(datapath)[0] + "_formatted.csv"

    format_type = options['analysis'].split('_')[1]
    misc.data_read_write(datapath, out_path, format_type, **options)


def _do_analysis(options):
    """
    Do analysis for a single run, as specified by options.

    Parameters
    ----------
    options : dict
        Option names and values for analysis

    """

    module = _function_location(options)
    core_results = _call_analysis_function(options, module)

    if module == 'emp' and ('models' in options.keys()):
        fit_results = _fit_models(options, core_results)
    else:
        fit_results = None

    _save_results(options, module, core_results, fit_results)


def _function_location(options):
    # TODO: Add spec and misc modules
    # This relies on the assumption that there are no duplicate member names
    # in the different modules.
    func_name = options['analysis'].split('.')[0]  # Ignore method if present
    emp_members = [x[0] for x in inspect.getmembers(emp)]
    mod_members = [x[0] for x in inspect.getmembers(mod)]
    if func_name in emp_members:
        module = 'emp'
    elif func_name in mod_members:
        module = 'mod'
    else:
        raise ValueError, ("No analysis of type '%s' is available" %
                           options['analysis'])
    return module


def _call_analysis_function(options, module):
    """
    Call function from module and get result, using inputs from options

    Parameters
    ----------
    options : dict
        Option names and values for analysis
    module : str
        Short name of module within macroeco containing analysis function

    Returns
    -------
    dataframe, array, value, list of tuples
        Functions from emp module return a list of tuples in which first
        element of the tuple gives a string describing the result and the
        second element giving the result of the analysis as a dataframe.
        Functions in other modules return dataframe, array, or value.

    """

    args, kwargs = _get_args_kwargs(options, module)
    return eval("%s.%s(*args, **kwargs)" % (module, options['analysis']))


def _get_args_kwargs(options, module):
    """
    Given an options (including analysis), and module, extract args and kwargs
    """

    if module == 'emp':
        options = _emp_extra_options(options)
    arg_names, kw_names = _arg_kwarg_lists(module, options['analysis'])

    # Create list of values for arg_names
    args = []
    for arg_name in arg_names:

        if arg_name == 'patch':  # For patch arg, append actual patch obj
            args.append(options['patch'])
            continue
        if arg_name == 'self':  # Ignore self from class methods
            continue
        if arg_name == 'k':  # scipy dists use k and x, we always use x
            arg_name = 'x'

        try:
            exec 'args.append(eval("%s"))' % options[arg_name]
        except SyntaxError: # eval failing because option is a string
            args.append(options[arg_name])
        except:
            raise ValueError, ("Value for required argument %s not provided"
                               % arg_name)

    # Create dict with vals for kw_names
    kwargs = {}
    for kw_name in kw_names:
        if kw_name in options.keys():  # If a value is given for this kwarg
            try:
                exec 'kwargs[kw_name] = eval("%s")' % options[kw_name]
            except SyntaxError:  # eval failing because value is a string
                kwargs[kw_name] = options[kw_name]
            except:
                raise ValueError, ("Value for optional argument %s is invalid"
                                   % kw_name)

    return args, kwargs


def _emp_extra_options(options):
    """
    Get special options patch, cols, and splits if analysis in emp module
    """

    # Check that metadata is valid
    metadata_path = os.path.normpath(os.path.join(options['param_dir'],
                                                  options['metadata']))
    if not os.path.isfile(metadata_path):
        raise IOError, ("Path to metadata file %s is invalid." %
                        metadata_path)
    options['metadata_path'] = metadata_path

    # Using subset if given, create and store patch
    subset = options.get('subset', '')
    options['patch'] = emp.Patch(metadata_path, subset)

    # If cols or splits not given in options, make empty strings
    if 'cols' not in options.keys():
        options['cols'] = ''
    if 'splits' not in options.keys():
        options['splits'] = ''

    return options


def _arg_kwarg_lists(module, analysis):

    # Get names of args and kwargs to method specified by analysis option
    exec ("arg_and_kwd_names, _, _, kw_defaults = "
          "inspect.getargspec(%s.%s)" % (module, analysis))
    if kw_defaults:  # If there are kwargs
        arg_names = arg_and_kwd_names[:-len(kw_defaults)]
        kw_names = arg_and_kwd_names[-len(kw_defaults):]
    else:  # If no kwargs
        arg_names = arg_and_kwd_names
        kw_names = []

    # Inspection for rv classes doesn't work since it uses args internally
    # Unless method is translate_args or fit_mle, appends shapes to args
    try:
        obj_meth = analysis.split('.')
        if obj_meth[1] not in ['fit_mle', 'translate_args']:
            arg_names += eval(module + '.' + obj_meth[0] + '.' +
                              "shapes.replace(' ','').split(',')")
        if obj_meth[1] == 'rvs':  # Inspection for size not working
            kw_names.append('size')
    except:
        pass

    return arg_names, kw_names


def _fit_models(options, core_results):
    """
    Fit models to empirical result from a function in emp module

    Parameters
    ----------
    options : dict
        Option names and values for analysis
    core_results : list of tuples
        Output of function in emp

    Returns
    -------
    list of dicts
        Each element in list corresponds to a subset. The dict has a key for
        each model given in options, and the value is a list of fitted
        parameters (tuple), values (array), comparison statistic names (list),
        and comparison statistic values (list).

    Notes
    -----
    To determine if the empirical result refers to a curve or a distribution,
    the result dataframe is inspected for a column 'x', which indicates a
    curve.

    """

    logging.info("Fitting models")
    models = options['models'].replace(' ', '').split(';')

    # TODO: Make work for 2D results, i.e., curves, comm_sep, o_ring
    # TODO: Make work for curves in general (check if 'x' present in core_res)
    fit_results = []
    for core_result in core_results:  # Each subset
        fit_result = {}
        for model in models:
            fits = _get_fits(core_result, model, options)
            values = _get_values(core_result, model, fits)
            stat_names, stats = _get_comparison_stat(core_result, values,
                                                     model, fits)
            fit_result[model] = [fits, values, stat_names, stats]
        fit_results.append(fit_result)

    return fit_results


def _get_fits(core_result, model, options):

    options_copy = {}
    for key, val in options.iteritems():
        if key not in ['patch']:  # Ignore patch since won't deepcopy
            options_copy[key] = copy.deepcopy(val)

    model_obj = eval('mod.' + model)
    if hasattr(model_obj, 'fit_mle'):
        options_copy['analysis'] = model + '.' + 'fit_mle'
        options_copy['data'] = core_result[1]['y'].values
    else:
        options_copy['analysis'] = model + '.' + 'fit_lsq'
        options_copy['x'] = core_result[1]['x'].values
        options_copy['y_obs'] = core_result[1]['y'].values
        options_copy['df'] = core_result[1]  # Entire result df, for mete_sar

    return _call_analysis_function(options_copy, 'mod')


def _get_values(core_result, model, fits):

    model_obj = eval('mod.' + model)
    if hasattr(model_obj, 'vals'):
        x = core_result[1]['x'].values  # Calc model at x values
        values = eval("mod.%s.vals(x, *fits)" % model)
    else:
        n = len(core_result[1])  # Calc model at data values
        values = eval("mod.%s.rank(n, *fits)" % model)

    return values


def _get_comparison_stat(core_result, values, model, fits):
    # Uses AIC for distributions, R2 one-to-one for curves

    try:  # Only curves have vals
        eval("mod.%s" % model + ".vals.__doc__")
        obs = core_result[1]['y'].values
        pred = values
        name = ['R2']
        stat = comp.r_squared(obs, pred, one_to_one=True)
    except AttributeError:
        obs = core_result[1]['y'].values
        name = ['AIC']
        stat = comp.AIC(obs, eval("mod.%s" % model + "(*fits)"))

    return name, stat


def _save_results(options, module, core_results, fit_results):
    """
    Save results of analysis as tables and figures

    Parameters
    ----------
    options : dict
        Option names and values for analysis
    module : str
        Module that contained function used to generate core_results
    core_results : dataframe, array, value, list of tuples
        Results of main analysis
    fit_results : list or None
        Results of comparing emp analysis to models, None if not applicable

    """

    logging.info("Saving all results")

    # Use custom plot format
    mpl.rcParams.update(misc.rcparams.ggplot_rc)

    # Make run directory
    os.makedirs(options['run_dir'])

    # Write core results
    _write_core_tables(options, module, core_results)

    # Write additional results if analysis from emp
    if module == 'emp':
        _write_subset_index_file(options, core_results)

    # Write model/data comparison if models were given
    if fit_results:
        models = options['models'].replace(' ','').split(';')
        for i, core_result in enumerate(core_results):
            _write_fitted_params(i, models, options, fit_results)
            _write_test_statistics(i, models, options, fit_results)
            _write_comparison_plot_table(i, models, options,
                                         core_results, fit_results)

def _write_core_tables(options, module, core_results):
    """
    Notes
    -----
    Depending on function that was called for analysis, core_results may be a
    list of tuples (empirical), a dataframe, an array, or a single value.

    For the list of tuples from empirical, the second element of each tuple is
    the raw result, and we write them all with the appropriate prefix. For
    dataframes, we write them. For arrays or single values, we convert to data
    frames and write them.

    """

    table_name = 'core_result.csv'
    single_file_path = os.path.join(options['run_dir'], table_name)

    if module == 'emp':  # List of tuples
        for i, core_result in enumerate(core_results):
            file_path = _get_file_path(i, options, table_name)
            core_result[1].to_csv(file_path, index=False, float_format='%.4f')

    elif type(core_results) == type(pd.DataFrame()):  # DataFrame
        core_results.to_csv(single_file_path, index=False, float_format='%.4f')

    else:  # Array or single value (atleast_1d corrects for unsized array)
        df = pd.DataFrame({'y': np.atleast_1d(core_results)})
        df.to_csv(single_file_path, index=False, float_format='%.4f')


def _get_file_path(spid, options, file_name):
    return os.path.join(options['run_dir'],
                        '%i_%s' % (spid+1, file_name))


def _write_subset_index_file(options, core_results):
    """
    Write table giving index of subsets, giving number and subset string
    """

    f_path = os.path.join(options['run_dir'], '_subset_index.csv')
    subset_strs = zip(*core_results)[0]
    index = np.arange(len(subset_strs)) + 1
    df = pd.DataFrame({'subsets': subset_strs}, index=index)
    df.to_csv(f_path)


def _write_fitted_params(spid, models, options, fit_results):
    # TODO: Consider converting to pandas, need to deal with variable length
    # TODO: Possibility - empty data frame max length, max width = nparams
    f = open(_get_file_path(spid, options, 'fitted_params.csv'), 'w')
    f.write("Model, Fit Parameters\n")

    for model in models:
        fit_result = fit_results[spid][model]
        mod_fits = str(fit_result[0])[1:-1]  # Drop parens around tuple
        f.write("%s,%s\n" % (model, mod_fits))
    f.close()


def _write_test_statistics(spid, models, options, fit_results):
    # TODO: Add delta test statistics columns
    # TODO: Make dataframe?
    f = open(_get_file_path(spid, options, 'test_statistics.csv'), 'w')

    # Gets stat name list from any element of result dict - same for all
    stat_names_list = next(fit_results[spid].itervalues())[2]
    stat_names_str = str(stat_names_list)[1:-1].strip("'")

    f.write("Model, %s\n" % stat_names_str)

    for model in models:
        fit_result = fit_results[spid][model]
        fit_stats = str(fit_result[3])[:]
        f.write("%s,%s\n" % (model, fit_stats))
    f.close()


def _write_comparison_plot_table(spid, models, options, core_results,
                                 fit_results):
    """
    Notes
    -----
    Only applies to analysis using functions from empirical in which models are
    also given.

    """
    # TODO: Clean up sorting, may not work if SAR x out of order, e.g.

    is_curve = 'x' in core_results[0][1]
    df = core_results[spid][1]
    df.rename(columns={'y': 'empirical'}, inplace=True)

    # If distribution, need to sort values so will match sorted rank in fits
    if not is_curve:
        x = np.arange(len(df)) + 1
        df = df.sort(columns='empirical')
        df.insert(0, 'x', x[::-1])

    # Add residual column for each model
    for model in models:
        fit_result = fit_results[spid][model]
        df[model] = fit_result[1]
        df[model + "_residual"] = df[model] - df['empirical']

    # If curve, sort now for plotting purposes
    if is_curve:
        df = df.sort(columns='x')

    # Set up file paths
    f_path = _get_file_path(spid, options, 'data_models.csv')
    p_path = _get_file_path(spid, options, 'data_models.pdf')

    # Save table
    df.to_csv(f_path, index=False, float_format='%.4f')  # Table

    # Save plot
    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax1.scatter(df['x'], df['empirical'], color='k')
    ax1.plot(df['x'], df[models])
    ax1.legend(models + ['empirical'], loc='best')
    ax1.set_xlabel('x')
    ax1.set_ylabel('value')

    ax2.hlines(0, np.min(df['x']), np.max(df['x']))
    ax2.plot(df['x'], df[[x + '_residual' for x in models]])
    ax2.legend(models + ['empirical'], loc='best')
    ax2.set_xlabel('x')
    ax2.set_ylabel('residual')
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(min(ax2.get_ylim()[0], -1), max(ax2.get_ylim()[1], 1))

    if options.get('log_y', None):
        ax1.set_yscale('log')
        ax2.set_yscale('symlog', linthreshy=1)
    if options.get('log_x', None):
        ax1.set_xscale('log')
        ax2.set_xscale('log')

    if not options.get('log_x', None) and not options.get('log_y', None):
        ax1.set_ylim(bottom=0)
        ax1.set_xlim(left=0)
        ax1 = _pad_plot_frame(ax1)
        ax2 = _pad_plot_frame(ax2)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()
    fig.savefig(p_path)

    plt.close('all')


def _pad_plot_frame(ax, pad=0.01):
    """
    Provides padding on sides of frame equal to pad fraction of plot
    """

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xr = xmax - xmin
    yr = ymax - ymin

    ax.set_xlim(xmin - xr*pad, xmax + xr*pad)
    ax.set_ylim(ymin - yr*pad, ymax + yr*pad)

    return ax


def _output_cdf_plot(core_result, spid, models, options, fit_results):
    """Function for plotting cdf"""

    # CDF
    x = core_result['y'].values
    df = emp.empirical_cdf(x)
    df.columns = ['x', 'empirical']

    def calc_func(model, df, shapes):
        return eval("mod.%s.cdf(df['x'], *shapes)" % model)

    plot_exec_str = "ax.step(df['x'], emp, color='k', lw=3);ax.set_ylim(top=1)"

    _save_table_and_plot(spid, models, options, fit_results, 'data_pred_cdf',
                         df, calc_func, plot_exec_str)


def output_pdf_plot(core_result, spid, models, options, fit_results):
    """ Function for plotting pdf/pmf """
    # PDF/PMF
    hist_bins = 11
    emp_hist, edges = np.histogram(core_result['y'].values, hist_bins,
                                   normed=True)
    x = (np.array(edges[:-1]) + np.array(edges[1:])) / 2
    df = pd.DataFrame({'x': x, 'empirical': emp_hist})

    def calc_func(model, df, shapes):
        try:
            return eval("mod.%s.pmf(np.floor(df['x']), *shapes)" % model)
        except:
            return eval("mod.%s.pdf(df['x'], *shapes)" % model)

    plot_exec_str = "ax.bar(df['x']-width/2, emp, width=width, color='gray')"

    _save_table_and_plot(spid, models, options, fit_results, 'data_pred_pdf',
                         df, calc_func, plot_exec_str)
