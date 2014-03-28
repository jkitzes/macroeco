"""
Set up logging
"""

import sys
import os
import traceback
import threading as thread

import twiggy
from twiggy import log
log = log.name('meco')
import decorator
import time

def setup_log(log_dir, file_name='_log.txt', clear=False):
    """
    Set up and return logger object
    """

    # Get path to log file and clear if requested
    log_path = os.path.join(log_dir, file_name)
    if clear and os.path.isfile(log_path):
        os.remove(log_path)

    # Get outputs and add emitters
    file_output, std_output = _logger_outputs(log_path)
    twiggy.addEmitters(('file', twiggy.levels.DEBUG, None, file_output),
                       ('stdout', twiggy.levels.INFO, None, std_output))

    # Get logger
    log = twiggy.log.name('meco')

    # Log uncaught exceptions (must occur after log declared)
    def log_uncaught(type1, value1, traceback1):
        tb_list = traceback.format_exception(type1, value1, traceback1)
        tb_str = ''.join(tb_list)
        log.options(suppress_newlines=False).critical('\n'+tb_str)
    sys.excepthook = log_uncaught

    return log


def _logger_outputs(log_path):

    # std_format - to ensure Macroeco Desktop shows logging, we just print
    class stdLineFormat(twiggy.formats.LineFormat):
        def __call__(self, msg):
            text = self.format_text(msg)
            print "{text}".format(**locals())
            return ""
    std_format = stdLineFormat(traceback_prefix='')

    # file_format - customized to show local time, etc
    conversion = twiggy.lib.converter.ConversionTable()
    conversion.add("time", _logger_better_time, "[{1}]".format)
    conversion.add("name", str, "{{{1}}}".format)
    conversion.add("level", str, "{1}".format)
    conversion.aggregate = ' '.join
    conversion.genericValue = str
    conversion.genericItem = "{0}={1}".format

    file_format = twiggy.formats.LineFormat(traceback_prefix='', separator=' ',
                                            conversion=conversion)

    # Set up outputs for file and stdout and create emitters
    file_output = twiggy.outputs.FileOutput(log_path, format=file_format)
    std_output = twiggy.outputs.StreamOutput(format=std_format)

    return file_output, std_output


def _logger_better_time(gmtime=None):
    return time.strftime("%Y/%m/%d %H:%M:%S %p", time.localtime())


def _thread_excepthook():
    """
    Make threads use sys.excepthook from parent process
    http://bugs.python.org/issue1230540
    """
    init_old = thread.Thread.__init__
    def init(self, *args, **kwargs):
        init_old(self, *args, **kwargs)
        run_old = self.run
        def run_with_except_hook(*args, **kw):
            try:
                run_old(*args, **kw)
            except (KeyboardInterrupt, SystemExit):

                raise
            except:
                sys.excepthook(*sys.exc_info())
        self.run = run_with_except_hook
    thread.Thread.__init__ = init


def inherit_docstring_from(cls):
    """
    This decorator modifies the decorated function's docstring by
    replacing occurrences of '%(super)s' with the docstring of the
    method of the same name from the class `cls`.

    If the decorated method has no docstring, it is simply given the
    docstring of cls method.

    Extracted from scipy.misc.doccer.

    """
    def _doc(func):
        cls_docstring = getattr(cls, func.__name__).__doc__
        func_docstring = func.__doc__
        if func_docstring is None:
            func.__doc__ = cls_docstring
        else:
            new_docstring = func_docstring % dict(super=cls_docstring)
            func.__doc__ = new_docstring
        return func
    return _doc


def doc_sub(*sub):
    """
    Decorator for performing substitutions in docstrings.

    Using @doc_sub(some_note, other_note) on a function with {0} and {1} in the
    docstring will substitute the contents of some_note and other_note for {0}
    and {1}, respectively.

    Decorator appears to work properly both with IPython help (tab completion
    and ?) and with Sphinx.

    """
    def dec(obj):
        obj.__doc__ = obj.__doc__.format(*sub)
        return obj
    return dec

def log_start_end(f):
    """
    Decorator to log start and end of function

    Use of decorator module here ensures that argspec will inspect wrapped
    function, not the decorator itself.
    http://micheles.googlecode.com/hg/decorator/documentation.html
    """
    def inner(f, *args, **kwargs):
        log.info('Starting %s' % f.__name__)
        res = f(*args, **kwargs)
        log.info('Finished %s' % f.__name__)
        return res
    return decorator.decorator(inner, f)
