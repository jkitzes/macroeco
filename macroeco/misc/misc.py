"""
Miscellaneous functions
"""

import sys
import os
import traceback
import threading as thread

import logging
import decorator
import time


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
        logging.info('Starting %s' % f.__name__)
        res = f(*args, **kwargs)
        logging.info('Finished %s' % f.__name__)
        return res
    return decorator.decorator(inner, f)


def check_parameter_file(filename):
    """
    Function does a rudimentary check whether the cols, splits and divs columns
    in the parameter files are formatted properly.

    Just provides a preliminary check. Will only catch basic mistakes

    Parameters
    ----------
    filename : str
        Path to parameters file

    Returns
    -------
    : list
        Contains the number of possible bad strings detected
    """

    # Load file
    with open(filename, "r") as fin:
        content = fin.read()

    # Check cols and splits strings

    bad_names = []
    line_numbers = []

    strs = ["cols", "splits", "divs"]

    for tstr in strs:

        start = content.find(tstr)

        while start != -1:

            cols_str = "".join(content[start:].split("\n")[0].split("=")[-1].split(" "))

            semis = cols_str.count(";")

            # Get line number
            line_end = content.find("\n", start)
            line_number = content[:line_end].count("\n") + 1

            if tstr == "divs":
                colons = cols_str.count(",")
            else:
                colons = cols_str.count(":")

            if colons != (semis + 1):
                bad_names.append(tstr)
                line_numbers.append(line_number)

            start = content.find(tstr, start + 1)

    return bad_names, line_numbers
