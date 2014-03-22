"""
Set up logging
"""

import twiggy
import traceback
import sys
import os
import threading as thread


def get_log(log_dir, clear=False):
    """
    Set up and return logger object
    """

    # Get path to log file and clear if requested
    log_path = os.path.join(log_dir,'log.txt')
    if clear and os.path.isfile(log_path):
        os.remove(log_path)
    
    # Get outputs and add emitters
    file_output, std_output = _logger_outputs(log_path)
    twiggy.addEmitters(('file', twiggy.levels.DEBUG, None, file_output), 
                       ('stdout', twiggy.levels.INFO, None, std_output))

    # Get logger
    # TODO: Once modules are in subdirs, change to __name__
    log = twiggy.log.name('meco')

    # Log uncaught exceptions (must occur after log declared)
    def log_uncaught(type1, value1, traceback1):
        tb_list = traceback.format_exception(type1, value1, traceback1)
        tb_str = ''.join(tb_list)
        log.options(suppress_newlines=False).critical('\n'+tb_str)
    sys.excepthook = log_uncaught

    # Make threads use sys.excepthook from parent process
    _installThreadExcepthook()

    return log


def _logger_outputs(log_path):

    # To ensure that Macroeco Desktop captures stdout, we just print it
    class stdLineFormat(twiggy.formats.LineFormat):
        def __call__(self, msg):
            text = self.format_text(msg)
            print "{text}".format(**locals())
            return ""
   
    # Choose formats for file and stdout
    file_format = twiggy.formats.LineFormat(traceback_prefix='')
    std_format = stdLineFormat(traceback_prefix='')

    # Set up outputs for file and stdout and create emitters
    file_output = twiggy.outputs.FileOutput(log_path, format=file_format)
    std_output = twiggy.outputs.StreamOutput(format=std_format)

    return file_output, std_output


def _installThreadExcepthook():
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
    docstring of `cls`s method.

    From scipy.misc.doccer

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
