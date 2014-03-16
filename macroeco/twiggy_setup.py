import twiggy
import traceback
import sys
import os
import threading as thread

# Output format for log file - remove traceback prefix
file_format = twiggy.formats.LineFormat(traceback_prefix='')

# Output format for terminal logging - only text message part
class stdLineFormat(twiggy.formats.LineFormat):
    def __call__(self, msg):
        text = self.format_text(msg)
        print "{text}".format(**locals())
        return ""
std_format = stdLineFormat(traceback_prefix='')

# Logger setup - returns logger object
def get_log(log_dir='/Users/jkitzes/Desktop/', clear=False):

    # Get path to log file - must be writable (ie, not inside pyinstaller app)
    log_path = os.path.join(log_dir,'log.txt')

    # Delete log file if requested
    if clear:
        try:
            os.remove(log_path)
        except OSError:
            pass

    # Set up outputs for file and stdout
    file_output = twiggy.outputs.FileOutput(log_path, format=file_format)
    std_output = twiggy.outputs.StreamOutput(format=std_format, 
                                             stream=sys.stdout)

    # Create emitters
    twiggy.addEmitters(('file', twiggy.levels.DEBUG, None, file_output), 
                       ('stdout', twiggy.levels.INFO, None, std_output))

    # Declare logger for macroeco
    # TODO: Once modules are in subdirs, change to __name__ to log module also
    log = twiggy.log.name('meco')

    return log

# Log uncaught exceptions
log = twiggy.log.name('meco')  # If below called before log def elsewhere
def log_uncaught(type1, value1, traceback1):
    tb_list = traceback.format_exception(type1, value1, traceback1)
    tb_str = ''.join(tb_list)
    log.options(suppress_newlines=False).critical('\n'+tb_str)
sys.excepthook = log_uncaught

# Use proper excepthook for threads also
def installThreadExcepthook():
    """
    Workaround for sys.excepthook thread bug
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
installThreadExcepthook()


