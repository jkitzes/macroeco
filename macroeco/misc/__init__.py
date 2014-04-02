from .misc import *
from .rcparams import *
from format_data import (format_columnar, format_dense, format_grid,
                        format_transect)
=======
"""
===============================
Misc (:mod:`macroeco.misc`)
===============================

This module contains miscellaneous functions that support the functions of
other modules of macroeco.

.. autosummary::
   :toctree: generated/

   setup_log
   inherit_docstring_from
   doc_sub
   log_start_end

"""

from .misc import (setup_log, _thread_excepthook, log_start_end,
                   inherit_docstring_from, doc_sub)
from .rcparams import ggplot_rc

_thread_excepthook()  # Make desktop app catch and log sys except from thread
