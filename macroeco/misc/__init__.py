"""
===============================
Misc (:mod:`macroeco.misc`)
===============================

This module contains miscellaneous functions that support the functions of
other modules of macroeco.

Support Functions
=================

.. autosummary::
   :toctree: generated/

   setup_log
   log_start_end
   inherit_docstring_from
   doc_sub

"""
"""

Data Formatting Functions
=========================

.. autosummary::
   :toctree: generated/

   data_read_write
   format_dense

"""

from .misc import (setup_log, _thread_excepthook, log_start_end,
                   inherit_docstring_from, doc_sub)
from .rcparams import ggplot_rc
from .format_data import (data_read_write, format_dense)

_thread_excepthook()  # Make desktop app catch and log sys except from thread
