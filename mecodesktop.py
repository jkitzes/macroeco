"""
MacroecoDesktop script for making standalone executable
"""

import sys as _sys
from macroeco import desktop

if len(_sys.argv) > 1:
    desktop(_sys.argv[1])
else:
    desktop()
