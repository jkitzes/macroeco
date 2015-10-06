"""
===============================================
Macroeco: Ecological pattern analysis in Python
===============================================

Macroeco provides a comprehensive set of functions for analyzing empirical
patterns in ecological data, predicting patterns using theory and models, and
comparing empirical patterns to theory.  Many major macroecological patterns
can be analyzed using this package, including the species abundance
distribution, the species and endemics area relationships, several measures of
beta diversity, and many others.

Macroeco can be used either as a scientific python Package or through a high-
level interface called MacroecoDesktop. Users new to Macroeco should begin by
reviewing the tutorials found below. Experienced Python programmers who wish to
use the ``macroeco`` Python package can ``pip install macroeco`` and refer to
the :ref:`using-macroeco` tutorial and the :ref:`reference` guide.

.. toctree::
   :maxdepth: 2

   tutorials
   reference
   about

"""

import sys as _sys

__version__ = '1.0.1'

import empirical
import models
import compare
import main
import misc

def desktop(param_path=False):
    if param_path:
        main.main(param_path)
    else:
        main.launch()
