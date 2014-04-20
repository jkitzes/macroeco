"""
===============================================
Macroeco: Ecological pattern analysis in Python
===============================================

Macroeco provides a comprehensive set of functions for analyzing empirical
patterns in data, predicting patterns using theory and models, and comparing
empirical results to theory. Many major macroecological patterns can be
analyzed using this package, including the species abundance distribution, the
species and endemics area relationships, several measures of beta diversity,
and many others.

Extensive documentation for macroeco, including tutorials and a reference
guide, are available at http://macroeco.org.

The package is organized into five submodules.

Empirical provides a Patch class for reading data and metadata from an
empirical census and functions that calculate empirical macroecological metrics
based on that data.

Models provides a set of distributions and curves that have been proposed by
basic theory to describe macroecological metrics.

Compare provides functions for comparing the empirical and modeled results.

Misc provides a set of miscellanous functions, including several that aid in
formatting census data for use by functions in the empirical module.

Main provides a programmatic interface to this package, known as Macroeco
Desktop, that allows a user to specify all of the parameters for an analysis in
a single parameters file, which is then executed, and results saved, with no
additional intervention needed.

Macroeco was developed at the University of California, Berkeley, by Justin
Kitzes and Mark Wilber. Additional contributors include Chloe Lewis and Ethan
White. The development of macroeco has been supported by the National Science
Foundataion, the Moore Foundation, and the Berkeley Institute for Global Change
Biology.

"""

import sys as _sys

__version__ = '0.3'

import empirical
import models
import compare
import main
import misc

def mecodesktop():
    if len(_sys.argv) > 1:
        param_path = _sys.argv[1]
        main.main(param_path)
    else:
        print "Macroeco Desktop must be called with path to parameters file"