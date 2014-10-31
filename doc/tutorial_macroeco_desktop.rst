.. _using-macroecodesktop:

=====================
Using MacroecoDesktop
=====================

This tutorial describes the basic usage of the the high-level MacroecoDesktop interface. Users who wish to use the ``macroeco`` Python package should refer to the :ref:`using-macroeco` tutorial.

This tutorial builds on the :ref:`first-steps-macroeco-desktop` tutorial, which should be completed first.

There are three basic types of analysis that can be completed using MacroecoDesktop: analysis of an empirical ecological pattern, fitting models to an empirical pattern, and exploration of a model of a macroecological pattern without empirical data. This tutorial describes these three types of analysis in turn.

Analyzing Empirical Patterns
============================

The first step in analyzing an empirical data set is to prepare the table and metadata file for the empirical data set as described in the :ref:`own-data` tutorial. It is generally most convenient to place the data table (usually in csv format) and the metadata file in the same folder.

The second step is to prepare a parameters file to describe the desired analysis. A parameter file has a section for each run that is part of an analysis. Each run is independent of the others, and multiple runs may be combined in a single parameter file for convenience (for example, several analyses may be run on a single data set, or a single metric may be calculated for many data sets).

An example of a run performing a species area analysis for the demo data set is shown below. ::

    [SAR]
    analysis = ssad

    metadata = ANBO.txt
    subset = column >= 2
    cols = spp_col:spp; count_col:count; x_col:row; y_col:column
    splits = row:2
    divs = 1,1; 1,2; 2,1; 2,2
    ear = False

Each run begins with a title that is enclosed in square brackets. The run title should not have any spaces in it.

The first section following the run name contains a single variable ``analysis`` which gives the name of an empirical ecological pattern to analyze for this data set. The available empirical analyses in Macroeco |version| are

.. currentmodule:: macroeco.empirical
.. autosummary::
   :toctree: generated/

   sad
   ssad
   sar
   comm_grid
   o_ring

The pages linked above from each analysis name describe the analysis and the different input parameters required to complete the analysis. Each of these input parameters can be specified here in the parameters file.

For example, examining the ``sar`` metric above shows that this analysis takes five input parameters:

* patch - a Patch object containing the empirical data
* cols - a string associating column headers in the data table with their uses
* splits - a string describing how and whether to split the data into multiple separate subsets before analysis
* divs - the number of divisions along coordinate columns
* ear - True or False, where True calculates an endemics area relationship and False calculates a species area relationship

Each of these five input parameters is provided directly in the run shown above, with the exception of the ``patch`` parameter, which is described slightly differently. Although the descriptions below apply to the ``sar`` metric, many of the same input parameters are used by the other analysis metrics.

In all empirical analyses in Macroeco, the first input parameter is a patch object. Instead of describing this object directly in MacroecoDesktop, the user first provides the ``metadata`` and, optionally, the ``subset`` input parameters.

The first parameter, ``metadata``, gives the relative path to the metadata file from the parameter file (if the parameter file and metadata file are in the same folder, this is just the name of the metadata file).

The second parameter, ``subset``, takes a subset of the empirical data for further analysis. Any logical mathematical statement beginning with a column name and ending with a number can be used here. For example, setting ``subset`` to ``year==2010; row < 2; spp=='cabr'`` would perform all subsequent calculations only for data in which the year column is 2010, the row column is greater than 2, and the species column is equal to 'cabr'. Multiple conditions are separated by semicolons. In the example run above, the SAR will be calculated only for columns 2 and 3 of the data.

The next input parameter for an SAR analysis is ``cols``, which is a string describing which column in the data table should be used for which "special columns" in analysis. The five possible special columns are

- spp_col - Unique species identifiers
- count_col - Number of individuals at a location
- x_col - x coordinate of location
- y_col - y coordinate of location
- energy_col - Energetic requirements of individual(s) at a location

Analyses that do not have a spatial component (like a species abundance distribution without subsets or splits) require only spp_col and count_col (if one exists - if not, each record is taken to represent one individual). Spatial analyses, such as the species-area relationship, also require x_col and y_col. Energy metrics require energy_col.

The ``cols`` parameter can also be set in the Description section of the metadata file, in which case it is not required here. If ``cols`` is set both in a run and in a metadata file, the value from the run takes precedence.

The next input parameter is ``splits``, which provides a convenient way to divide a data set into separate analyses. The value ``year:split; row:2:``, for example, would split the data set into unique years and also into two subplots along the row column, each of equal length. The value before the ``:`` is a column name, and the value after is either a number (if a numeric column is to be split into equal sized divisions) or the word "split" (if a column is to be split among all unique values).

This parameter is particularly useful if a column for plot ID, family name, functional group, etc. is present in the data table, in which case splitting on that column would perform an identical analysis for each different plot, family, group, etc. ``splits`` can also be used, for example, to split a plot into four subplots along two coordinate axes and perform a species area analysis within each one.

The next input parameter is ``divs``, which gives the number of divisions to perform along the x and y columns. For example, ``3:2;`` will divide a plot into six subplots, with three "columns" formed by splitting the x axis into three parts and two "rows" formed by splitting the y axis into two parts. Here, ``1,1; 1,2; 2,1; 2,2`` will analyze the species area relationship for the entire plot, half plots (split in both directions), and quarter plots.

The final input parameter, ``ear``, determines whether a species area or an endemics area relationship should be calculated. This is a boolean value that can be either True (endemics area relationship) or False (species area relationship).

Once the parameters file has been created and saved, it can be executed using MacroecoDesktop by following the instructions at the end of the :ref:`first-steps-macroeco-desktop` tutorial.

A sample parameter file containing runs that complete many of the above empirical analyses can be found in :ref:`recipes`.


Fitting and Comparing Models of Empirical Patterns
==================================================

NacroecoDesktop can also be used to fit models to empirical data patterns, analyze the goodness of fit of these models, and to compare the fits of multiple models. This process is identical to that described above for analyzing empirical patterns, except that one additional set of input parameters is added to a run in the parameters file. ::

    [SAR]
    analysis = ssad

    metadata = ANBO.txt
    subset = column >= 2
    cols = spp_col:spp; count_col:count; x_col:row; y_col:column
    splits = row:2
    divs = 1,1; 1,2; 2,1; 2,2
    ear = False

    models = power_law
    log_x = true
    log_y = true

The third portion of this run begins with the input parameter ``models``, which can be set equal to one or several of the models within the ``macroeco`` package. If the metric is a curve, such as the species area relationship, the following models may be used.

.. currentmodule:: macroeco.models
.. autosummary::
   :toctree: generated/

   power_law
   mete_sar
   mete_sar_iterative
   mete_ear

If the metric is a probability distribution, the following models may be used (note that some are discrete and some continuous).

.. autosummary::
   :toctree: generated/

   expon
   expon_uptrunc
   lognorm
   geom
   geom_uptrunc
   nbinom
   cnbinom
   logser_uptrunc
   plnorm
   plnorm_ztrunc

More information about these models can be found by clicking on their names above. Some of these models have additional optional parameters that can be provided here (see the Methods section of the page for each individual model).

Two special input parameters, ``log_x`` and ``log_y``, are used to log transform the x and y axes of the output plots created by MacroecoDesktop.

As another example, the run below will calculate a species abundance distribution for the demo data set, fit both a lognormal and upper-truncated logseries distribution to the empirical data, and compare their fits. ::

    [SAD]
    analysis = sad

    metadata = ANBO.txt

    models = logser_uptrunc; lognorm
    log_y = True

Note that no subsets or splits are given here, so that the entire data table is used for the analysis. The ``cols`` parameter is also not given, and the value of this parameter from the metadata file is used as a result.

Exploring Models
================

Finally, MacroecoDesktop may also be used to explore the behavior of models without specific reference to empirical data. Given a set of model parameters, the "y" values of curves may be calculated for any "x" values, and the probability density, cumulative density, random variates, and many other values may be calculated for probability distributions.

To see the possible options for exploring models, choose a model from the lists above and refer to the Methods section of that page. Any model and any method may be used with MacroecoDesktop so long as all of the input parameters required by that method are provided. Note that although the ``loc`` and ``scale`` parameters are listed for distributions, these are not used by Macroeco and should not be entered in a parameters file.

For example, the parameter file below contains runs that calculate the pmf of a geometric distribution with a known shape parameter ``p``, calculate the ``p`` parameter of the upper-truncated geometric distribution from the distribution mean and aggregation parameter ``k``, fit the parameters of a lognormal distribution to a small data set, and draw 20 random variables from a conditioned negative binomial distribution. ::

    [Geom-pmf]
    analysis = geom.pmf

    x = 0,1,2,3,4,5
    p = 0.5

    [GeomUptrunc-p]
    analysis = geom_uptrunc.translate_args

    mu = 5
    b = 20

    [Lognorm-fit]
    analysis = lognorm.fit_mle

    data = 2,2,5,8,4,3

    [Cnbinom-random]
    analysis = cnbinom.rvs

    mu = 10
    k_agg = 2
    b = 15
    size = 10

