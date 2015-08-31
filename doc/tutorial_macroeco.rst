.. _using-macroeco:

==============
Using macroeco
==============

This tutorial describes the basic usage of the ``macroeco`` Python package. Users who wish to use the high-level MacroecoDesktop interface should refer to the :ref:`using-macroecodesktop` tutorial.

The `models` subpackage
============================

The `models` subpackage contains a number of common statistical distributions and curves used in macroecological analyses.  A full list of the available models and their documentation can be found at :doc:`models`.

All the statistical distributions contained in `models` inherit from the `rv_discrete` or `rv_continuous` classes defined in the package `scipy.stats.distributions`. Therefore, all the distribution can be used in the same way as any distributions defined in `scipy.stats`.  For example, a Fisher logseries distribution, a common distribution used for species abundance distributions, is defined by one parameter `p` which can take values between 0 and 1. This distribution could defined as follows

    >>> import macroeco.models as md
    >>> logser_dist = md.logser(p=0.9)

The probability of observing a species with one individual is

    >>> logser_dist.pmf(1)
    0.39086503371292664

The probability of observing and species with 10 individuals or less (given by the cumulative distribution function) is

    >>> logser_dist.cdf(10)
    0.9201603889810761

Similarly, the following function call gives the same results as above

    >>> md.logser.pmf(1, 0.9)
    0.39086503371292664
    >>> md.logser.cdf(10, 0.9)
    0.9201603889810761

The distributions contained in the `models` package also contain some special functions such as `rank` which gives the rank abundance distribution for a given distribution. For example, the rank abundance distribution for the logseries distribution given above with 30 species in the community is given by

    >>> md.logser.rank(30, 0.9)
    array([  1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,
         1.,   2.,   2.,   2.,   2.,   2.,   3.,   3.,   3.,   4.,   4.,
         5.,   5.,   6.,   7.,   8.,  10.,  13.,  21.])

A community with 30 species following this logseries distribution is expected to have 15 species with one individual.

The `models` subpackage also contains objects known as curves. These consist of macroecological curves such as species area relationships (SAR) and endemic area relationships (EAR).  Currently, there are 4 implemented curves

* Power Law
* METE SAR/EAR with direct upscaling and downscaling
* METE SAR with iterative upscaling and downscaling

The METE SAR is a particular SAR that is described at length in the book **Maximum Entropy and Ecology: A Theory of Abundance, Distribution, and Energetics** by John Harte. It can be used to upscale and downscale species richness given knowledge of the total number of species (`S`) and the total number of individual (`N`) at some base area.

To estimate the expected number of species present in an area that is double the size of the base area and half the size of the base area given `S = 30` and `N = 1000` we can use the following code

    >>> # Number of species in base area
    >>> S = 30

    >>> # Number of individuals in base area
    >>> N = 1000

    >>> # A list of habitat areas
    >>> areas = [1, 2, 0.5]
    >>> md.mete_sar_iterative.vals([1, 2, 0.5], S, N, approx=True)
    array([ 30.        ,  36.15434332,  23.46676609])

For the parameter `areas`, the first number in the list (1 in this example) is always the base area (e.g. 50 ha, 2.5 m^2, 300 in^2), and the following numbers are additional areas at which to calculate species richness (2 and 0.5 in this example). Using the argument `approx=True` significantly speeds up the calculation and will tend to given very similar answers to `approx=False`.

    >>> md.mete_sar_iterative.vals([1, 2, 0.5], S, N, approx=False)
    array([ 30.        ,  36.15164239,  23.46953897])


Additional subpackages
=========================

In addition to the `models` package, the ``macroeco`` package contains two other main subpackages of interest:

* `empirical` - loads data tables and performs empirical analysis of macroecological metrics, such as the species abundance distribution and species area relationship (:doc:`empirical`)

* `compare` - provides utility functions for comparing the fit of models to empirical metrics, such as AIC weights and r-squared statistics (:doc:`compare`)

A common workflow involves loading data, calculating an empirical metric, fitting one or more models to the empirical metric, and evaluating the fit of the model to the metric.

A simple species abundance distribution analysis
================================================

The following example shows a simple species abundance distribution analysis for the demo data.

First, the ``Patch`` class from the empirical subpackage is used to create a Patch object that holds the data table and a metadata dictionary describing the data. ``Patch`` requires a path, absolute or relative, to a metadata file as a mandatory argument (see :ref:`own-data` for information on creating a metadata file for a new data set).

    >>> import macroeco as meco
    >>> pat = meco.empirical.Patch('~/Desktop/ANBO.txt')

The empirical subpackage contains a number of functions that operate on patch objects and return macroecological metrics. Here we'll use the function ``sad`` to calculate a species abundance distribution. The first argument is the patch object to use, the second is a string specifying which column has the species names (`spp_col`) and which, if any, has a count of individuals at a particular location (`count_col`), and the third is a string specifying how to split the data (see :doc:`empirical` for more information on input arguments).

    >>> sad = meco.empirical.sad(pat, 'spp_col:spp; count_col:count', '')

All functions for macroecological metrics return their results as a list of tuples. Each tuple has two elements

1. A string describing how the data were split

2. A result table with a column ``y`` (for univariate distributions like the species abundance distribution) or columns ``y`` and ``x`` (for curves such as a species area relationship) giving the results of the analysis.

Since the data were not split in this example, the list has only one tuple.  The result is

    >>> sad
    [('',        spp     y
    0    arsp1     2
    1     cabr    31
    2   caspi1    58
    3     chst     1
    4    comp1     5
    5     cran     4
    6     crcr    65
    7    crsp2    79
    8     enfa     1
    9     gnwe    41
    10   grass  1110
    11   lesp1     1
    12    magl     1
    13    mesp     6
    14    mobe     4
    15    phdi   210
    16   plsp1     1
    17    pypo    73
    18    sasp     2
    19    ticr   729
    20   unsh1     1
    21   unsp1    18
    22   unsp3     1
    23   unsp4     1)]

where the first element of the tuple is `''` (an empty string because no split occurred) and the second element in the tuple is a `pandas` DataFrame with two columns: 1) the species ID (`spp`) and 2) the abundance of each species (`y`).  The DataFrame itself can easily be extracted

    >>> sad_df = sad[0][1]

where we recognize that the DataFrame is the second element (index 1) of the first tuple in the list (index 0).  This notation will make more sense when we consider splitting the data below.

Any number of distributions from the models subpackage can be fit to the resulting empirical metric. The code below fits the two parameters of the upper truncated logseries distribution and uses the function ``AIC`` from the compare subpackage to calculate the AIC for this distribution and data.

    >>> # Fit the uptruncated logseries distribution to the empirical SAD
    >>> p, b = meco.models.logser_uptrunc.fit_mle(sad_df['y'])

    >>> # Look at the two parameters of the distribution
    >>> p, b
    (0.9985394369365049, 2445.0)

The upper truncated logseries distribution has two parameters: `p` which is a parameter that is a function of the mean of the distribution and `b` which is the maximum value the distribution can take.  Using these parameters, we can get an AIC value to determine the "goodness of fit" of the upper truncated logseries distribution to the empirical data.

    >>> # Get the AIC value
    >>> logser_aic = meco.compare.AIC(sad_df['y'], meco.models.logser_uptrunc(p, b))
    >>> logser_aic
    208.61902087378027

If you are using the `ipython` environment you can see the arguments that meco.compare.AIC takes using `meco.compare.AIC?`.  In short, the function takes in the data (in this case the species abundance distribution) and fitted model object and returns the AIC value.  Of course, AICs aren't very useful by themselves so let's compare the logseries fit to a broken distribution, another classic theoretical SAD.  This is equivalent to a zero-truncated negative binomial distrbution with aggregation parameter `k` equal to 1.

    >>> # Get Broken Stick AIC
    >>> mu, sigma = meco.models.lognorm.fit_mle(sad_df['y'])
    >>> broken_stick_aic = meco.compare.AIC(sad_df['y'], meco.models.nbinom_ztrunc(np.mean(sad_df['y']), 1))
    >>> broken_stick_aic
    274.27490655552322

We can see that the lower AIC for the logseries suggests that this is a more appropriate model for this SAD.

We could also visually compare these models using their rank abundance distributions.  We first generate the rank abundance distributions for the fitted logseries and the broken stick distributions and then plot it against the empirical data.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> logser_rad = meco.models.logser_uptrunc.rank(len(sad_df), p, b)
    >>> broken_stick_rad = meco.models.nbinom_ztrunc.rank(len(sad_df), np.mean(sad_df['y']), 1)

    >>> # Plot the empirical data. Note that [::-1] reverses the order of a vector
    >>> plt.semilogy(np.sort(sad_df['y'])[::-1])

    >>> # Plot the RAD of the models
    >>> plt.semilogy(logser_rad[::-1])
    >>> plt.semilogy(broken_stick_rad[::-1])
    >>> plt.show()

A simple species-area relationship analysis
===========================================

We can also analyze species-area relationships (SAR)s using `macroeco`. To get an empirical SAR from the ANBO data we use the function `meco.empirical.sar`.  As described in the documentation, this function takes 4 key arguments

1. `patch`: The empirical Patch object

2. `cols`: A semicolon-separated column string that identifies the species column (i.e. `spp_col`, which column contains the species names), the count column (i.e. `count_col`, which column contains the species counts), the x column (i.e. `x_col`, which column specifies the spatial location of an individual in the x direction),  and the y column (i.e. `y_col`, which column specifies the spatial location of an individual in the y direction).  For example, this string for the ANBO data would be `'spp_col:spp; count_col:count; x_col:row; y_col:column'` because the column that contains the species names is `spp`, the column that contains the counts is `count`, the column that contain the spatial location of an individual in the x direction is `row` and the the column that contains the spatial location of an individual in the y direction is `column`.  For the SAR analysis, `x_col` and `y_col` must be specified.

3. `splits`: A string specifying whether the analysis should be run on different subsets of the data. For example, if one had a column `year` specifying different years that the community census was completed the string `year:split` would run the analysis on each year separately. `split` is a key word described in the documentation.

4. `divs`: A semicolon-separated string that describes how to successively divide the patch along the `x_col` and `y_col` dimensions. For example, the string `'1,2; 2,2; 2,4'` will calculate the average species richness at three areas. The first areas (1,2) will be made by splitting the y column in half.  The second areas (2, 2) will be made by splitting the x column and the y column in half.  The third areas (2, 4) will be made by splitting the x column in half and the y column in fourths.

The ANBO plot is 4m x 4m = 16 m^2. The get the SAR for the areas 1m^2, 2m^2, 4m^2, 8m^2, and 16m^2 we use the following code.

    >>> sar = meco.empirical.sar(pat, 'spp_col:spp; count_col:count; x_col:row; y_col:column', "", "1,1; 1,2; 2,1; 2,2; 2,4; 4,2; 4,4")
    >>> sar
    [('',    div  n_individs    n_spp   x        y
      0  1,1   2445.0000  24.0000  16  24.0000
      1  1,2   1222.5000  18.5000   8  18.5000
      2  2,1   1222.5000  17.0000   8  17.0000
      3  2,2    611.2500  13.5000   4  13.5000
      4  2,4    305.6250  10.1250   2  10.1250
      5  4,2    305.6250  10.5000   2  10.5000
      6  4,4    152.8125   7.5625   1   7.5625)]

The output of the SAR function is a list of tuples where each tuple is a particular split.  Because we did not split the data (i.e. the `split` parameter was `''`), we have one tuple.  The second item in this tuple is a `pandas` DataFrame that contains the key results of the analysis

    >>> sar_table = sar[0][1]
    >>> sar_table
       div  n_individs    n_spp   x        y
    0  1,1   2445.0000  24.0000  16  24.0000
    1  1,2   1222.5000  18.5000   8  18.5000
    2  2,1   1222.5000  17.0000   8  17.0000
    3  2,2    611.2500  13.5000   4  13.5000
    4  2,4    305.6250  10.1250   2  10.1250
    5  4,2    305.6250  10.5000   2  10.5000
    6  4,4    152.8125   7.5625   1   7.5625


The column `div` gives the divisions specified in the function call. The column `n_individs` specifies the average number of individuals across the cells made from the given division. `n_spp` gives the average species across the cells made from the given division. `x` gives the absolute area of the given division. `y` gives the same information as `n_spp` and is included for easy plotting.

For plotting, one might want to combine like areas to a single value and then plot.

    >>> # Combine similar areas
    >>> combined_sar = sar_table.groupby('x').mean().reset_index()

    >>> # Plot the SAR
    >>> plt.plot(combined_sar['x'], combined_sar['y'], '-o')
    >>> plt.xlabel("Area")
    >>> plt.ylabel("Species")
    >>> plt.show()

If we want to compare the empirical SAR to a power law SAR and a METE SAR we can first fit each of these curves to the data.  To fit the METE SAR, we only need the total number of species (`n_spp`) and total number of individuals (`n_individs`) at the base scale (i.e. at `div = 1,1`). We could either look at the table at see that `n_spp` at `div = 1,1` is 24 and `n_individs` is 2445 or pass in the data frame to the `fit_lsq` method of the `mete_sar` curve

    >>> # Fit the METE SAR
    >>> S0, N0 = meco.models.mete_sar_iterative.fit_lsq(sar_table)

    >>> # Get the predicted values from the fitted METE SAR
    >>> pred_mete = meco.models.mete_sar_iterative.vals(combined_sar['x'][::-1], S0, N0, approx=True)

We can fit a power law SAR using similar notation

    >>> # Fit the power law
    >>> c, z = meco.models.power_law.fit_lsq(combined_sar['x'], combined_sar['y'])

    >>> # Get the predicted value from the fitted power law
    >>> pred_power_law = meco.models.power_law.vals(combined_sar['x'][::-1], c, z)

and then compare these theoretical SARs to the empirical SAR

    >>> plt.plot(combined_sar['x'][::-1], np.array([pred_power_law, pred_mete]).T)

Clearly the power law SAR provides a better fit to the data than the METE SAR.  We can confirm this quantitatively using the R^2 value from the one to one line when we compare observed and predicted values.  If the predicted SAR is a perfect fit to the observed SAR, the predicted values will exactly equal the observed values (i.e. fall along the one to one line).

    >>> r2_mete = meco.compare.r_squared(combined_sar['y'][::-1], pred_mete, one_to_one=True, log_trans=True)
    >>> r2_mete
    0.64318610592964909

    >>> r2_power_law = meco.compare.r_squared(combined_sar['y'][::-1], pred_power_law, one_to_one=True, log_trans=True)
    >>> r2_power_law
    0.99939083620342017

The R^2 close to 1 for the power law supports the result from the plots that the power law is a better fitting SAR to the data.

A simple spatial analysis
==========================

We also might want to analyze the spatial patterns of individuals in the plot. We can get the spatial patterns of all the species in plot by using the `meco.empirical.ssad` function.

The SSAD is a species-level spatial abundance distribution.  In other words, how are the individuals of a species distributed in space? The empirical SSAD function has three arguments. The first is the Patch object, the second is the `cols` string, and the third is the split string specifying how to grid a given landscape.

For example, the split string `'row:4; column:4'` says to divide the column `row` into 4 equally spaced sections and divide the column `column` into 4 equally spaced sections.  This gives a grid with 16 equally sized cells.

We can do this for the ANBO data

    >>> all_spp_ssads = meco.empirical.ssad(pat, 'spp_col:spp; count_col:count', 'row:4; column:4')

`all_spp_ssads` is a list with 24 tuples where each tuple contains two items.  The first item is a string giving a species name and the second item is a data frame giving the abundance of the given species in each of the 16 cells.

    >>> all_spp_ssads[0]
        ('arsp1',     y
     0   0
     1   0
     2   0
     3   0
     4   0
     5   0
     6   0
     7   0
     8   1
     9   0
     10  0
     11  0
     12  0
     13  1
     14  0
     15  0)

We can loop through all of the species in `all_spp_ssads` and fit a finite negative binomial distribution to each species. The `k` parameter of this distribution specifies how aggregated a species is in space with `k` approaching 0 being very aggregated and `k` approaching infinity being binomially distributed.::


    # Store the results
    k_res = {}

    # Loop through all species
    for spp_name, data in all_spp_ssads:

        # Fit negative
        k_param = meco.models.cnbinom.fit_mle(data['y'], k_array=np.linspace(0.01, 5, num=1000))[1]

        # Get total abundance for a given species
        total_abund = data['y'].sum()

        # Store k parameter and total abundance for each species
        k_res[spp_name] = (k_param, total_abund)


The dictionary `k_res` contains the `k` parameter and total abundance for each species in the ANBO data.


A more complex analysis
=========================

One of the major benefits of `macroeco` is that you can explore how macroeco logical patterns vary across scale and/or for different subsets of your data. For example, what if we wanted to compare how an SAD changed across scale?  We will again use the ANBO data to illustrate this example

The ANBO metadata is given below ::

    [Description]
    name = Anzo Borrego
    author = Mary Ellen Harte and John Harte
    description = Vegetation census conducted at Anza-Borrego Desert State Park. Site in Indian Valley at N 32' 52.091", W 116' 14.447". Elevation 1195 feet. Census was conducted on a 4 m x 4 m grid, with 16 grid cells each 1 m2 in area.
    citation = Unpublished

    datapath = ANBO.csv
    cols = spp_col:spp; count_col: count; x_col: row; y_col: column

    [year]
    description = Year of census

    [cell]
    description = Unique cell identifier, from 0 to 15 (total of 16 cells)

    [row]
    description = Row of cell in gridded plot
    min = 0
    max = 3
    step = 1

    [column]
    description = Column of cell in gridded plot
    min = 0
    max = 3
    step = 1

    [spp]
    description = Name of species

    [count]
    description = Number of individuals of a species in a cell

and tells us that the ANBO plant census was conducted on a 4m x 4m grid where each cell was 1m x 1m.

To examine how the SAD changes across scale, we will look at the average SAD in each 1m x 1m, 2m x 2m cell, 2m x 4m, 4m x 2m, and 4m x 4m plot.::

    import macroeco as meco

    pat = meco.empirical.Patch("/Users/mqwilber/Desktop/ANBO.txt")

    # Get the empirical SAD in each 1m x 1m cell
    splits1 = "row:4; column:4"

    # Get the empirical SAD in each 2m x 1m cell
    splits2 = "row:2; column:4"

    # Get the empirical SAD in 4 2m x 2m cells (upper left , upper right, lower left, lower right)
    splits3 = "row:2; column:2"

    # Get the empirical SAD in left half and right half 4m x 2m cells
    splits4 = "row:1; column:2"

    # Get the SAD for the full plot
    splits5 = "row:1; column:1"

    all_splits = [splits1, splits2, splits3, splits4, splits5]

    # Store all the empirical SAD results
    results = []

    for split in all_splits:
        results.append(meco.empirical.sad(pat, 'spp_col:spp; count_col:count', splits=split))

The parameter `results` stores the empirical SAD results across scales. For example, `results[0]` is a list of length 16 that has the SAD for each cell in the plot.

    >>> len(results[0])
    16
    >>> results[0][0]
    ('row>=-0.5; row<0.5; column>=-0.5; column<0.5',       spp    y
     1    cabr    2
     3    chst    1
     5    cran    1
     6    crcr    3
     10  grass   42
     15   phdi    8
     16  plsp1    1
     17   pypo    3
     19   ticr  140
     23  unsp4    1)

Notice that the first item in the tuple contains a description of the cell from which the SAD was calculated: the cell with a row coordinate between -0.5 and 0.5 and a column coordinate between -0.5 and 0.5 (the upper left hand cell...IMAGE?). The minus coordinate is a result of accounting for the `step` parameter (i.e. minimum unit between spatial coordinates) specified in the metadata file.  In this case, the `step` parameter is 1 because the minimum distance between two spatial points is 1m.  The second item in the tuple is the SAD for that cell.

Now we fit the SAD to a zero-truncated negative binomial distribution where the shape parameter `k` can give us some insight into how the SAD changes across scale.  When `k` approaches 0 the SAD is close to logseries and when `k` is 1 the SAD follows a Broken Stick distribution. ::

    # Fit the SAD

    import numpy as np
    import matplotlib.pyplot as plt

    # Store the average ks
    avg_ks = []

    for tres, split_str in zip(results, all_splits):

        within_scale_ks = []

        for split in tres:

            within_scale_ks.append(meco.models.nbinom_ztrunc.fit_mle(split[1]['y'])[1])

        avg_ks.append(np.mean(within_scale_ks))

    # Plot the results
    areas = [1, 2, 4, 8, 16]
    plt.plot(areas, avg_ks, '-o')
    plt.xlabel("Area")
    plt.ylabel("k of zero-truncated NBD")

`k` is clearly decreasing with increasing scale.








