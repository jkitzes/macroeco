.. _recipes:

=======================
MacroecoDesktop Recipes
=======================

To provide a "jump start" on setting up analyses for MacroecoDesktop, the sample parameter file below contains a variety of runs that perform different types of calculations on the demo dataset provided with Macroeco. This file, or individual runs from this file (consisting of a run title in square brackets and all subsequent lines until the next run title), can be copied and pasted into parameters files and modified as needed.

The lines beginning with the ``#`` symbol are comments. They are purely for information and are ignored by MacroecoDesktop. In some cases below, lines containing variables are prefaced by the ``#`` symbol, indicating that they are "commented out" and will not affect the analysis. Removing the ``#`` at the start of these lines will have the effect described in the associated comment for that line. ::

    # The runs below provide examples of empirical data analysis, some with
    # model comparisons.

    # A simple species abundance distribution for the full plot
    [SAD]
    analysis = sad

    metadata = ANBO.txt

    models = logser_uptrunc; lognorm
    log_y = True  # Log transform the y axis of output plots

    # Four separate SAD's for the four quadrants of the plot
    # cols is only required if it is not set in the metadata file
    [SAD4]
    analysis = sad

    metadata = ANBO.txt
    #cols = spp_col:spp; count_col:count; x_col:row; y_col:column
    splits = row:2; column:2
    clean = True  # Remove species with 0 individuals from SADs

    models = logser_uptrunc; lognorm
    log_y = True  # Log transform the y axis of output plots

    # Empirical spatial abundance distribution for all 16 cells
    [SSAD]
    analysis = ssad

    metadata = ANBO.txt
    splits = row: 4; column: 4

    # Species area relationship
    [SAR ANBO]
    analysis = sar

    metadata = ANBO.txt
    divs = 1,1;1,2;2,1;2,2;2,4;4,4

    models = mete_sar_iterative
    #ear = True  # Endemics area relationship instead of species area
    log_y = True
    log_x = True

    # Gridded commonality, calculating Sorensen index for each pair of cells
    [Commonality]
    analysis = comm_grid

    metadata = ANBO.txt
    #subset = row>=2;column>=2  # Use only cells in rows 2-3 and columns 2-3
    cols = spp_col:spp; count_col:count; x_col:row; y_col:column
    #splits = row:2  # Perform analysis once for rows 0-1 and again for 2-3
    divs = 2,2
    #metric = Jaccard  # Use Jaccard instead of Sorensen index

    models = power_law

    # O ring measure of distance decay
    # This measure is best suited to point count census data
    [Oring]
    analysis = o_ring

    metadata = ANBO.txt
    cols = spp_col:spp; count_col:count; x_col:row; y_col:column
    spp = 'crcr'
    bin_edges = 0, 1, 2, 3, 4

    # The runs below provide examples of model exploration

    # pmf of geometric distribution
    [Geom-pmf]
    analysis = geom.pmf

    x = range(10)  # x values from 0 to 9
    p = 0.5

    # Shape parameter of upper truncated geometric distribution
    [GeomUptrunc-p]
    analysis = geom_uptrunc.translate_args

    mu = 5
    b = 20

    # Fit parameters of lognormal to a small data set
    [Lognorm-fit]
    analysis = lognorm.fit_mle

    data = 2,2,5,8,4,3


    # Draw random variates from a conditioned negative binomial distribution
    [Cnbinom-random]
    analysis = cnbinom.rvs

    mu = 10
    k_agg = 2
    b = 15
    size = 10
