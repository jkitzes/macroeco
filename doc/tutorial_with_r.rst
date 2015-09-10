===========================
MacroecoDesktop for R users
===========================

Users who primarily work in R can access the functionality of Macroeco through the command line MacroecoDesktop interface.

First, install a working copy of MacroecoDesktop by following the installation instructions in :ref:`installation`. Windows and Linux users will need to install a Python environment and the ``macroeco`` package, while Mac users can instead install the standalone MacroecoDesktop program. Follow the :ref:`first-steps-macroeco-desktop` tutorial to create the "new_parameters.txt" file and ensure that your copy of MacroecoDesktop is working properly.

For all platforms and installation options, the basic idea will be to call MacroecoDesktop from an R script using the command line interface, wait for the analysis to complete, and then read in any output tables saved by MacroecoDesktop that will be used for further analysis.

As an example, the script below completes the following steps:

* Writes a "new_parameters.txt" file describing a desired MacroecoDesktop analysis
* Uses the ``system`` command within R to execute the MacroecoDesktop analysis specified in "new_parameters.txt"
* Reads in the resulting data tables
* Plots gridded distance decay data with a best fit power law curve
* Prints out the R2 value for the power law fit to the data::

    param_dir <- "~/Desktop/demo/"
    param_file <- "new_parameters.txt"

    cat("
    [DistanceDecay]

    analysis = comm_grid

    metadata = ANBO.txt
    cols = spp_col: spp; count_col: count; y_col: row; x_col: column
    divs = 4,4
    models = power_law
    ",file=paste(param_dir,param_file,sep=""), sep="\n")

    system(paste("mecodesktop ", param_dir, param_file, sep=""))

    data_models <- read.csv(paste(param_dir, "results/DistanceDecay/1_data_models.csv", sep=""))
    test_statistics <- read.csv(paste(param_dir, "results/DistanceDecay/1_test_statistics.csv", sep=""))

    plot(data_models$x, data_models$empirical)
    lines(data_models$x, data_models$power_law)

    test_statistics$R2

Mac users who installed the standalone MacroecoDesktop program should replace ``"mecodesktop "`` above with ``"/Applications/MacroecoDesktop.app/Contents/MacOS/mecodesktop "``.
