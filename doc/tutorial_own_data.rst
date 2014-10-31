.. _own-data:

==============
Preparing Data
==============

Both data tables and metadata must be provided to MacroecoDesktop and the package ``macroeco`` for empirical analyses. Data should be placed in a csv file following the basic structure described below. Metadata must also be prepared to describe features of the data table that cannot be inferred from the table itself (for example, the minimum and maximum values of the extent of a census, as these may be smaller and larger, respectively, than the minimum and maximum coordinates of recorded individuals).

.. note
   To avoid the possibility of errors, the names of the data table and metadata files should not contain any spaces. Additionally, the column headers within the data table must not contain any spaces.

Preparing Data Tables
---------------------

Data should be prepared as a csv (comma separated values) file. The first row should contain column names, and each subsequent row should refer to a single record, most commonly a combination of a species identifier and a coordinate location. For point census data, each record will identify a single individual, while gridded census data will generally have an additional "count" column that gives the number of individuals of a species found in a grid cell.

Other columns, such as those identify genera, plot ID, etc. may also be included. The ``splits`` argument used by the empirical data analysis functions can be easily used to divide the analysis according to the values found in any provided column. For example, splitting a species area analysis on a column containing a plot ID will perform a separate species area analysis for each unique plot.

The demo data file ANBO.csv provides an examples of a correctly formatted data table file.

Preparing a Metadata File
-------------------------

Macroeco requires a metadata file to be provided along with each data table file (both MacroecoDesktop and ``macroeco`` require the user to provide the location of a metadata file, not the data table itself). The metadata file contains basic descriptive information about the data table as well as parameter values that are necessary for empirical data analysis.

The format of the metadata file is very similar to the parameters files used to describe analyses for MacroecoDesktop. A metadata file has an initial section called Description, followed by a section containing information for each column in the data table.

The metadata file ANBO.txt is shown here. ::

    [Description]
    name = Anzo Borrego
    author = Mary Ellen Harte and John Harte
    description = Vegetation census conducted at Anza-Borrego Desert State Park. Site in Indian Valley at N 32' 52.091", W 116' 14.447". Elevation 1195 feet. Census was conducted on a 4 m x 4 m grid, with 16 grid cells each 1 m2 in area.
    citation = Unpublished

    datapath = ANBO.csv
    cols = spp_col:spp; count_col: count; x_col; row: y_col; column

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

The initial section, Description, begins with a number of variables providing basic information on the data table.

The ``datapath`` variable in this section gives the location of the data table file relative to the metadata file. If the data file and metadata file are in the same directory, as is usually the case, then datapath is simply the name of the data table file.

The ``cols`` variable here provides an opportunity to identify the columns in the data table that indicate different values used in empirical data analysis. The four special columns shown here, which are common to most data tables, are

* spp_col - the species identifier
* count_col - optional column with the number of individuals of a species at a point (if count_col is not given, each row is taken to represent a single individual)
* x_col - the x coordinate of the record location
* y_col - the y coordinate of the record location

The value of ``cols`` can also be specified separately in individual runs in MacroecoDesktop or when calling individual functions in ``macroeco``. The values given here in the metadata file are the defaults which are used if ``cols`` is not specified in these other locations.

The remaining sections each refer to a column in the data table. Each section begins with a short description of the data in that column. Additionally, numeric columns (any column that can be split or subset by a numeric value) must have a minimum and maximum value and a step size giving the precision of the census. These are most commonly used with coordinate columns where, for example, the min and max values give the extent of the census and the step gives the minimum distance between two individuals.

The demo metadata file ANBO.txt contains the metadata shown above.

Using Data Files with Macroeco
------------------------------

Once the data and metadata files are prepared, they can be used with both MacroecoDesktop and ``macroeco``.

In MacroecoDesktop, each run that involves empirical data analysis must contain the variable ``metadata_path``, which should indicate the path of the metadata file relative to the parameters file. If the parameters file and the data file are in the same folder, this is simply the name of the metadata file.

In ``macroeco``, the absolute path to the metadata file (or the relative path from the present working directory) is a required argument to the Patch class of the empirical subpackage. Patch objects are required for all empirical pattern analysis using the functions in empirical.

