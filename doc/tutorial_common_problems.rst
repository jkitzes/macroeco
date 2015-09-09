.. _common-problems:

================================
Common problems (and solutions!)
================================

The majority of errors that we have encountered when using `macroeco` and MacroecoDesktop have to do with improperly formatted strings passed to the various arguments to the `empirical` functions (:doc:`empirical`). See below for some common errors and examples.

If you are running MacroecoDesktop and something is not working, check the log file `_log.txt` that is saved in the `results` folder. This will contain information on the analysis. If you see a message similar to the following ::

    Possible formatting error(s) in parameters.txt: parameters ['divs'] on lines [15]

you likely have an improperly formatted `cols`, `divs`, or `splits` parameter in your parameter file. Double check the parameter on the given line to make sure it is properly formatted.

Common Errors
=============

Here are some errors that result from improperly formatted strings and how to fix them.

1. `IndexError: list index out of range`

This is probably the most common error you will get when using `macroeco`.  It probably occurred because `cols` string is not properly formatted.  For example, check out the following code

    >>> import macroeco as meco
    >>> pat = meco.empirical.Patch("ANBO.txt")
    >>> sad = meco.empirical.sad(pat, cols="spp_col;spp; count_col:count", splits="")
    IndexError: list index out of range

Notice that the `cols` string `"spp_col;spp; count_col:count"` is not properly formatted. This string should contain a semi-colon separated list with colons with colons word pairs. In this case, change the first `;` to `:` and things will work

    >>> sad = meco.empirical.sad(pat, cols="spp_col:spp; count_col:count", splits="")

You might have also gotten this error if you were using the `empirical.sar` function and have a typo in the `divs` string. For example,

    >>> sar = meco.empirical.sar(pat, cols="spp_col:spp; count_col:count; x_col:row; y_col:column", splits="", divs="1,1; 1:2")

    IndexError: list index out of range

The `divs` list should read `"1,1; 1,2"` (switching out the colon for a comma).

Note that if you are using MacroecoDesktop, both of these errors will output a warning message to your `_log.txt` file as described above.


2. `ValueError: need more than 1 value to unpack`

This is another common error you may see and it is exactly the same as the one described above, but in this case there is likely a problem with your `splits` string.  For example, check out the following code

    >>> import macroeco as meco
    >>> pat = meco.empirical.Patch("ANBO.txt")
    >>> sad = meco.empirical.sad(pat, cols="spp_col:spp; count_col:count", splits="year:split; column,4")
    ValueError: need more than 1 value to unpack

Notice that the `splits` string is not formatted properly.  This should string should contain colon pairs separated by semicolons. Notice the `,` after `column`.  Make the following change and it will work

    >>> sad = meco.empirical.sad(pat, cols="spp_col:spp; count_col:count", splits="year:split; column:4")

Note that if you are using MacroecoDesktop, both of these errors will output a warning message to your `_log.txt` file as described above.

3. `KeyError`

This probably happened because one of the column names in you strings is misnamed.  For example,

    >>>  sad = meco.empirical.sad(pat, cols="spp_col:spp; count_col:coun", splits="year:split; column:4")
    KeyError: 'coun'

In this case, I misspelled `count` as `coun` and `macroeco` could not find it.

More subtle errors: Name collision
==================================








