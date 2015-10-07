.. _common-problems:

================================
Common problems (and solutions!)
================================

The majority of errors that we have encountered when using `macroeco` and MacroecoDesktop have to do with improperly formatted strings given as arguments to the `empirical` functions (:doc:`empirical`). See below for some common errors and examples.

If you are running MacroecoDesktop and something is not working, check the log file `_log.txt` that is saved in the `results` folder. This will contain information on the analysis. If you see a message similar to the following ::

    Possible formatting error(s) in parameters.txt: parameters ['divs'] on lines [15]

you likely have an improperly formatted `cols`, `divs`, or `splits` parameter in your parameter file. Double check the parameter on the given line to make sure it is properly formatted.

    .. note::

        **What is proper formatting?** Check out :ref:`using-macroeco` and :ref:`using-macroecodesktop` for some examples and :doc:`empirical` for method documentation. Here are some quick rules for `subset`, `cols`, `splits`, and `divs`.

        1. All parameters are semi-colon separated strings in the form `"item1; item2; item3"`

        2. For `cols` and `splits` each "item" is a colon separated pair: `cols="spp_col:spp; count_col:count"` or `splits="year:split; row:2"`

        3. For `splits` each item is comparison operator: `splits = "row>=3; year==2010"`

        4. For `divs` each item is a comma separated pair: `divs = "1, 1"` or `divs = "1, 1; 2, 2; 2, 1"`


Common Errors
=============

Here are some errors that result from improperly formatted strings and how to fix them.

1. `IndexError: list index out of range`

This is probably the most common error you will get when using `macroeco`.  It probably occurred because `cols` string is not properly formatted.  For example, check out the following code

    >>> import macroeco as meco
    >>> pat = meco.empirical.Patch("ANBO.txt")
    >>> sad = meco.empirical.sad(pat, cols="spp_col;spp; count_col:count", splits="")
    IndexError: list index out of range

The `cols` string `"spp_col;spp; count_col:count"` is not properly formatted. This string should contain a semi-colon separated list with colon separated pairs. In this case, change the first `;` to `:` and things will work

    >>> sad = meco.empirical.sad(pat, cols="spp_col:spp; count_col:count", splits="")

You might have also gotten this error if you were using the `empirical.sar` function and have a typo in the `divs` string. For example,

    >>> sar = meco.empirical.sar(pat, cols="spp_col:spp; count_col:count; x_col:row; y_col:column", splits="", divs="1,1; 1:2")
    IndexError: list index out of range

The `divs` list should read `"1,1; 1,2"` (switching out the colon for a comma).  Similarly, if you only want to calculate the species richness at one area, the `divs` string should read `"1,1"`.

.. note:: If you are using MacroecoDesktop, both of these errors will output a warning message to your `_log.txt` file as described above.


2. `ValueError: need more than 1 value to unpack`

This is another common error you may see and it is exactly the same as the one described above, but in this case there is likely a problem with your `splits` string.  For example, check out the following code

    >>> import macroeco as meco
    >>> pat = meco.empirical.Patch("ANBO.txt")
    >>> sad = meco.empirical.sad(pat, cols="spp_col:spp; count_col:count", splits="year:split; column,4")
    ValueError: need more than 1 value to unpack

Notice that the `splits` string is not formatted properly.  This should string should contain colon pairs separated by semicolons. Notice the `,` after `column`.  Make the following change and it will work

    >>> sad = meco.empirical.sad(pat, cols="spp_col:spp; count_col:count", splits="year:split; column:4")

If you are using MacroecoDesktop, both of this error will output a warning message to your `_log.txt` file as described above.

3. `KeyError`

This probably happened because one of the column names in the strings is misnamed.  For example,

    >>>  sad = meco.empirical.sad(pat, cols="spp_col:spp; count_col:coun", splits="year:split; column:4")
    KeyError: 'coun'

In this case, we misspelled `count` as `coun` and `macroeco` could not find it. This could happen in your `cols` parameters or your `splits` parameter.

Similarly, if you misname a column in `subset` you will also see a `KeyError`. For example

    >>> import macroeco as meco
    >>> pat = meco.empirical.Patch("ANBO.txt", subset="rw>2; year==2010")
    KeyError: "Column 'rw' not found"

`rw` needs to be changed to `row` (because that is what it is called in ANBO.txt)

.. warning:: If you are running MacroecoDesktop and your analysis never "completes" and you don't see any error displayed, check the log file or the MacroecoDesktop output widow for a message similar to `Possible formatting error(s) in parameters.txt: parameters ['divs'] on lines [15]`









