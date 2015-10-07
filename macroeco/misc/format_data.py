import numpy as np
import pandas as pd

def data_read_write(data_path_in, data_path_out, format_type, **kwargs):
    """
    General function to read, format, and write data.

    Parameters
    ----------
    data_path_in : str
        Path to the file that will be read
    data_path_out : str
        Path of the file that will be output
    format_type : str
        Either 'dense', 'grid', 'columnar', or 'transect'
    kwargs
        Specific keyword args for given data types. See Notes

    Notes
    -----

    'Dense Parameters'

    non_label_cols : str
        Comma separated list of non label columns. ex. "lat, long, tree"
    sep : str
        The delimiter for the dense data. Default, ","
    na_values : int, float, str
        Value to be labeled as NA. Default, ""

    See misc.format_dense() for additional keyword parameters
    """

    if format_type == "dense":

        # Set dense defaults
        kwargs = _set_dense_defaults_and_eval(kwargs)

        # Try to parse non label columns appropriately
        try:
            nlc = [nm.strip() for nm in kwargs['non_label_cols'].split(",")]
            kwargs.pop('non_label_cols', None)
        except KeyError:
            raise KeyError("'non_label_cols' is a required keyword dense data")

        # Read data with dense specific keywords
        arch_data = pd.read_csv(data_path_in, sep=kwargs['delimiter'],
                                na_values=kwargs['na_values'])

        form_data = format_dense(arch_data, nlc, **kwargs)

    elif format_type == "grid":
        pass
    elif format_type == "stacked":
        pass
    elif format_type == "transect":
        pass
    else:
        raise NameError("%s is not a supported data format" % format_type)

    form_data.to_csv(data_path_out, index=False)


def format_dense(base_data, non_label_cols, **kwargs):
    """
    Formats dense data type to stacked data type.

    Takes in a dense data type and converts into a stacked data type.

    Parameters
    ----------
    data : DataFrame
        The dense data
    non_label_cols : list
        A list of columns in the data that are not label columns
    label_col : str
        Name of the label column in the formatted data. Default, "label"
    count_col : str
        Name of the count column in the formatted data. Default, "count"
    nan_to_zero : bool
        Set all nans to zero. Default, False
    drop_na : bool
        Drop all columns with nan in the dataset. Default, False

    Returns
    -------
    : DataFrame
        A formatted DataFrame in the stacked format


    Notes
    -----
    Example of Dense Data conversion

    >>> import pandas as pd
    >>> dense_data = pd.DataFrame({'row' : [1,2,1,2], 'column' : [1,1,2,2],
        'labelA': [1,0,3,4], 'labelB' : [3,2,1,4]})

    >>> dense_data
           column  labelA  labelB  row
    0       1       1       3       1
    1       1       0       2       2
    2       2       3       1       1
    3       2       4       4       2

    [4 rows x 4 columns]

    # labelA and labelB might be species names. 'row' and 'column'
    # are non-species names so pass these in as non_label_cols

    >>> stacked_data = format_dense(dense_data, ['row', 'column'])
    >>> stacked_data
       row  column   label  count
    0    1       1  labelA      1
    1    1       1  labelB      3
    2    2       1  labelA      0
    3    2       1  labelB      2
    4    1       2  labelA      3
    5    1       2  labelB      1
    6    2       2  labelA      4
    7    2       2  labelB      4

    [8 rows x 4 columns]


    """

    kwargs = _set_dense_defaults_and_eval(kwargs)

    # Stack data in columnar form.
    indexed_data = base_data.set_index(keys=non_label_cols)
    columnar_data = indexed_data.stack(dropna=False)
    columnar_data = columnar_data.reset_index()

    # Rename columns
    num = len(non_label_cols)
    columnar_data.rename(columns={0: kwargs['count_col'], 'level_%i' % num:
        kwargs['label_col']}, inplace=True)

    # Set nans to zero?
    if kwargs['nan_to_zero']:
        ind = np.isnan(columnar_data[kwargs['count_col']])
        columnar_data.loc[ind, kwargs['count_col']] = 0
        columnar_data.reset_index(inplace=True, drop=True)

    # Drop nans?
    if kwargs['drop_na']:
        columnar_data = columnar_data.dropna(how="any")
        columnar_data.reset_index(inplace=True, drop=True)

    return columnar_data


def _set_dense_defaults_and_eval(kwargs):
    """
    Sets default values in kwargs if kwargs are not already given.

    Evaluates all values using eval

    Parameters
    -----------
    kwargs : dict
        Dictionary of dense specific keyword args

    Returns
    -------
    : dict
        Default, evaluated dictionary

    """

    kwargs['delimiter'] = kwargs.get('delimiter', ',')
    kwargs['na_values'] = kwargs.get('na_values', '')
    kwargs['nan_to_zero'] = kwargs.get('nan_to_zero', False)
    kwargs['drop_na'] = kwargs.get('drop_na', False)
    kwargs['label_col'] = kwargs.get('label_col', 'label')
    kwargs['count_col'] = kwargs.get('count_col', 'count')

    for key, val in kwargs.iteritems():
        try:
            kwargs[key] = eval(val)
        except:
            kwargs[key] = val

    return kwargs

def format_stacked():
    """
    """
    pass

def format_transect():
    """
    """
    pass

def format_grid():
    """
    """
    pass
