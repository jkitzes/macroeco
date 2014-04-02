import numpy as np
import pandas as pd

def format_columnar():
    """
    """
    pass


def format_dense(data_path, non_label_cols, evaluate=False, **kwargs):
    """
    Formats dense data type to columnar data type.

    Takes in a dense data type and converts into a stacked data type.

    Parameters
    ----------
    data_path : str
        A path to the dense data
    non_label_cols : list
        A list of columns in the data that are not label columns
    evaluate : bool
        If True, eval values in kwargs
    delimiter : str
        The delimiter for the dense data. Default, ","
    na_values : int, float, str
        Value to be labeled as NA. Default, ""
    item_col : str
        Name of the item column in the formatted data. Default, "label"
    count_col : str
        Name of the count column in the formatted data. Default, "count"
    nan_to_zero : bool
        Set all nans to zero. Default, False
    drop_na : bool
        Drop all columns with nan in the dataset. Default, False

    Notes
    -----
    Examples of Dense Data conversion...TODO


    """
    kwargs = set_defaults_and_eval(kwargs, evaluate)

    base_data = pd.read_csv(data_path, sep=kwargs['delimiter'],
                    na_values=kwargs['na_values'])

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
        columnar_data[kwargs['count_col']][ind] = 0

    # Drop nans?
    if kwargs['drop_na']:
        columnar_data = columnar_data.dropna(how="any")

    return columnar_data


def set_defaults_and_eval(kwargs, evaluate):
    """
    Sets default values in kwargs if kwargs are not already given
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

def format_transect():
    """
    """
    pass

def format_grid():
    """
    """
    pass
