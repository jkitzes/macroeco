
import numpy as np
import pandas as pd

def format_columnar():
    """
    """
    pass


def format_dense(data_path, non_spp_cols, delimiter=",", na_values="",
        item_col="spp", count_col="count", nan_to_zero=False, drop_na=False):
    """
    Formats dense data type to columnar data type.

    Takes in a dense data type and converts into a stacked data type.

    Parameters
    ----------
    data_path : str
        A path to the dense data
    non_spp_cols : list
        A list of columns in the data that are not species columns
    delimiter : str
        The delimiter for the dense data. Default, ","
    na_values : int, float, str, or list
        Values to be labeled as NA. Default, ""
    item_col : str
        Name of the item column in the formatted data. Default, "spp"
    count_col : str
        Name of the count column in the formatted data. Default, "count"
    nan_to_zero : bool
        Set all nans to zero. Default, False
    drop_na : bool
        Drop all columns with nan in the dataset. Default, False

    Notes
    -----
    Examples of Dense Data conversion


    """

    # Default arguments
    base_data = pd.read_csv(data_path, sep=delimiter,
                    na_values=na_values)

    # Stack data in columnar form
    indexed_data = base_data.set_index(keys=non_spp_cols)
    columnar_data = indexed_data.stack()
    columnar_data = columnar_data.reset_index()

    # Set nans to zero?
    if nan_to_zero:
        columnar_data[np.isnan(columnar_data)] = 0

    # Drop nans?
    if drop_na:
        columnar_data = columnar_data.dropna(how="any")

    # Rename columns
    num = len(non_spp_cols)
    columnar_data.rename(columns={0: count_col, 'level_%i' % num:
        item_col}, inplace=True)

    return columnar_data


def format_transect():
    """
    """
    pass

def format_grid():
    """
    """
    pass
