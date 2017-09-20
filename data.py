def load_data_by_type(file_name, dtype, classes):
    """Load the data in file name with the corresponding type from a list of
    rregistered classes.

    Parameters:
        file_name (str): data file name.
        dtype (str): data type.
        classes (dict): dictionary relating each class to an object.
    """
    if dtype not in classes:
        raise TypeError('Type %s does not exist' % dtype)
    return classes[dtype](file_name)
