""" Module containing useful functions for string parsing """


def str_to_num(stringNum):
    """
    Convert string to either int or float with priority on int.
    
    Examples:
        "10.000" returns float
        "5" returns int
        "abc" raises error
        "1e25" returns float
    
    Parameters
    ----------
    stringNum : str
        String representing a numeric value
    
    Returns
    -------
    int if string represents a valid integer, otherwise float
    
    Raises
    ------
    ValueError : if stringNum does not represent a numeric value
    TypeError : if stringNum is not string, bytes-like object, or number
    """
    try:
        return int(stringNum)
    except(ValueError):
        return float(stringNum)
