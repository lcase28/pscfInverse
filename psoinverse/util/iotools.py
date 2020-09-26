""" Module containing useful functions for string parsing """
import re
import pathlib

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

def checkPath(root):
    """
    Checks if path exists, creates it if not.
    Checks if path is a directory, if not takes the
    parent directory.
    Fully resolves path.
    
    Parameters
    ----------
    root : pathlib.Path or string
        The directory path to be resolved.
    
    Returns
    -------
    resolvedRoot : pathlib.Path
        The resolved path to the specified directory.
    flag : bool
        True if the path resolved without error.
        False if path resolution threw one of:
        FileNotFoundError, FileExistsError, RuntimeError,
        which capture the expected throws from the Path
        object during resolution.
    """
    if not isinstance(root, pathlib.Path):
        pathlib.Path(root)
    root = root.resolve()
    if root.exists() and root.is_dir():
        return root.resolve(), True
    # path either does not exist or is not directory
    try:
        if not root.exists():
            root.mkdir(parents=True)
        if not root.is_dir():
            root = root.parent
        return root.resolve(), True
    except (FileNotFoundError, FileExistsError, RuntimeError):
        return None, False
    
def wordsGenerator(stringIterable):
    """
    Split a string iterable into individual words.
    
    Words, here, defined as groups of characters separated by whitespace, 
    or groups of characters (including whitespace) fully enclosed by double-
    or single-quotes.
    
    example: words(["This is 'a test'", "'one word' oneword two word."])
        returns a generator producing: ["This", "is", "'a test'", "'one word'", "oneword", "two", "word."]
    
    Paramters
    ---------
    stringIterable : iterable of string objects
        The "stream" of strings to convert to words, such as a file.
    """ 
    lineStream = iter(stringIterable)
    splitPattern = re.compile("(\s+|\\\".*?\\\"|'.*?')")
    for line in lineStream:
        line = line.strip()
        #wordlist = [w for w in re.split("(\s+|\\\".*?\\\"|'.*?')", line) if w.strip()]
        wordlist = [w for w in splitPattern.split(line) if w.strip()]
        for word in wordlist:
            yield word.strip()

def writeCsvLine(fname, data, writeStyle='a'):
    styleList = ['a','w','x']
    if writeStyle not in styleList:
        raise(ValueError("writeCsvLine requires writeStyle in {}; gave {}".format(styleList,writeStyle)))
    strline = ','.join(map(str,data))
    strline += '\n'
    with open(fname,writeStyle) as f:
        f.write(strline)

