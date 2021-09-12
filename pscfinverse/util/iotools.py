""" Module containing useful functions for string parsing """
import re
import pathlib

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

def writeCsvLine(fname, data, writeStyle='a'):
    styleList = ['a','w','x']
    if writeStyle not in styleList:
        raise(ValueError("writeCsvLine requires writeStyle in {}; gave {}".format(styleList,writeStyle)))
    strline = ','.join(map(str,data))
    strline += '\n'
    with open(fname,writeStyle) as f:
        f.write(strline)

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
        wordlist = [w for w in splitPattern.split(line) if w.strip()]
        for word in wordlist:
            yield word.strip()

class FileParser:
    """ Class to parse a file word-by-word with type-specific access. """
    
    def __init__(self, fname):
        """ Initialize FileParser.
        
        Parameters
        ----------
        fname : string, Pathlib.path
            The path to the file to be parsed.
        """
        self._fname = fname
        self._file_obj = open(fname, 'r')
        self._words = wordsGenerator(self._file_obj)
        self._terminated = False
        self._last_word = None
        self._next_word = next(self._words)
        self._count = 0
    
    def __enter__(self):
        return self
    
    def __exit__(self, *args, **kwargs):
        self._file_obj.close()
        self._terminated = True
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if self._terminated:
            raise(StopIteration)
        self._last_word = self._next_word
        try:
            self._next_word = next(self._words)
        except StopIteration:
            self._next_word = None
            self._terminated = True
        #TRACER.trace("Token '{}' read from file.".format(self._last_word),TraceLevel.ECHO)
        print("Token '{}' read from file.".format(self._last_word))
        return self._last_word
    
    def next_string(self):
        """ Increment iterator and return next value as string. """
        out = next(self)
        #TRACER.trace("Token '{}' returned as string.".format(out),TraceLevel.DETAIL)
        return out
    
    def next_int(self):
        """ Increment iterator, but return next value as int.
        
        Raises ValueError if next value is not interpretable as int.
        """
        nextStr = next(self)
        out = int(nextStr)
        #TRACER.trace("Token '{}' returned as int {}.".format(nextStr,out),TraceLevel.DETAIL)
        return out
    
    def next_float(self):
        """ Increment iterator, but return next value as float.
        
        Raises ValueError if next value is not interpretable as float.
        """
        nextStr = next(self)
        out = float(nextStr)
        #TRACER.trace("Token '{}' returned as float {}.".format(nextStr,out),TraceLevel.DETAIL)
        return out
    
    def next_number(self):
        """ Increment iterator, but return next value as either int or float.
        
        Raises ValueError if next value is not interpretable as a number
        """
        nextStr = next(self)
        out = str_to_num(nextStr)
        #TRACER.trace("Token '{}' returned as number {}.".format(nextStr,out),TraceLevel.DETAIL)
        return out
    
    def next_try_number(self):
        """ Return next value as number if possible, else return string. """
        nextstr = next(self)
        try:
            out = str_to_num(nextstr)
            #TRACER.trace("Token '{}' returned as number {}.".format(nextStr,out),TraceLevel.DETAIL)
        except ValueError:
            out = nextstr
            #TRACER.trace("Token '{}' returned as string.".format(out),TraceLevel.DETAIL)
        return out
    
    def __str__(self):
        formstr = "< FileParser of '{}' at {}:{} >"
        return formstr.format(self._fname, self._count, self._next_word)


