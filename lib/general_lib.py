"""
General purpose library of things.
"""
import os
import gzip
import operator
import itertools
import types
import errno
from collections import OrderedDict, Callable, namedtuple

__author__ = "Ian Fiddes"


class DefaultOrderedDict(OrderedDict):
    """
    Source: http://stackoverflow.com/a/6190500/562769
    """
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
           not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                               OrderedDict.__repr__(self))


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def skip_header(path):
    """
    The attributes file produced by the pipeline has a header. Skip it. Return a open file handle pointing to line 2.
    """
    f_h = open(path)
    _ = f_h.next()
    return f_h


def opener(filename):
    """
    Transparently opens a file that is gzipped or not
    """
    f = open(filename, 'rb')
    if f.read(2) == '\x1f\x8b':
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.close()
        return open(filename, 'r')


def classes_in_module(module):
    """
    http://stackoverflow.com/questions/5520580/how-do-you-get-all-classes-defined-in-a-module-but-not-imported
    """
    md = module.__dict__
    return [
        md[c] for c in md if (
            isinstance(md[c], type) and md[c].__module__ == module.__name__
        )
    ]


def functions_in_module(module):
    """
    http://stackoverflow.com/questions/5520580/how-do-you-get-all-classes-defined-in-a-module-but-not-imported
    """
    md = module.__dict__
    return [
        md[c] for c in md if (
            isinstance(md[c], types.FunctionType) and md[c].__module__ == module.__name__
        )
    ]


def format_ratio(numerator, denominator):
    """
    Convenience function that converts two numbers, integer or no, to a ratio
    """
    if denominator == 0:
        return float("nan")
    return float(numerator) / denominator


def combine_dicts(a, b, op=operator.add):
    """
    http://stackoverflow.com/questions/11011756/is-there-any-pythonic-way-to-combine-two-dicts-adding-values-for-keys-that-appe
    """
    return dict(a.items() + b.items() + [(k, op(a[k], b[k])) for k in b.viewkeys() & a.viewkeys()])


def tokenize_stream(stream):
    """
    Iterator through a tab delimited file, returning lines as list of tokens
    """
    for line in stream:
        if line != '' and not line.startswith("#"):
            tokens = line.rstrip().split("\t")
            yield tokens


def dict_to_named_tuple(d, name):
    return namedtuple(name, d.keys())(**d)


def merge_dicts(list_of_dicts):
    """
    This will merge a list of dicts. Any duplicate keys will end up with the last value seen.
    """
    return reduce(lambda a, d: a.update(d) or a, list_of_dicts, {})
