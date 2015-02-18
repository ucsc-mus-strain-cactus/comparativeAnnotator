"""
General purpose library of things.

Author: Ian Fiddes
With contributions from Dent Earl

"""


import os, argparse, sys, gzip, operator


class FullPaths(argparse.Action):
    """
    Expand user- and relative-paths
    https://gist.github.com/brantfaircloth/1443543
    """
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def opener(filename):
    """opens files that may be gzip files"""
    f = open(filename, 'rb')
    if (f.read(2) == '\x1f\x8b'):
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.close()
        return open(filename, 'r')


def classesInModule(module):
    """
    http://stackoverflow.com/questions/5520580/how-do-you-get-all-classes-defined-in-a-module-but-not-imported
    """
    md = module.__dict__
    return [
        md[c] for c in md if (
            isinstance(md[c], type) and md[c].__module__ == module.__name__
        )
    ]


def formatRatio(numerator, denominator):
    if denominator == 0:
        return float("nan")
    return float(numerator)/denominator
    

def DirType(d):
    """
    given a string path to a directory, D, verify it can be used.
    """
    d = os.path.abspath(d)
    if not os.path.exists(d):
        os.mkdir(d)
    if not os.path.isdir(d):
        raise argparse.ArgumentTypeError('DirType:%s is not a directory' % d)
    if os.access(d, os.R_OK):
        return d
    else:
        raise argparse.ArgumentTypeError('DirType:%s is not a readable dir' % d)


def FileType(f):
    """
    given a string path to a file, F, verify it can be used.
    """
    f = os.path.abspath(f)
    if not os.path.exists(f):
        raise argparse.ArgumentTypeError('FileType:%s does not exist' % f)
    if not os.path.isfile(f):
        raise argparse.ArgumentTypeError('FileType:%s is not a regular file' % f)
    if os.access(f, os.R_OK):
        return f
    else:
        raise argparse.ArgumentTypeError('FileType:%s is not a readable file' % f)

def combineDicts(a, b, op=operator.add):
    """
    http://stackoverflow.com/questions/11011756/is-there-any-pythonic-way-to-combine-two-dicts-adding-values-for-keys-that-appe
    """
    return dict(a.items() + b.items() + [(k, op(a[k], b[k])) for k in b.viewkeys() & a.viewkeys()])