"""
General purpose library of things.

Author: Ian Fiddes
With contributions from Dent Earl

"""


import os
import argparse


class FullPaths(argparse.Action):
        """
        Expand user- and relative-paths
        https://gist.github.com/brantfaircloth/1443543
        """
        def __call__(self, parser, namespace, values, option_string=None):
                setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def formatRatio(numerator, denominator):
        if denominator == 0:
                return float("nan")
        return float(numerator)/denominator


def DirType(d):
    """ given a string path to a directory, D, verify it can be used.
    """
    d = os.path.abspath(d)
    if not os.path.exists(d):
        raise argparse.ArgumentTypeError('DirType:%s does not exist' % d)
    if not os.path.isdir(d):
        raise argparse.ArgumentTypeError('DirType:%s is not a directory' % d)
    if os.access(d, os.R_OK):
        return d
    else:
        raise argparse.ArgumentTypeError('DirType:%s is not a readable dir' % d)


def FileType(f):
    """ given a string path to a file, F, verify it can be used.
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