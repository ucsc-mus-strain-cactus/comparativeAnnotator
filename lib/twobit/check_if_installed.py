try:
    from Cython.Distutils import build_ext
except ImportError:
    print "Please install Cython. (pip install Cython)"
    raise

try:
    import _twobit
except ImportError:
    print "Not installed. Run 'python setup_twobit.py install' with sudo (or on a personal python installation)"
