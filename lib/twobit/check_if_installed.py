import subprocess

try:
    from Cython.Distutils import build_ext
except ImportError:
    print ("Please install Cython. (pip install Cython)")
    raise

try:
    import _twobit
except ImportError:
    print ("Installing twobit module")
    sys.exit(3)