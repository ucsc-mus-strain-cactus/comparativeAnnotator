#!/usr/bin/python
# A setup script to install the _twobit module
# ---
#

from distutils.core import setup
from distutils.extension import Extension

extensions = []
extensions.append(Extension("_twobit", ["_twobit.pyx"]))

def main():
    setup(name="twobit",
        ext_modules=extensions,
        cmdclass={'build_ext': build_ext})
      
if __name__ == "__main__":
    main()
