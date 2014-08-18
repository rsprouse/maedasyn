Currently the C code that goes with this package is not available in the repository.
If you have access to it, copy *.c and *.h to a folder named 'c' in the same directory
as setup.py. After that, do:

  python setup.py install

This will compile the package with Cython and install it to your local python.


The 'maedasyn' script is a small example of how to use the package. Note that if you
run the script from the package's root directory python might get confused and not
find your installed package. Just run the script from another directory.
