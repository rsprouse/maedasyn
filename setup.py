#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

# The numpy.get_include() is often needed for compiling cython modules, but
# synth.pyx doesn't seem to need it, and the include_dirs line caused an
# import error in which
#
#   import maedasyn.synth
#
# resulted in the error
#
#   ImportError: No module named synth
#
# There was a workaround for this situation. This worked:
#
#   import numpy as np
#   import maedasyn.synth
#
# Removing the include_dirs line makes the workaround unnecessary and makes
# it possible to do import maedasyn.synth without importing numpy first.
# 
# If we ever need to revert to using include_dirs in order for compilation
# to succeed, we'll want to recall the workaround.
ext_modules = [
  Extension(
    name="maedasyn.synth",
    sources=["maedasyn/synth.pyx"],
    libraries = ["m"],
#    include_dirs=[numpy.get_include()],
    language="c",
  )
]

setup(
  name = 'maedasyn',
  cmdclass = {'build_ext': build_ext},
  #ext_package='klatt_wrap',
  ext_modules = ext_modules,
  packages = ['maedasyn'],
  scripts = [
    'scripts/maedasyn_wx'
  ],
  classifiers = [
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',
    'Topic :: Multimedia :: Sound/Audio :: Speech'
  ]
)
