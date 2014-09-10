#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [
  Extension(
    name="maedasyn.synth",
    sources=["maedasyn/synth.pyx"],
    libraries = ["m"],
    include_dirs=[numpy.get_include()],
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
