#!/usr/bin/env python3

from distutils.core import setup
setup(name='brefito',
      version='1.0',
      py_modules=[],
      scripts=['scripts/brefito'],
      packages=['workflow'],
      package_dir={'workflow': 'workflow'},
      package_data={'workflow': ['envs/*.yml', 'rules/*.smk']},
      )