#!/usr/bin/env python3

from distutils.core import setup
setup(name='brefito',
      version='1.0',
      py_modules=[],
      scripts=['scripts/brefito'],
      packages=['brefito'],
      package_dir={'brefito': 'workflow'},
      package_data={'brefito': ['envs/*.yml', 'rules/*.smk']},
      )