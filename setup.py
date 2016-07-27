#! /usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from pkg_resources import parse_requirements, RequirementParseError

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('LICENSE') as f:
    license_str = f.read()

try:
    with open('requirements.txt') as f:
        ireqs = parse_requirements(f.read())
except RequirementParseError:
    raise
requirements = [str(req) for req in ireqs]


setup(name='TCGAIntegrator',
      version='0.1.0',
      description='Automation of Firebrowse API for creating integrated TCGA datasets',
      author='Lee Cooper',
      author_email='lee.cooper@emory.edu',
      url='https://github.com/cooperlab/TCGAIntegrator',
      packages=['tcgaintegrator'],
      package_dir={'tcgaintegrator': 'tcgaintegrator'},
      include_package_data=True,
      install_requires=requirements,
      license=license_str,
      zip_safe=False,
      keywords='tcgaintegrator',
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'License :: OSI Approved :: Apache Software License',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2',
          'Topic :: Scientific/Engineering :: Artificial Intelligence',
          'Topic :: Software Development :: Libraries :: Python Modules',
      ],
)
