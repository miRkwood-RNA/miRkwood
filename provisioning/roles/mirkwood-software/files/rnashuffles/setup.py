#!/usr/bin/python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

try:
    import shuffles
    version = shuffles.__version__
except ImportError:
    version = 'Undefined'

classifiers = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Intended Audience :: Developers',
    'Operating System :: Unix',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
]

packages = ['shuffles']
requires = ['pexpect', 'argparse']
entry_points = {
        'console_scripts': [
            'RNAshuffles = shuffles.Randfold:main',
            ]
        }


setup(name='shuffles',
      version=version,
      author='Bonsai Bioinformatics',
      author_email='bonsai-software@univ-lille1.fr',
      url='http://bioinfo.lifl.fr',
      classifiers=classifiers,
      packages=packages,
      install_requires=requires,
      entry_points=entry_points,
      )
