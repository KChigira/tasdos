#!/usr/bin/env python
from setuptools import setup, find_packages
from tasdos.__init__ import __version__

setup(name='tasdos',
      version=__version__,
      description='Downstream analysis for targetted amplicon sequencing.',
      author='Koki Chigira',
      author_email='s211905s@st.go.tuat.ac.jp',
      url='https://github.com/KChigira/tasdos/',
      license='MIT',
      packages=find_packages(),
      install_requires=[
        'pandas',
        'matplotlib',
      ],
      entry_points={'console_scripts': [
            'tasdosA = tasdos.tasdosa:main',
            'tasdosB = tasdos.tasdosb:main',
            'tasdosC = tasdos.tasdosc:main',
            ]
      }
    )
