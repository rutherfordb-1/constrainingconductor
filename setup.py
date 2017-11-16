""""constrainingConducttor: A package for organizing gromacs pulling simulations 

"""

from setuptools import setup
from setuptools.command.test import test as TestCommand
import sys

try:
    import mdtraj
except ImportError:
    print('Building and running msibi requires mdtraj. See '
          'http://mdtraj.org/latest/installation.html for help!')
    sys.exit(1)



setup(name='constrainingconductor',
      version='0.1',
      description='',
      url='https://github.com/ahy3nz/constrainingconductor',
      author='Alexander Yang',
      author_email='alexander.h.yang@vanderbilt.edu',
      #license='MIT',
      packages=['constrainingconductor'],
      #install_requires=requirements,
      zip_safe=False
      #test_suite='tests',
      #cmdclass={'test': PyTest},
      #extras_require={'utils': ['pytest']},
)
