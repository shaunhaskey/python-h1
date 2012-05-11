from distutils.core import setup

execfile('h1/version.py')

setup(name='h1',
      version=__version__,
      packages=['h1']
      )
