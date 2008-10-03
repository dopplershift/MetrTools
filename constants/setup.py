from distutils.core import setup
from constants import __version__

setup(name = 'MetConstants',
      version = str(__version__),
      py_modules = ['constants'],
      author = 'Ryan May',
      author_email = 'rmay31@gmail.com',
      platforms = ['Linux', 'UNIX', 'Windows', 'MacOSX'],
      description = 'Collection of meteorological constants',
      url = 'http://weather.ou.edu/~rmay/research.html',
      )

