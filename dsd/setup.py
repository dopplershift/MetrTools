from distutils.core import setup
from dsd import __version__

setup(name = 'dsd',
      version = str(__version__),
      author = 'Ryan May',
      author_email = 'rmay31@gmail.com',
      platforms = ['Linux', 'UNIX', 'Windows', 'MacOSX'],
      description = 'Functions for doing drop-size distribution calculations',
      url = 'http://weather.ou.edu/~rmay/research.html',
      packages = ['dsd']
      )

