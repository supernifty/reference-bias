try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Methods for quantifying reference bias in resequencing techniques',
    'author': 'Peter Georgeson',
    'url': 'https://github.com/supernifty/reference-bias',
    'version': '0.1',
    'install_requires': ['nose'],
    'packages': ['bias'],
    'scripts': [ 'bin/calculate_bias.py' ],
    'name': 'reference-bias'
}

setup(**config)
