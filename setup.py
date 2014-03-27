try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name = 'macroeco',
    version= '0.3',
    description = 'Analysis of ecological patterns in Python',
    author = 'Justin Kitzes, Mark Wilber, Chloe Lewis',
    url = 'https://github.com/jkitzes/macroeco',
    packages = ['macroeco'],
    license = 'BSD',
)
