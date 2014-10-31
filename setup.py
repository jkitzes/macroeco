from setuptools import setup, find_packages
from macroeco import __version__

setup(
    name = 'macroeco',
    version = __version__,
    packages = find_packages(),
    entry_points = {'console_scripts': ['mecodesktop=macroeco:mecodesktop',],},
    package_data = {'': ['*.txt', '*.csv']},

    author = 'Justin Kitzes and Mark Wilber',
    author_email = 'jkitzes@berkeley.edu',
    description = 'Ecological pattern analysis in Python',
    long_description = open('README.rst').read(),
    license = 'BSD',
    keywords = ('ecology biology environment conservation biodiversity '
                'informatics data science'),
    url = 'http://github.com/jkitzes/macroeco',

    classifiers = [
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",],

    install_requires = [
        'numpy>=1.6',
        'scipy>=0.12',
        'pandas>=0.13',
        'matplotlib',
        # 'shapely',  # Do not force install if user doesn't have
        'configparser',
        'decorator',
        'twiggy'],
)

# python setup.py sdist bdist_egg upload -r https://testpypi.python.org/pypi