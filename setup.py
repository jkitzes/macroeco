from setuptools import setup, find_packages

setup(
    name = 'macroeco',
    version = '1.0.1',
    packages = find_packages(),
    # entry_points = {'console_scripts': ['mecodesktop=macroeco:mecodesktop',],},
    package_data = {'': ['*.txt', '*.csv']},

    author = 'Justin Kitzes and Mark Wilber',
    author_email = 'jkitzes@berkeley.edu',
    description = 'Ecological pattern analysis in Python',
    long_description = open('README.rst').read(),
    license = 'BSD',
    keywords = ('ecology macroecology biology environment biodiversity '
                'informatics data science'),
    url = 'http://github.com/jkitzes/macroeco',

    classifiers = [
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
    ],

    install_requires = [
        'numpy>=1.6',
        'scipy>=0.12',
        'pandas>=0.14',
        'matplotlib>=1.3',
        'mpmath>=0.19',
        'configparser',
        'decorator',
        # 'shapely',  # Do not force install if user doesn't have
        # 'wxpython', 
    ],
)

# python setup.py sdist bdist_egg upload -r https://testpypi.python.org/pypi
