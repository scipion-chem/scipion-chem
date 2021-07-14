"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='scipion-chem',
    version='3.0.0',
    description='Base Scipion plugin defining objects and protocols for CHEMoinformatics',
    long_description=long_description,
    url='https://github.com/scipion-chem/scipion-chem',
    author='Carlos Oscar Sorzano',
    author_email='scipion@cnb.csic.es',
    keywords='scipion pdb scipion-3',
    packages=find_packages(),
    install_requires=[requirements],
    include_package_data=True,
    package_data={
       'pwchem': ['tool.png'],
    },
    entry_points={
        'pyworkflow.plugin': 'pwchem = pwchem'
    }
)
