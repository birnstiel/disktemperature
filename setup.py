"""
Setup file for the disk temperature package.
"""
from setuptools import setup, find_packages
import os

setup(
    name='disktemperature',
    use_scm_version=True,
    description='Simple protoplanetary disk temperature estimates and auxiliary functions',
    long_description=open(os.path.join(
        os.path.dirname(__file__), 'README.md')).read(),
    url='http://www.til-birnstiel.de',
    author='Til Birnstiel',
    author_email='birnstiel@me.com',
    packages=find_packages(),
    license='GPLv3',
    include_package_data=True,
    setup_requires=['setuptools_scm'],
    install_requires=[
        'astropy',
        'scipy',
        'numpy',
        'matplotlib',
        'tbhmie'],
    dependency_links=['git+ssh://git@github.com/birnstiel/tbhmie/tarball/master#egg=tbhmie-0.0.3'],
    zip_safe=False)
