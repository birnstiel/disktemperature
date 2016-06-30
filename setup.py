from setuptools import setup
import os

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name='disktemperature',
    version='0.1',
    description='Simple protoplanetary disk temperature estimates and auxiliary functions',
    long_descrioption=read('README.md'),
    url='http://www.til-birnstiel.de',
    author='Til Birnstiel',
    author_email='birnstiel@me.com',
    license='GPLv3',
    packages=['disktemperature'],
    include_package_data=True,
    install_requires=[
        'astropy',
        'scipy',
        'numpy',
        'matplotlib'
        ],
    zip_safe=False)
