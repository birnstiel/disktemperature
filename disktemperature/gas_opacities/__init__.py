"""
Initializaition file for sub-package `gas_opacities`.
"""
from setuptools_scm import get_version as __get_version
__version__ = __get_version(root='../..', relative_to=__file__)

from . import semenov_opacity

__all__ = ['semenov_opacity']
