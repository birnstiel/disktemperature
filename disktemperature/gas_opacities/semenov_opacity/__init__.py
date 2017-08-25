"""
Initializaition file for sub-sub-package `semenov_opacities`.
"""
from setuptools_scm import get_version as __get_version
__version__ = __get_version(root='../../..', relative_to=__file__)

from . import semenov_opacs

__all__ = ['semenov_opacs']
