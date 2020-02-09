from setuptools import setup
import codecs
import os
import re


# parse version from init
# from: https://github.com/pypa/pip/blob/master/setup.py
here = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    with codecs.open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup(name='structure_factor_tools',
      version=find_version("structure_factor_tools", "__init__.py"),
      description='structure_factor_tools: Python library of tools to calculate the structure factor for simple persovkite systems and relate these to atomic site locations',
      long_description='structure_factor_tools: Python library of tools to calculate the structure factor for simple persovkite systems and relate these to atomic site locations',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Physics',
      ],
      keywords=[
          'microscopy',
          'STEM',
          'TEM',
          'fast pixelated detector',
          'fpd',
          'data storage',
          'EMD',
          'hdf5',
          'data analysis',
          'differential phase contrast',
          'DPC',
          'segmented DPC',
          'pixelated DPC',
          'virtual detector',
          'lattice analyis',
      ],
      url='',
      author='Tom Macgregor',
      author_email='t.macgregor.1@reasearch.gla.ac.uk',
      license='GPL v3',
      packages=['structure_factor_tools'],
      package_data={'structure_factor_tools': ['perceptually_uniform_cmap.npy']},
      install_requires=[
          'numpy',
          'scipy',
          'scikit-image',
          'matplotlib',
          'h5py',
          'tqdm',
          'hyperspy',
          'pandas',
      ],
      include_package_data=True,
      zip_safe=True)
