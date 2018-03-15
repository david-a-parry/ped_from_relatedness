try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup

setup(
    name = "ped_from_relatedness",
    packages = ["ped_from_relatedness"],
    version = "0.0.1",
    description = "Reconstruct or test small pedigress from relatedness data",
    author = "David A. Parry",
    author_email = "david.parry@igmm.ed.ac.uk",
    url = "https://github.com/gantzgraf/ped_from_relatedness",
    license='GPLv3',
    install_requires=[
          'pysam',
          'parse_vcf',
      ],
    scripts = ["bin/ped_from_relatedness"],
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
)
