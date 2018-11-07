#!/usr/bin/env python

# Prefer setuptools over distutils
from setuptools import setup, find_packages

setup(
    name='variantFormatter',
    version='0.0.1',
    description='Accurate conversion of Pseudo VCF variants into the HGVS format and mapping between reference sequences',
    long_description=open('README.txt').read(),
    url='',
    author='Peter J. Causey-Freeman',
    author_email='pjf9@leicester.ac.uk',
	package_data={"variantFormatter": ["configuration/*.ini"],},
    packages=find_packages(),
    include_package_data=True,
    license="https://www.gnu.org/licenses/agpl-3.0.en.html",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Audience
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Specify the Python versions
        'Programming Language :: Python :: 2.7',
    ],
 
    # What does your project relate to?
      keywords=[
          "bioinformatics",
          "computational biology",
          "genome variants",
          "genome variation",
          "genomic variants",
          "genomic variation",
          "genomics",
          "hgvs",
          "HGVS",
          "sequencevariants",
      ],

	# List run-time dependencies here.  These will be installed by pip when the project is installed.
    install_requires=[
        "hgvs == 1.1.3", # This will install BioPython
		"biocommons.seqrepo == 0.4.2",
		"configparser == 3.5.0",
    ],
)