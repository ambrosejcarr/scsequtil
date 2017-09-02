from setuptools import setup

CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Natural Language :: English",
    "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: Implementation :: PyPy",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

setup(
    name='scsequtil',
    version='0.0.3',
    description='Utilities for large scale single cell data processing',
    author='Ambrose J. Carr',
    author_email='mail@ambrosejcarr.com',
    package_dir={'': 'src'},
    packages=['scsequtil', 'scsequtil/plot', 'scsequtil/test'],
    install_requires=[
        'numpy',
        'pysam',
        'matplotlib>=2.0.0',
        'nose2',
        'jinja2',
        'weasyprint'
    ],
    classifiers=CLASSIFIERS,
    include_package_data=True
)
