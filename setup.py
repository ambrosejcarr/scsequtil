from setuptools import setup

setup(
    name='util',
    description='Utilities for large scale single cell data processing',
    author='Ambrose J. Carr',
    author_email='mail@ambrosejcarr.com',
    package_dir={'': 'src'},
    packages=['util', 'util/plot', 'util/test'],
    install_requires=[
        'numpy',
        'pysam',
        'matplotlib>=2.0.0',
        'nose2'
    ],
    include_package_data=True
)
