from setuptools import setup

setup(
    name='scsequtil',
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
    include_package_data=True
)
