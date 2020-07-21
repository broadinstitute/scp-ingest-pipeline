from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='scp-ingest-pipeline',
    version='1.3.10',
    description='ETL pipeline for single-cell RNA-seq data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/broadinstitute/scp-ingest-pipeline',
    author='Single Cell Portal team',
    author_email='scp-support@broadinstitute.zendesk.com',
    install_requires=[
        'google-cloud-storage',
        'google-cloud-bigquery',
        'requests',
        'numpy',
        'scipy',
        'jsonschema',
        'pandas',
        'pandocfilters',
        'colorama',
        'dataclasses',
        'mypy_extensions',
        'pymongo',
        'loompy',
        'backoff',
        'opencensus',
        'opencensus-context',
        'opencensus-ext-stackdriver',
        'google-cloud-trace',
        'grpcio',
    ],
    packages=find_packages(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: BSD License',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.7',
    ],
    python_requires='>=3.7',
)
