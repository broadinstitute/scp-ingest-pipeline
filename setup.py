from setuptools import setup, find_packages

setup(
    name='scp_ingest_pipeline',
    version='0.1',
    description='ETL pipeline for single-cell RNA-seq data',
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
        # 'loompy',
        'backoff',
	'opencensus',
    ],
    packages=find_packages(),
)
