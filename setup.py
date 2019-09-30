from setuptools import setup, find_packages

setup(
    name='scp_ingest_pipeline',
    version='0.4',
    description='ETL pipeline for single-cell RNA-seq data',
    author='Single Cell Portal team',
    author_email='scp-support@broadinstitute.zendesk.com',
    install_requires=[
        'google-cloud-firestore',
        'google-cloud-storage',
        'requests',
        'numpy',
        'scipy',
        'jsonschema',
        'pandas',
        'pandocfilters',
    ],
    packages=find_packages(),
)
