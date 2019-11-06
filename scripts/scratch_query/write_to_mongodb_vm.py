import os
from pymongo import MongoClient


def get_mongo_client():
    user = os.environ['MONGO_USER']
    password = os.environ['MONGO_PASSWORD']
    host = os.environ['MONGO_HOST']
    # scp_env = os.environ['SCP_ENV']

    client = MongoClient(
        host,
        username=user,
        password=password,
        authSource='single_cell_portal_development',
        authMechanism='SCRAM-SHA-1',
    )
    return client


db = get_mongo_client()['single_cell_portal_development']

genes = db.genes
gene = {'gene': 'HBB'}
gene_mongo_id = genes.insert_one(gene).inserted_id

print(gene_mongo_id)
