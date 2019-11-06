import os
from pymongo import MongoClient


def get_mongo_db():
    user = os.environ['MONGO_USER']
    password = os.environ['MONGO_PASSWORD']
    host = os.environ['MONGO_HOST']
    scp_env = os.environ['SCP_ENV']

    db_name = 'single_cell_portal_' + scp_env

    client = MongoClient(
        host,
        username=user,
        password=password,
        authSource=db_name,
        authMechanism='SCRAM-SHA-1',
    )
    return client[db_name]


db = get_mongo_db()

genes = db.genes
gene = {'gene': 'HBB'}
gene_mongo_id = genes.insert_one(gene).inserted_id

print(gene_mongo_id)
