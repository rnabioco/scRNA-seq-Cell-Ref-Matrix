from pymongo import MongoClient

client = MongoClient('mongo "mongodb+srv://cluster0-fqznl.mongodb.net/test" --username sid')

db = client.get_database('cyberhood')

user = db.user