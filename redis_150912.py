import redis
redis = redis.Redis(host='localhost', port=6379, db=0) 
# Data Type : String Value
redis.set("name", "zedo") 
print(redis.get("name"))

# Data Type : Integer Value
redis.set("counter", 1) 
print(redis.get("counter")) # 1
redis.incr("counter") 
print(redis.get("counter")) # 2
redis.decr("counter")
print(redis.get("counter")) #1

# List : possible to duplicate values
redis.rpush("members", "r1") 
redis.rpush("members", "r2")
redis.lpush("members", "l1")
redis.lpush("members", "l2")
print(redis.lrange("members", 0, 0))
print(redis.lrange("members", 0, 1))
print(redis.lrange("members", 0, 2))
print(redis.llen("members"))
print(redis.lrange("members",0, redis.llen("members")-1))
print(redis.lindex("members",3))
print(redis.rpop("members"))
print(redis.lpop("members"))
print(redis.llen("members"))
print(redis.lrange("members",0, redis.llen("members")-1))
redis.delete("members") 

#Sets  : impossible to duplicate values
redis.sadd("members", "s1")
redis.sadd("members", "s1")
redis.sadd("members", "s2")
redis.sadd("members", "s3")
redis.sadd("members", "s4")
redis.sadd("members", "s5")
print(redis.smembers("members"))
redis.delete("members")

# Using JSON Format
import json
import redis
import codecs

r = redis.StrictRedis(host='localhost', port=6379, db=0)
txn= [
    {
        "FACTORY_CODE" : "P7",
        "SHOP_CODE" : "TFT",
        "LOT_ID": "P7T003F00",
        "CST_ID": "AST001",
        "TRANSACTION_INFO": {
            "SERIAL_NO": "201509121703010302",
            "USER_ID":"617836"
        }
    },
    {
        "FACTORY_CODE" : "P7",
        "SHOP_CODE" : "TFT",
        "LOT_ID": "P7T003F01",
        "CST_ID": "AST001",
        "TRANSACTION_INFO": {
            "SERIAL_NO": "201509121703010342",
            "USER_ID":"617836"
        }
    },
    {
        "FACTORY_CODE" : "P7",
        "SHOP_CODE" : "TFT",
        "LOT_ID": "P7T003F03",
        "CST_ID": "AST001",
        "TRANSACTION_INFO": {
            "SERIAL_NO": "201509121703010442",
            "USER_ID":"617836"
        }
    }
]


json_txn = json.dumps(txn)
r.set('txn', json_txn)
unpacked_txn = json.loads(r.get('txn').decode('utf-8-sig'))
txn == unpacked_txn

# Using JSON Format
import json
import redis
import codecs

r = redis.StrictRedis(host='localhost', port=6379, db=0)

txn1= {
        "FACTORY_CODE" : "P7",
        "SHOP_CODE" : "TFT",
        "LOT_ID": "P7T003F00",
        "CST_ID": "AST001",
        "TRANSACTION_INFO": {
            "SERIAL_NO": "201509121703010302",
            "USER_ID":"617836"
        }
    }
txn2 = {
        "FACTORY_CODE" : "P7",
        "SHOP_CODE" : "TFT",
        "LOT_ID": "P7T003F01",
        "CST_ID": "AST001",
        "TRANSACTION_INFO": {
            "SERIAL_NO": "201509121703010342",
            "USER_ID":"617836"
        }
    }
txn3 = {
        "FACTORY_CODE" : "P7",
        "SHOP_CODE" : "TFT",
        "LOT_ID": "P7T003F03",
        "CST_ID": "AST001",
        "TRANSACTION_INFO": {
            "SERIAL_NO": "201509121703010442",
            "USER_ID":"617836"
        }
    }

r.set('v1', json.dumps(txn1))
unpacked_txn = json.loads(r.get('v1').decode('utf-8-sig'))
r.lpush('v', unpacked_txn)

r.set('v2', json.dumps(txn2))
unpacked_txn = json.loads(r.get('v2').decode('utf-8-sig'))
r.lpush('v', unpacked_txn)

r.set('v3', json.dumps(txn3))
unpacked_txn = json.loads(r.get('v3').decode('utf-8-sig'))
r.lpush('v', unpacked_txn)
print(r.llen("v"))
print(r.lrange("v",0, r.llen("v")-1))
#r.delete('v')
