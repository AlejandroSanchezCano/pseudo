
'''
import numpy as np
import h5py

with h5py.File("/home/asanche2/asanche2/ppi/data/per-residue.h5", "r") as file:
    print(f"number of entries: {len(file.items())}")
    for sequence_id, embedding in file.items():
        print(
            f"  id: {sequence_id}, "
            f"  embeddings shape: {embedding.shape}, "
            f"  embeddings mean: {np.array(embedding).mean()}"
        )
'''
from structure import Structure
from timeit import timeit
import utils
import msgpack
import pickle

# PICKLE
pic = "utils.unpickling('/home/asanche2/asanche2/ppi/data/Structures/A0A0A0KC22.pkl')"
time_pickle = timeit(pic, number = 10, globals=globals())

# MSG
s = Structure('A0A6J1JL93')
s.__dict__['structure'] = 'structure'

m = msgpack.packb(s.__dict__, use_bin_type=True)
with open('aa', 'w') as handle:
    pickle.dump(
                obj = str(m),
                file = handle, 
                protocol = pickle.HIGHEST_PROTOCOL
                )
    
def msss():
    m = pickle.load('aa')
    msgpack.unpackb(m, raw=False)
msg = 'msss()'
time_msg = timeit(msg, number = 10, globals=globals())

print(time_pickle, time_msg)
