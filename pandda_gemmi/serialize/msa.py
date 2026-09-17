from rich import print as rprint
import json

def unflatten(d):
    flat = {}
    for k, v in d.items():
        if type(k) is tuple:
            subkeys = []
            for j in range(len(k)):
                subkey = k[j]
                subdict = flat
                for key in subkeys:
                    subdict = subdict[key]

                if subkey not in subdict:
                    subdict[subkey] = {}

                subkeys.append(subkey)
                if len(k) == len(subkeys):
                    subdict[subkey] = v
        else:
            flat[k] = v

    return flat


def flatten(d, depth=3):
    if depth == 1:
        return d
    flattened = {}
    for k, v in d.items():
        if type(v) is dict:
            flattened_v = flatten(v, depth=depth-1)
            # print(flattened_v)
            for k_2, v_2 in flattened_v.items():
                if type(k_2) is tuple:
                    flattened[tuple([k,] + [_k for _k in k_2])] = v_2
                else:
                    flattened[tuple([k, k_2])] = v_2
        else:
            d[k] = v

    return flattened

def serialize_msa(msa):
    unflattened_inner = {}
    for k, v in msa.items():
        unflattened_inner[k] = unflatten(v)
    unflattened = unflatten(unflattened_inner)
    return unflattened

def unserialize_msa(msa):
    flattened_1 = flatten(msa,depth=2)
    
    flattened_2 = {}
    for k, v in flattened_1.items():
        flattened_2[k] = flatten(v, depth=2)

    return flattened_2
    ...

def output_msa(msa, msa_file):
    # Covert tuple keys to 
    serialized_msa = serialize_msa(msa)
    with open(msa_file, 'w') as f:
        json.dump(serialized_msa, f)

def read_msa(msa_file):
    with open(msa_file, 'r') as f:
        serialized_msa = json.load(f)
    msa = unserialize_msa(serialized_msa)
    return msa