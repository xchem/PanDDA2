

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

def serialize_residue_assignments(residue_assignments):

    unflattened = unflatten(residue_assignments)
    return unflattened

def unserialize_residue_assignments(residue_assignments):
    flattened_1 = flatten(residue_assignments,depth=3)
    
    return flattened_1


def read_residue_assignments(sequence_assignment_file):
    with open(sequence_assignment_file, 'r') as f:
        serialized_assignment_file = json.load(f)
    residue_assignments = unserialize_residue_assignments(serialized_assignment_file)
    return residue_assignments


def output_residue_assignments(residue_assignments, sequence_assignment_file):
    # Covert tuple keys to 
    serialized_msa = serialize_residue_assignments(residue_assignments)
    with open(sequence_assignment_file, 'w') as f:
        json.dump(serialized_msa, f)
