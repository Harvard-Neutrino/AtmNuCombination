from params import *
import math 

# calculates how many task ID's we need
def num_id():
    return len(t23l) * len(m31l)


def id_to_param(idx):
    t23id = math.floor(idx / 5)
    m31id = idx % 5
    t23val = t23l[t23id]
    m31val = m31l[m31id]
    return t23val, m31val 