import os
from params import *
import math 

# calculates how many task ID's we need
def num_id():
    return len(t23l) * len(m31l)


def id_to_param(idx):
    t23id = math.floor(idx / len(m31l))
    m31id = idx % len(m31l)
    t23val = t23l[t23id]
    m31val = m31l[m31id]
    return t23val, m31val 

def local():
    totidx = len(t23l) * len(m31l)
    print("index should be ", totidx)
    for i in range(totidx):
        file_name = "{}.txt".format(i)
        complete_name = os.path.join(dir_name, file_name)
        file1 = open(complete_name, "a")
        file1.close()