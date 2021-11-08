import os
from params import *
import math 

# calculates how many task ID's we need
def num_id():
    return len(t23l) * len(m31l)

def id_to_2did(idx):
    t23id = math.floor(idx / len(m31l))
    m31id = idx % len(m31l)
    return t23id, m31id

def id_to_param(idx):
    t23id, m31id = id_to_2did(idx)
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
    
def read_output():
    res = np.ndarray(shape = (len(m31l), len(t23l)))
    for filename in os.listdir(dir_name):
        if filename.endswith(".txt"):
            with open(os.path.join(dir_name, filename), "r") as file1:
                try:
                    temp = [line.rstrip('\n') for line in file1]
                    # print(temp)
                    # print(int(float(temp[0])))
                    # i = 0
                    # Lines = file1.readLines()
                    # for line in Lines:
                    #     temp[i] = float(line)
                    thidx, dmidx = id_to_2did(int(float(temp[0])))
                    res[dmidx][thidx] = float(temp[1])
                except:
                    print(filename)
    return res 
