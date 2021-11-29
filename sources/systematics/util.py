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

def fisherid(t1, t2):
    return len(t1) * len(t2)

def fisher_id_to_2did(idx, t1, t2):
    t1id = math.floor(idx / len(t2))
    t2id = idx % len(t2)
    return t1id, t2id

def fisher_id_to_param(idx, t1, t2):
    t1id, t2id = fisher_id_to_2did(idx, t1, t2)
    t1val = t1[t1id]
    t2val = t2[t2id]
    return t1val, t2val 

def local():
    totidx = len(t23l) * len(m31l)
    print("index should be ", totidx)
    for i in range(totidx):
        file_name = "{}.txt".format(i)
        complete_name = os.path.join(dir_name, file_name)
        file1 = open(complete_name, "a")
        file1.close()
    
def fisher_local(t1, t2 = np.array([0])):
    totidx = len(t1) * len(t2)
    print("index should be ", totidx)
    for i in range(totidx):
        file_name = "{}.txt".format(i)
        complete_name = os.path.join(dir_name, file_name)
        file1 = open(complete_name, "a")
        file1.close()
    
def read_output(dir_name = dir_name):
    res = np.ndarray(shape = (len(m31l), len(t23l)))
    # print(m31l)
    # print(res.shape)
    for filename in os.listdir(dir_name):
        if filename.endswith(".txt"):
            print(filename)
            with open(os.path.join(dir_name, filename), "r") as file1:
                try:
                    temp = [line.rstrip('\n') for line in file1]
                    # print(temp)
                    # print(temp)
                    # print(int(float(temp[0])))
                    # i = 0
                    # Lines = file1.readLines()
                    # for line in Lines:
                    #     temp[i] = float(line)
                    thidx, dmidx = id_to_2did(int(float(temp[0])))
                    res[dmidx][thidx] = float(temp[1])
                except:
                    # print("excepted")
                    index = int(filename.split(".")[0])
                    # print(index)
                    thidx, dmidx = id_to_2did(int(float(index)))
                    res[dmidx][thidx] = 0
                    continue
    print(res)
    return res 

def process(chils, thels):
    processed_chi = []
    processed_the = []
    for i in range(len(chils)):
        if chils[i] != 0: # or i >= len(chils) / 2:
            processed_chi.append(chils[i])
            processed_the.append(thels[i])
    return processed_chi, processed_the