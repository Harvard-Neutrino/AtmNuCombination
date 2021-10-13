import os

for file in os.listdir('.'):
    if file.endswith('out'):
        x = file.split('_')
        index = file.split('_')[-1][0]
        print("output index is ", index)
        with open(file) as f:
            for line in f:
                print(line)
