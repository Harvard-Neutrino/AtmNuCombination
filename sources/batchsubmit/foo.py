import sys

def foo(x, y):
    return x ** 2 + y ** 2

dm = float(sys.argv[1])
th = float(sys.argv[2])

print(foo(dm, th))