import numpy as np


def list_to_float(l):
    for i in range(0, len(l)):
        l[i] = float(l[i])


def list_to_int(l):
    for i in range(0, len(l)):
        l[i] = int(l[i])


def list_to_str(l):
    for i in range(0, len(l)):
        l[i] = str(l[i])


class Myobj2:
    def __init__(self, d):
        self.d = d


class Myobj:
    def __init__(self, c, d):
        self.a = 1
        self.b = [1, 2, (1, 2), "1, 2"]
        self.c = c
        self.d = Myobj2(d)


def create_l(mo):
    nl = np.zeros(9)
    for i in range(0, 9):
        mo.d = i * 10
        nl[i] = mo.d
    return nl


if __name__ == "__main__":

    print('part 1')
    ul = [1, 2, 3, 4, 5]
    print(ul)
    list_to_float(ul)
    print(ul)
    list_to_int(ul)
    print(ul)
    list_to_str(ul)
    print(ul)

    print('part 2')

    umo = Myobj(1, -10)
    ul = create_l(umo)
    print(ul)
    umo.d = -10
    print(ul)