import numpy as np
import math


def dot(m):
    result = m[0].dot(m[1])
    for i in range(2, len(m)):
        result = result.dot(m[i])
    return result


def zebra(m):
    m = D1(m)
    temp = []
    for i in range(len(m)):
        if i % 2 == 1:
            temp.append(m[i])
    return np.array(temp)


def zip(m, k):
    for i in range(k):
        m = zebra(m)
    return m


def array_zebra(arr):
    list = []
    for i in range(len(arr)):
        list.append(zebra(arr[i]))
    return np.array(list)


def PC(arr, k, vector):
    return dot([arr, vector[0:k].T])


def array_PC(arr, k, vector):
    list = []
    for i in range(len(arr)):
        list.append(PC(arr[i], k, vector))
    return np.array(list)


def array_zip(arr, k):
    for i in range(k):
        arr = array_zebra(arr)
    return arr


def D1(m):
    return m.flatten()


def D2(m):
    return np.array([m])


def inv(m):
    return np.linalg.pinv(m)


def sum(m):
    result = np.add(m[0], m[1])
    for i in range(2, len(m)):
        result = np.add(result, m[i])
    return result


def reverse(m):
    return np.linalg.pinv(m)


def array_distance(start, distance, end):
    result = []
    rg = (end - start) / distance
    for i in range(int(rg) + 1):
        result.append(int(start + i * distance))
    return result


def array_repeat(value, repeat):
    result = []
    for i in range(int(repeat)):
        result.append(value)
    return result


def cleanNaN(m):
    for i in range(len(m)):
        if math.isnan(m[i]):
            break
    return np.delete(m, np.s_[i:], 0)


def applyFunction(a, b):
    return np.vectorize(b)(a)
