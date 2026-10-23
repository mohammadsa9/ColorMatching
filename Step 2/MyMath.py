import numpy as np
import math


def dot(m):
    result = m[0].dot(m[1])
    for i in range(2, len(m)):
        result = result.dot(m[i])
    return result


def sum(m):
    result = np.add(m[0], m[1])
    for i in range(2, len(m)):
        result = np.add(result, m[i])
    return result


def D1(m):
    return m.flatten()


def D2(m):
    return np.array([m])


def inv(m):
    return np.linalg.pinv(m)


def reverse(m):
    return np.linalg.pinv(m)


def array_distance(start, distance, end):
    result = []
    rg = (end-start)/distance
    for i in range(int(rg)+1):
        result.append(int(start + i*distance))
    return result


def array_repeat(value, repeat):
    result = []
    for i in range(int(repeat)):
        result.append(value)
    return result


def cleanNaN(m):
    return m[~np.isnan(m)]


def applyFunction(a, b):
    return np.vectorize(b)(a)
