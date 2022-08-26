import numpy as np
import math


def dot(m):
    result = m[0].dot(m[1])
    for i in range(2, len(m)):
        result = result.dot(m[i])
    return result


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


def PC(arr, k, vector, R_mean):
    start = 0
    arr = D1(arr)
    arr = sum([arr, -1 * R_mean.T[0]])
    return dot([arr, vector[start : start + k].T])


def lift(arr):
    result = np.power(arr, 2)
    result = np.sum(result)
    result = np.power(result, 0.5)
    return np.hstack((arr, result))


def array_lift(arr):
    list = []
    for i in range(len(arr)):
        list.append(lift(arr[i]))
    return np.array(list)


def array_PC(arr, k, vector, R_mean):
    list = []
    for i in range(len(arr)):
        list.append(PC(arr[i], k, vector, R_mean))
    return np.array(list)


def array_zip(arr, k):
    for i in range(k):
        arr = array_zebra(arr)
    return arr


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


def RMS(arr1, arr2):
    arr1 = np.array(arr1)
    arr2 = np.array(arr2)
    arr1 = D1(arr1)
    arr2 = D1(arr2)
    result = np.subtract(arr1, arr2)
    result = np.power(result, 2)
    result = np.sum(result)
    result = np.power(result, 0.5)
    return result
