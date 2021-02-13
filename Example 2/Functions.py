import numpy as np
import MyMath as mm
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay


def find_KOVERS(r):
    return (pow((1 - r), 2)) / (2*r)


def find_r(k_on_s):
    return (1 + k_on_s) - pow((pow(k_on_s, 2) + 2*k_on_s), 0.5)


class Viewer:
    def __init__(object, xbar, ybar, zbar):
        object.xbar = xbar
        object.ybar = ybar
        object.zbar = zbar


class LightSource:
    def __init__(object, E):
        object.E = E


class Observation:
    def __init__(object, LightSource, Viewer, R=1):
        object.LightSource = LightSource
        object.Viewer = Viewer
        object.E = LightSource.E
        object.R = R
        object.xbar = Viewer.xbar
        object.ybar = Viewer.ybar
        object.zbar = Viewer.zbar
        object.k = object.getK()
        if not np.isscalar(R):
            object.light = Observation(LightSource, Viewer, 1)
        else:
            object.light = object

    def getK(object):
        k = 100/np.sum(object.E.dot(object.ybar.transpose()))
        return k

    def getX(object):
        if np.isscalar(object.R):
            x = object.k * np.sum(object.E.dot(object.xbar.transpose()))
        else:
            x = object.k * \
                np.sum(object.E.dot(object.xbar.transpose()).dot(object.R))
        return x

    def getY(object):
        if np.isscalar(object.R):
            y = object.k * np.sum(object.E.dot(object.ybar.transpose()))
        else:
            y = object.k * \
                np.sum(object.E.dot(object.ybar.transpose()).dot(object.R))
        return y

    def getZ(object):
        if np.isscalar(object.R):
            z = object.k * np.sum(object.E.dot(object.zbar.transpose()))
        else:
            z = object.k * \
                np.sum(object.E.dot(object.zbar.transpose()).dot(object.R))
        return z

    def getL(object):
        return 116 * pow(object.getY()/object.light.getY(), 1/3) - 16

    def getA(object):
        return 500 * (pow(object.getX()/object.light.getX(), 1/3) - pow(object.getY()/object.light.getY(), 1/3))

    def getB(object):
        return 200 * (pow(object.getY()/object.light.getY(), 1/3) - pow(object.getZ()/object.light.getZ(), 1/3))


class Compare:
    def __init__(object, obj1, obj2):
        object.mat1 = obj1
        object.mat2 = obj2

    def delta_E(obj):
        temp = pow(obj.mat1.getL()-obj.mat2.getL(), 2) + pow(obj.mat1.getA() -
                                                             obj.mat2.getA(), 2) + pow(obj.mat1.getB()-obj.mat2.getB(), 2)
        return pow(temp, 0.5)

    def RMS(obj):
        result = np.subtract(obj.mat1.R, obj.mat2.R)
        result = np.power(result, 2)
        result = np.sum(result)
        result = np.power(result, 0.5)
        return result

    def delta_T(obj):
        return np.array([[obj.mat1.getX()-obj.mat2.getX()], [obj.mat1.getY()-obj.mat2.getY()], [obj.mat1.getZ()-obj.mat2.getZ()]])


class MyDelaunay:
    def __init__(obj, points):
        obj.points = points
        obj.tri = Delaunay(points)

    def possible(obj, point):
        obj.s = obj.tri.find_simplex(point)
        if obj.s == -1:
            return False
        else:
            return True

    def locate(obj, point):
        obj.possible(point)
        temp_point = obj.points[obj.tri.simplices][obj.s]
        return temp_point

    def getSource(obj, munsell_R):
        RR = np.array(munsell_R)
        return RR[obj.tri.simplices][obj.s]


def Interploration(calc, OBS, munsell_R):
    light_source = OBS.LightSource
    viewer = OBS.Viewer
    R_sample = OBS.R

    sample = Observation(light_source, viewer, R_sample)
    sample_XYZ = [sample.getX(), sample.getY(), sample.getZ()]

    checker = calc.possible(sample_XYZ)

    A = calc.locate(sample_XYZ).T
    A = np.vstack((A, [1, 1, 1, 1]))

    sample_XYZ.append(1)
    B = np.array([sample_XYZ]).T

    Variables = mm.inv(A).dot(B)
    R_calc = Variables.T.dot(calc.getSource(munsell_R))
    R_calc = R_calc.T

    return R_calc, checker
