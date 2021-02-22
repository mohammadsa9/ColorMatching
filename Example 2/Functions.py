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


class CIE:
    def __init__(obj):
        obj.X, obj.Y, obj.Z
        obj.A, obj.B, obj.L

    def getX(obj):
        return obj.X

    def getY(obj):
        return obj.Y

    def getZ(obj):
        return obj.Z

    def getA(obj):
        return obj.A

    def getB(obj):
        return obj.B

    def getL(obj):
        return obj.L

    def setX(obj, x):
        obj.X = x

    def setY(obj, y):
        obj.Y = y

    def setZ(obj, z):
        obj.Z = z

    def setA(obj, a):
        obj.A = a

    def setB(obj, b):
        obj.B = b

    def setL(obj, l):
        obj.L = l


class Observation:
    def __init__(object, LightSource, Viewer, R=1):
        object.LightSource = LightSource
        object.Viewer = Viewer
        object.E = mm.D1(LightSource.E)
        object.R = R
        object.xbar = mm.D1(Viewer.xbar)
        object.ybar = mm.D1(Viewer.ybar)
        object.zbar = mm.D1(Viewer.zbar)
        object.k = object.getK()
        if not np.isscalar(R):
            object.light = Observation(LightSource, Viewer, 1)
        else:
            object.light = object

    def getK(object):
        k = 100/(np.sum((object.E)*(object.ybar)))
        return k

    def getX(object):
        if np.isscalar(object.R):
            R = 1
        else:
            R = mm.D1(object.R)
        x = object.k * np.sum((object.E)*(object.xbar) * (R))
        return x

    def getY(object):
        if np.isscalar(object.R):
            R = 1
        else:
            R = mm.D1(object.R)
        y = object.k * np.sum((object.E)*(object.ybar) * (R))
        return y

    def getZ(object):
        if np.isscalar(object.R):
            R = 1
        else:
            R = mm.D1(object.R)
        z = object.k * np.sum((object.E)*(object.zbar) * (R))
        return z

    def getL(object):
        return 116 * object.f(object.getY()/object.light.getY()) - 16

    def getA(object):
        return 500 * (object.f(object.getX()/object.light.getX()) - object.f(object.getY()/object.light.getY()))

    def getB(object):
        return 200 * (object.f(object.getY()/object.light.getY()) - object.f(object.getZ()/object.light.getZ()))

    def f(obj, t):
        sigma = 6/29
        if t > pow(sigma, 3):
            return pow(t, 1/3)
        else:
            return (t/(3*pow(sigma, 2))) + (4/29)


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
        eps = np.finfo(float).eps
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

    def getResult(obj, target, source):

        A = obj.locate(target).T
        one = mm.array_repeat(1, len(target)+1)
        A = np.vstack((A, one))

        B = np.hstack((target, [1]))
        B = np.array([B]).T

        Variables = mm.inv(A).dot(B)
        result = Variables.T.dot(obj.getSource(source))

        return result


def pow2(x):
    return pow(x, 2)


def distance(x1, x2):
    A = np.array(x1[0])
    x2 = [x2]
    B = np.vstack((x2, x2, x2, x2))
    a = mm.sum([A, -1*B])
    b = mm.applyFunction(a, pow2)
    c = np.sum(b)
    return pow(c, 0.5)


def myInterploration(OBS, munsell_XYZ, munsell_R):
    light_source = OBS.LightSource
    viewer = OBS.Viewer
    R_sample = OBS.R[0]
    sample = Observation(light_source, viewer, R_sample)
    sample_XYZ = [sample.getX(), sample.getY(), sample.getZ()]
    XYZ_BEST = np.array(
        [munsell_XYZ[0], munsell_XYZ[0], munsell_XYZ[0], munsell_XYZ[0]])
    R_BEST = np.array(
        [munsell_R[0], munsell_R[0], munsell_R[0], munsell_R[0]])
    mindis = distance(
        [munsell_XYZ[0], munsell_XYZ[0], munsell_XYZ[0], munsell_XYZ[0]], sample_XYZ)
    max = 500
    for x in range(max):
        for y in range(x, max):
            for z in range(y, max):
                for d in range(z, max):
                    dis = distance(
                        [[munsell_XYZ[x], munsell_XYZ[y], munsell_XYZ[z], munsell_XYZ[d]]], sample_XYZ)
                    if dis <= mindis:
                        XYZ_BEST = np.array(
                            [munsell_XYZ[x], munsell_XYZ[y], munsell_XYZ[z], munsell_XYZ[d]])
                        R_BEST = np.array(
                            [munsell_R[x], munsell_R[y], munsell_R[z], munsell_R[d]])
        A = XYZ_BEST.T
        one = mm.array_repeat(1, 4)
        A = np.vstack((A, one))

        B = np.hstack((sample_XYZ, [1]))
        B = np.array([B]).T

        Variables = mm.inv(A).dot(B)
        result = Variables.T.dot(R_BEST)
        return result


def Interploration(calc, OBS, munsell_R):
    sample_XYZ = [OBS.getX(), OBS.getY(), OBS.getZ()]
    checker = calc.possible(sample_XYZ)
    R_calc = calc.getResult(sample_XYZ, munsell_R).T
    return R_calc, checker


def revInterploration(calc, OBS, munsell_XYZ):
    R_sample = OBS.R
    R_s = R_sample.T
    checker = calc.possible(R_s)
    XYZ_calc = calc.getResult(R_s, munsell_XYZ)
    return XYZ_calc[0], checker
