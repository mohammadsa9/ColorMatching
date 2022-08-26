import numpy as np
from scipy.spatial import Delaunay
from . import MyMath as mm
from pyhull.delaunay import DelaunayTri
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import LinearNDInterpolator


def find_KOVERS(r):
    return (pow((1 - r), 2)) / (2 * r)


def find_r(k_on_s):
    return (1 + k_on_s) - pow((pow(k_on_s, 2) + 2 * k_on_s), 0.5)


def dif(r):
    return (-2 * pow(r, 2)) / (1 - pow(r, 2))


class Dye:
    def __init__(object, sample_num, size):
        object.num = sample_num
        object.size = size
        object.KoverS = []

    def setC(object, c):
        object.c = mm.D2(mm.D1(np.array(c))).T

    def setSub(obj, sub):
        obj.sub = mm.D2(mm.D1(np.array(sub))).T

    def setR(obj, R):
        obj.ref = []
        for i in range(obj.num):
            obj.ref.append(R[i])

    def getKOVERS(obj):
        KSCOEF = obj.c
        Res = []
        for x in range(obj.size):
            OBS = []
            for y in range(obj.num):
                OBS.append(find_KOVERS(obj.ref[y][x]) - find_KOVERS(obj.sub[x]))
            OBS = np.array(OBS)
            KOVERS = (
                np.linalg.pinv(np.matmul(KSCOEF.transpose(), KSCOEF))
                .dot(KSCOEF.transpose())
                .dot(OBS)
            )
            Res.append(KOVERS[0])
            obj.KoverS = np.array(Res)
        return obj.KoverS

    def getR(obj):
        return mm.applyFunction(obj.KoverS, find_r)


class Viewer:
    def __init__(object, xbar, ybar, zbar):
        object.xbar = mm.D2(mm.D1(np.array(xbar))).T
        object.ybar = mm.D2(mm.D1(np.array(ybar))).T
        object.zbar = mm.D2(mm.D1(np.array(zbar))).T


class LightSource:
    def __init__(object, E):
        object.E = mm.D2(mm.D1(np.array(E))).T


class Observation:
    def __init__(object, LightSource, Viewer, R=1, name="", color="black"):
        object.LightSource = LightSource
        object.Viewer = Viewer
        object.E = LightSource.E
        object.R = R
        object.name = name
        object.color = color
        object.xbar = Viewer.xbar
        object.ybar = Viewer.ybar
        object.zbar = Viewer.zbar
        object.k = object.getK()
        if not np.isscalar(R):
            object.light = Observation(LightSource, Viewer, 1)
            # Fix R
            object.R = mm.D2(mm.D1(np.array(object.R))).T

        else:
            object.light = object

    def getK(object):
        k = 100 / (np.sum((object.E) * (object.ybar)))
        return k

    def getX(object):
        if np.isscalar(object.R):
            R = 1
        else:
            R = object.R
        x = object.k * np.sum((object.E) * (object.xbar) * (R))
        return x

    def getY(object):
        if np.isscalar(object.R):
            R = 1
        else:
            R = object.R
        y = object.k * np.sum((object.E) * (object.ybar) * (R))
        return y

    def getZ(object):
        if np.isscalar(object.R):
            R = 1
        else:
            R = object.R
        z = object.k * np.sum((object.E) * (object.zbar) * (R))
        return z

    def getxy(object):
        X = object.getX()
        Y = object.getY()
        Z = object.getZ()
        return [X / (X + Y + Z), Y / (X + Y + Z)]

    def getXYZ(object):
        X = object.getX()
        Y = object.getY()
        Z = object.getZ()
        return np.array([X, Y, Z])

    def getL(object):
        return 116 * object.f(object.getY() / object.light.getY()) - 16

    def getA(object):
        return 500 * (
            object.f(object.getX() / object.light.getX())
            - object.f(object.getY() / object.light.getY())
        )

    def getB(object):
        return 200 * (
            object.f(object.getY() / object.light.getY())
            - object.f(object.getZ() / object.light.getZ())
        )

    def getAB(object):
        A = object.getA()
        B = object.getB()
        return [A, B]

    def f(obj, t):
        sigma = 6 / 29
        if t > pow(sigma, 3):
            return pow(t, 1 / 3)
        else:
            return (t / (3 * pow(sigma, 2))) + (4 / 29)


class Compare:
    def __init__(object, obj1, obj2):
        object.mat1 = obj1
        object.mat2 = obj2
        object.mat1.R = mm.D1(object.mat1.R)
        object.mat2.R = mm.D1(object.mat2.R)

    def delta_E(obj):
        temp = (
            pow(obj.mat1.getL() - obj.mat2.getL(), 2)
            + pow(obj.mat1.getA() - obj.mat2.getA(), 2)
            + pow(obj.mat1.getB() - obj.mat2.getB(), 2)
        )
        return pow(temp, 0.5)

    def RMS(obj):
        result = np.subtract(obj.mat1.R, obj.mat2.R)
        result = np.power(result, 2)
        result = np.sum(result)
        result = np.power(result, 0.5)
        return result

    def delta_T(obj):
        return np.array(
            [
                [obj.mat1.getX() - obj.mat2.getX()],
                [obj.mat1.getY() - obj.mat2.getY()],
                [obj.mat1.getZ() - obj.mat2.getZ()],
            ]
        )

    def GFC(object):
        R = mm.D1(object.mat1.R)
        Rhat = mm.D1(object.mat2.R)
        GFC = np.inner(R.T, Rhat.T) / (
            ((sum(R.T ** 2)) ** (0.5)) * (sum(Rhat.T ** 2)) ** (0.5)
        )
        return GFC


class Mixture:
    def __init__(object, r_sub):
        object.r_sub = mm.D2(mm.D1(np.array(r_sub))).T
        object.result = mm.applyFunction(r_sub, find_KOVERS)

    def add(obj, c, KOVERS):
        KOVERS = mm.D2(mm.D1(np.array(KOVERS))).T
        obj.result = mm.sum([obj.result, c * KOVERS])

    def getKOVERS(obj):
        return obj.result

    def getR(obj):
        return mm.applyFunction(obj.result, find_r)

    def clear(obj):
        obj.result = mm.applyFunction(obj.r_sub, find_KOVERS)


class MyDelaunay:
    def __init__(obj, points, source):
        obj.Delaunay = LinearNDInterpolator(points, source, -1)
        obj.Nearest = NearestNDInterpolator(points, source, -1)

    def getInterpolation(obj, target, SP_Compare):
        result = obj.Delaunay(target)
        if result[0][0] == -1:
            raise Exception("out of gamut")
        return result

    def getBestInterpolation(obj, target, SP_Compare):
        result1 = obj.Delaunay(target)
        if result1[0][0] == -1:
            raise Exception("out of gamut")
        result2 = obj.Nearest(target)
        if SP_Compare.RMS(result2) <= SP_Compare.RMS(result1):
            return result2
        return result1

    def getNearest(obj, target, SP_Compare):
        result = obj.Nearest(target)
        if result[0][0] == -1:
            raise Exception("out of gamut")
        return result


class BestDelaunay:
    def __init__(obj, points, sources):
        obj.Delaunay = []
        for i in range(len(points)):
            obj.Delaunay.append(MyDelaunay(points[i], sources[i]))

    def getBestInterpolation(obj, target, SP_Compare):
        best = []
        best_RMS = 100
        for i in range(len(obj.Delaunay)):
            try:
                result = obj.Delaunay[i].getBestInterpolation(target, SP_Compare)
                if best_RMS > SP_Compare.RMS(result):
                    best_RMS = SP_Compare.RMS(result)
                    best = result
            except Exception as e:
                continue
        if best_RMS == 100:
            raise Exception("out of gamut")
        return best

    def getInterpolation(obj, target, SP_Compare):
        best = []
        best_RMS = 100
        for i in range(len(obj.Delaunay)):
            try:
                result = obj.Delaunay[i].getInterpolation(target, SP_Compare)
                if best_RMS > SP_Compare.RMS(result):
                    best_RMS = SP_Compare.RMS(result)
                    best = result
            except Exception as e:
                continue
        if best_RMS == 100:
            raise Exception("out of gamut")
        return best


class SpecialCompare:
    def __init__(obj, blue, red, yellow, Mix):
        obj.blue_KOVERS = blue
        obj.red_KOVERS = red
        obj.yellow_KOVERS = yellow
        obj.Mix = Mix
        obj.mat1 = []

    def setR1(obj, R1):
        obj.mat1 = mm.D1(R1)

    def RMS(obj, arr):
        obj.Mix.clear()
        obj.Mix.add(arr[0][0], obj.blue_KOVERS)
        obj.Mix.add(arr[0][1], obj.red_KOVERS)
        obj.Mix.add(arr[0][2], obj.yellow_KOVERS)
        obj.mat2 = mm.D1(obj.Mix.getR())
        result = np.subtract(obj.mat1, obj.mat2)
        result = np.power(result, 2)
        result = np.sum(result)
        result = np.power(result, 0.5)
        return result


class MyDelaunay2:
    def __init__(obj, points, opt=""):
        obj.points = points
        obj.tri = Delaunay(
            points, furthest_site=False, incremental=False, qhull_options=opt
        )

    def possible(obj, point):
        eps = np.finfo(float).eps
        # 0.009
        obj.s = obj.tri.find_simplex(point, bruteforce=False, tol=0)
        if obj.s == -1:
            # return True
            raise Exception("Out of gamut")
            return False
        else:
            return True

    def locate(obj, point):
        obj.possible(point)
        temp_point = obj.points[obj.tri.simplices][obj.s]
        return temp_point

    def getSource(obj, points):
        RR = np.array(points)
        return RR[obj.tri.simplices][obj.s]

    def getDetail(obj, point, points):
        s = obj.tri.find_simplex(point, bruteforce=False)
        RR = np.array(points)
        return RR[obj.tri.simplices][s]

    def getVariables(obj, target):
        A = obj.locate(target).T
        one = mm.array_repeat(1, len(target) + 1)
        A = np.vstack((A, one))

        B = np.hstack((target, [1]))
        B = np.array([B]).T

        Variables = mm.inv(A).dot(B)

        return Variables

    def getResult(obj, target, source):
        A = obj.locate(target).T
        one = mm.array_repeat(1, len(target) + 1)
        A = np.vstack((A, one))

        B = np.hstack((target, [1]))
        B = np.array([B]).T

        Variables = mm.inv(A).dot(B)
        result = Variables.T.dot(obj.getSource(source))
        # print("R1", result)
        # Another way

        # why = LinearNDInterpolator(obj.tri, source)
        # why = NearestNDInterpolator(obj.tri, source)
        # result = why(target)

        # print("R2", result)
        return result


"""
End Of File
"""

# Not needed anymore:


def findC1(all_KOVERS, delta_KOVERS):
    temp = mm.dot([all_KOVERS.transpose(), all_KOVERS])
    temp = mm.reverse(temp)
    c_all = mm.dot([temp, all_KOVERS.transpose(), delta_KOVERS])
    return c_all


def findC2(STD, r_sub, C_First, all_KOVERS, maxRMS):
    r_std = STD.R
    viewer = STD.Viewer
    light_source = STD.LightSource

    all_E = []
    blue_KOVERS = np.array([all_KOVERS.T[0]]).T
    red_KOVERS = np.array([all_KOVERS.T[1]]).T
    yellow_KOVERS = np.array([all_KOVERS.T[2]]).T

    k_std = mm.applyFunction(r_std, find_KOVERS)

    Mix = Mixture(r_sub)
    Mix.add(C_First[0], blue_KOVERS)
    Mix.add(C_First[1], red_KOVERS)
    Mix.add(C_First[2], yellow_KOVERS)
    R_First = Mix.getR()
    D = mm.applyFunction(r_std, dif)
    D = np.diagflat(D)
    CC = np.array([C_First[0], C_First[1], C_First[2]])
    pi = all_KOVERS
    T = np.hstack((viewer.xbar, viewer.ybar, viewer.zbar)).transpose()
    E_Diag = np.diagflat(light_source.E)
    r_new = []
    ESTN = Observation(light_source, viewer, R_First)
    STD = Observation(light_source, viewer, r_std)
    for i in range(1000):
        comp = Compare(STD, ESTN)
        D_E = comp.delta_E()
        if D_E <= maxRMS:
            break
        all_E.append(D_E)
        delta_t = comp.delta_T()
        TEDPI = mm.dot([T, E_Diag, D, pi])
        delta_c = mm.dot([mm.reverse(TEDPI), delta_t])
        CC = np.add(CC, delta_c)
        Mix.clear()
        Mix.add(CC[0], blue_KOVERS)
        Mix.add(CC[1], red_KOVERS)
        Mix.add(CC[2], yellow_KOVERS)
        k_new = Mix.getKOVERS()
        r_new = Mix.getR()
        ESTN = Observation(light_source, viewer, r_new)
    return [CC, all_E, i]