import numpy as np
from scipy.spatial import Delaunay
import MyMath as mm


def find_KOVERS(r):
    return (pow((1 - r), 2)) / (2*r)


def find_r(k_on_s):
    return (1 + k_on_s) - pow((pow(k_on_s, 2) + 2*k_on_s), 0.5)


def dif(r):
    return (-2*pow(r, 2))/(1-pow(r, 2))


class Dye:

    def __init__(object, sample_num, size):
        object.num = sample_num
        object.size = size

    def setC(object, c):
        object.c = c

    def setSub(obj, sub):
        obj.sub = sub

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
                OBS.append(find_KOVERS(
                    obj.ref[y][x]) - find_KOVERS(obj.sub[x]))
            OBS = np.array(OBS)
            KOVERS = np.linalg.pinv(
                np.matmul(KSCOEF.transpose(), KSCOEF)).dot(KSCOEF.transpose()).dot(OBS)
            Res.append(KOVERS[0])
        return np.array(Res)


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
        k = 100/(np.sum((object.E)*(object.ybar)))
        return k

    def getX(object):
        if np.isscalar(object.R):
            R = 1
        else:
            R = object.R
        x = object.k * np.sum((object.E)*(object.xbar) * (R))
        return x

    def getY(object):
        if np.isscalar(object.R):
            R = 1
        else:
            R = object.R
        y = object.k * np.sum((object.E)*(object.ybar) * (R))
        return y

    def getZ(object):
        if np.isscalar(object.R):
            R = 1
        else:
            R = object.R
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


class Mixture:
    def __init__(object, r_sub):
        object.r_sub = r_sub
        object.result = mm.applyFunction(r_sub, find_KOVERS)

    def add(obj, c, KOVERS):
        obj.result = mm.sum([obj.result, c*KOVERS])

    def getKOVERS(obj):
        return obj.result

    def getR(obj):
        return mm.applyFunction(obj.result, find_r)

    def clear(obj):
        obj.result = mm.applyFunction(obj.r_sub, find_KOVERS)


class MyDelaunay:
    def __init__(obj, points, opt=''):
        obj.points = points
        obj.tri = Delaunay(points, furthest_site=False,
                           incremental=False, qhull_options=opt)

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
    k_est = Mix.getKOVERS()
    R_First = Mix.getR()
    """
    alpha1 = np.subtract(r_std, R_First)
    alpha2 = np.subtract(k_std, k_est)
    D = np.diagflat(alpha1/alpha2)
    """
    D = mm.applyFunction(r_std, dif)
    D = np.diagflat(D)
    CC = np.array([C_First[0], C_First[1], C_First[2]])
    pi = all_KOVERS
    T = np.hstack((viewer.xbar, viewer.ybar,
                   viewer.zbar)).transpose()
    E_Diag = np.diagflat(light_source.E)
    r_new = []
    ESTN = Observation(light_source, viewer, R_First)
    STD = Observation(light_source, viewer, r_std)
    for i in range(1000):
        #print("E", i, delta_E(ESTN, STD))
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


"""

def find_delta_t(r1, r2):
    return np.array([[r1.getX()-r2.getX()], [r1.getY()-r2.getY()], [r1.getZ()-r2.getZ()]])


def delta_E(obj1, obj2):
    temp = pow(obj1.getL()-obj2.getL(), 2) + pow(obj1.getA() -
                                                 obj2.getA(), 2) + pow(obj1.getB()-obj2.getB(), 2)
    return pow(temp, 0.5)


def RMS(r1, r2, mode_1="R", mode_2="R"):
    if mode_1 != "R":
        r1 = np.vectorize(find_KOVERS)(r1)
    if mode_2 != "R":
        r2 = np.vectorize(find_KOVERS)(r2)

    result = np.subtract(r1, r2)
    result = np.power(result, 2)
    result = np.sum(result)
    result = np.power(result, 0.5)
    return result




def inv_matrix(m):
    a, b = m.shape
    if a != b:
        raise ValueError("Only square matrices are invertible.")
    i = np.eye(a, a)
    return np.linalg.lstsq(m, i)[0]
    
class Spectrum:
    def __init__(object, r_est, LightSource):
        object.r_est = r_est
        object.E = LightSource.E
        object.xbar = LightSource.xbar
        object.ybar = LightSource.ybar
        object.zbar = LightSource.zbar
        object.k = LightSource.getK()
        object.Source = LightSource

    def getX(object):
        x = object.k * \
            np.sum(object.E.dot(object.xbar.transpose()).dot(object.r_est))
        return x

    def getY(object):
        y = object.k * \
            np.sum(object.E.dot(object.ybar.transpose()).dot(object.r_est))
        return y

    def getZ(object):
        z = object.k * \
            np.sum(object.E.dot(object.zbar.transpose()).dot(object.r_est))
        return z

    def getL(object):
        return 116 * pow(object.getY()/object.Source.getY(), 1/3) - 16

    def getA(object):
        return 500 * (pow(object.getX()/object.Source.getX(), 1/3) - pow(object.getY()/object.Source.getY(), 1/3))

    def getB(object):
        return 200 * (pow(object.getY()/object.Source.getY(), 1/3) - pow(object.getZ()/object.Source.getZ(), 1/3))


class getKOVERS:

    def __init__(object, r1, r2, r3, r4, r5, r6, r7, sub, c, size):
        object.r1 = r1
        object.r2 = r2
        object.r3 = r3
        object.r4 = r4
        object.r5 = r5
        object.r6 = r6
        object.r7 = r7
        object.sub = sub
        object.c = c
        object.size = size

    def setSub(object, sub):
        object.sub = sub

    def setC(object, c):
        object.c = c

    def get_r(obj):
        KSCOEF = obj.c
        Res = []
        for x in range(obj.size):
            OBS = np.array([find_KOVERS(obj.r1[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r2[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r3[x]) - find_KOVERS(obj.sub[x]),
                            find_KOVERS(obj.r4[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r5[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r6[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r7[x]) - find_KOVERS(obj.sub[x])])
            KOVERS = np.linalg.pinv(
                np.matmul(KSCOEF.transpose(), KSCOEF)).dot(KSCOEF.transpose()).dot(OBS)
            Res.append(find_r(KOVERS[0]))
        return np.array(Res)

    def get_k_on_s(obj):
        KSCOEF = obj.c
        Res = []
        for x in range(obj.size):
            OBS = np.array([find_KOVERS(obj.r1[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r2[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r3[x]) - find_KOVERS(obj.sub[x]),
                            find_KOVERS(obj.r4[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r5[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r6[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r7[x]) - find_KOVERS(obj.sub[x])])
            KOVERS = np.linalg.pinv(
                np.matmul(KSCOEF.transpose(), KSCOEF)).dot(KSCOEF.transpose()).dot(OBS)
            Res.append(KOVERS[0])
        return np.array(Res)


class getC:

    def __init__(object, r1, r2, r3, r_std, size):
        object.r1 = r1
        object.r2 = r2
        object.r3 = r3
        object.r_std = r_std
        object.size = size

"""
