import numpy as np
import MyMath as mm


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
            # OBS = np.array([find_KOVERS(obj.r1[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r2[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r3[x]) - find_KOVERS(obj.sub[x]),
            #                find_KOVERS(obj.r4[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r5[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r6[x]) - find_KOVERS(obj.sub[x]), find_KOVERS(obj.r7[x]) - find_KOVERS(obj.sub[x])])
            OBS = np.array(OBS)
            KOVERS = np.linalg.pinv(
                np.matmul(KSCOEF.transpose(), KSCOEF)).dot(KSCOEF.transpose()).dot(OBS)
            Res.append(KOVERS[0])
        return np.array(Res)


def find_KOVERS(r):
    return (pow((1 - r), 2)) / (2*r)


def find_r(k_on_s):
    return (1 + k_on_s) - pow((pow(k_on_s, 2) + 2*k_on_s), 0.5)


def inv_matrix(m):
    a, b = m.shape
    if a != b:
        raise ValueError("Only square matrices are invertible.")
    i = np.eye(a, a)
    return np.linalg.lstsq(m, i)[0]


class Viewer:
    def __init__(object, xbar, ybar, zbar):
        object.xbar = xbar
        object.ybar = ybar
        object.zbar = zbar


class LightSource:
    def __init__(object, E, Viewer):
        object.E = E
        object.xbar = Viewer.xbar
        object.ybar = Viewer.ybar
        object.zbar = Viewer.zbar

    def getK(object):
        k = 100/np.sum(object.E.dot(object.ybar.transpose()))
        return k

    def getX(object):
        k = object.getK()
        x = k * np.sum(object.E.dot(object.xbar.transpose()))
        return x

    def getY(object):
        k = object.getK()
        y = k * np.sum(object.E.dot(object.ybar.transpose()))
        return y

    def getZ(object):
        k = object.getK()
        z = k * np.sum(object.E.dot(object.zbar.transpose()))
        return z


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


def find_delta_t(r1, r2):
    return np.array([[r1.getX()-r2.getX()], [r1.getY()-r2.getY()], [r1.getZ()-r2.getZ()]])


def findC1(all_KOVERS, delta_KOVERS):
    temp = mm.dot([all_KOVERS.transpose(), all_KOVERS])
    temp = mm.reverse(temp)
    c_all = mm.dot([temp, all_KOVERS.transpose(), delta_KOVERS])
    return c_all


def findC2(light_source, C_First, all_KOVERS, r_std, r_sub, maxRMS):
    all_E = []
    blue_KOVERS = np.array([all_KOVERS.T[0]]).T
    red_KOVERS = np.array([all_KOVERS.T[1]]).T
    yellow_KOVERS = np.array([all_KOVERS.T[2]]).T

    k_std = mm.applyFunction(r_std, find_KOVERS)
    k_sub = mm.applyFunction(r_sub, find_KOVERS)

    KOVERS_First = mm.sum([C_First[0]*blue_KOVERS, C_First[1] *
                           red_KOVERS, C_First[2]*yellow_KOVERS, k_sub])
    R_First = mm.applyFunction(KOVERS_First, find_r)
    k_est = np.vectorize(find_KOVERS)(R_First)

    alpha1 = np.subtract(r_std, R_First)
    alpha2 = np.subtract(k_std, k_est)
    D = np.diagflat(alpha1/alpha2)
    CC = np.array([C_First[0], C_First[1], C_First[2]])
    pi = all_KOVERS
    T = np.hstack((light_source.xbar, light_source.ybar,
                   light_source.zbar)).transpose()
    E_Diag = np.diagflat(light_source.E)
    r_new = []
    ESTN = Spectrum(R_First, light_source)
    STD = Spectrum(r_std, light_source)
    for i in range(1000):
        #print("E", i, delta_E(ESTN, STD))
        D_E = delta_E(ESTN, STD)
        if D_E <= maxRMS:
            break
        all_E.append(D_E)
        delta_t = find_delta_t(STD, ESTN)
        TEDPI = mm.dot([T, E_Diag, D, pi])
        delta_c = mm.dot([mm.reverse(TEDPI), delta_t])
        CC = np.add(CC, delta_c)
        k_new = mm.sum([CC[0]*blue_KOVERS, CC[1] *
                        red_KOVERS, CC[2] * yellow_KOVERS, k_sub])
        r_new = mm.applyFunction(k_new, find_r)
        ESTN = Spectrum(r_new, light_source)
    return [CC, all_E, i]


def autolabel(rects, ax, xpos='center', p=6):
    """
    Attach a text label above each bar in *rects*, displaying its height.
    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """
    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0, 'right': 1, 'left': -1}

    for rect in rects:
        height = rect.get_height()
        height = round(height, p)
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(offset[xpos]*3, 3),  # use 3 points offset
                    textcoords="offset points",  # in both directions
                    ha=ha[xpos], va='bottom')


"""
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
