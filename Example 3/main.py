import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from pyhull.delaunay import DelaunayTri
# My Libraries
from Spectrum import *
from MyPlot import *
import MyMath as mm

"""
Initial Data
"""
start_wave = 400
end_wave = 700
data_size = 31
distance = (end_wave - start_wave)/(data_size-1)

blue_sample_num = 7
yellow_sample_num = 7
red_sample_num = 7

maxRMS = 0.001
precise = 6
"""
Creating Wave lengths array for plots
Creating Example Data for checking plots
"""
# [400, 410, 420, ..., 700]
wave_length = mm.array_distance(400, distance, 700)

# [1, 1, 1, ..., 1]
example_data = mm.array_repeat(1, data_size)

"""
Start Getting Data From Excel file
"""
data = pd.read_excel('data.xls')

extract_data = pd.DataFrame(data, columns=['c'])
c = extract_data.to_numpy()
c = mm.cleanNaN(c)

R_red = []
for i in range(1, red_sample_num+1):
    extract_data = pd.DataFrame(data, columns=['r'+str(i)])
    R_red.append(extract_data.to_numpy())

R_yellow = []
for i in range(1, yellow_sample_num+1):
    extract_data = pd.DataFrame(data, columns=['y'+str(i)])
    R_yellow.append(extract_data.to_numpy())

R_blue = []
for i in range(1, blue_sample_num+1):
    extract_data = pd.DataFrame(data, columns=['b'+str(i)])
    R_blue.append(extract_data.to_numpy())


extract_data = pd.DataFrame(data, columns=['Rsub'])
R_sub = extract_data.to_numpy()
k_sub = mm.applyFunction(R_sub, find_KOVERS)

extract_data = pd.DataFrame(data, columns=['Rstd'])
R_std = extract_data.to_numpy()
k_std = mm.applyFunction(R_std, find_KOVERS)

extract_data = pd.DataFrame(data, columns=['xbar'])
xbar = extract_data.to_numpy()

extract_data = pd.DataFrame(data, columns=['ybar'])
ybar = extract_data.to_numpy()

extract_data = pd.DataFrame(data, columns=['zbar'])
zbar = extract_data.to_numpy()

viewer = Viewer(xbar, ybar, zbar)

extract_data = pd.DataFrame(data, columns=['D65'])
E_D65 = extract_data.to_numpy()
light_source = LightSource(E_D65)

"""
Getting Data is Finished
"""

# initial object to find K OVER S for Blue Dye
BBB = Dye(blue_sample_num, data_size)
BBB.setR(R_blue)
BBB.setC(c)
BBB.setSub(R_sub)


# initial object to find K OVER S for Red Dye
RRR = Dye(red_sample_num, data_size)
RRR.setR(R_red)
RRR.setC(c)
RRR.setSub(R_sub)


# initial object to find K OVER S for Yellow Dye
YYY = Dye(yellow_sample_num, data_size)
YYY.setR(R_yellow)
YYY.setC(c)
YYY.setSub(R_sub)

blue_KOVERS = BBB.getKOVERS()
red_KOVERS = RRR.getKOVERS()
yellow_KOVERS = YYY.getKOVERS()

# Method 1
all_KOVERS = np.hstack((blue_KOVERS, red_KOVERS, yellow_KOVERS))
delta_KOVERS = mm.sum([k_std, -1*k_sub])
C_First = findC1(all_KOVERS, delta_KOVERS)


First = Mixture(R_sub)
First.add(C_First[0], blue_KOVERS)
First.add(C_First[1], red_KOVERS)
First.add(C_First[2], yellow_KOVERS)
KOVERS_First = First.getKOVERS()
R_First = First.getR()

EST = Observation(light_source, viewer, R_First)
STD = Observation(light_source, viewer, R_std)
compare_1 = Compare(EST, STD)
RMS_First = compare_1.RMS()
DeltaE_First = compare_1.delta_E()

# Method 2
Data3 = findC2(STD, R_sub, C_First, all_KOVERS, maxRMS)
C_Last = Data3[0]
all_E = Data3[1]
num_tried = Data3[2]


Last = Mixture(R_sub)
Last.add(C_Last[0], blue_KOVERS)
Last.add(C_Last[1], red_KOVERS)
Last.add(C_Last[2], yellow_KOVERS)
KOVERS_Last = Last.getKOVERS()
R_Last = Last.getR()

ESTN = Observation(light_source, viewer, R_Last)
compare_2 = Compare(STD, ESTN)
RMS_Last = compare_2.RMS()
DeltaE_Last = compare_2.delta_E()

# Method 3 Interpolation using XYZ
pr = 10
Dis1 = np.linspace(0*C_First[0], 2*C_First[0], pr)
Dis2 = np.linspace(0*C_First[1], 2*C_First[1], pr)
Dis3 = np.linspace(0*C_First[2], 2*C_First[2], pr)

XYZ_Lookup = []
C_Lookup = []

Mix = Mixture(R_sub)
for x in range(pr):
    for y in range(pr):
        for z in range(pr):
            Mix.clear()
            Mix.add(Dis1[x][0], blue_KOVERS)
            Mix.add(Dis2[y][0], red_KOVERS)
            Mix.add(Dis3[z][0], yellow_KOVERS)
            Temp = Observation(light_source, viewer, Mix.getR())
            XYZ_Lookup.append([Temp.getX(), Temp.getY(), Temp.getZ()])
            C_Lookup.append([Dis1[x][0], Dis2[y][0], Dis3[z][0]])

XYZ_Lookup = np.array(XYZ_Lookup)
C_Lookup = np.array(C_Lookup)

Temp = Observation(light_source, viewer, R_std)
Find = [Temp.getX(), Temp.getY(), Temp.getZ()]
calc = MyDelaunay(XYZ_Lookup)
res = calc.getResult(Find, C_Lookup)
C_Inter1 = res[0]

Mix.clear()
Mix.add(C_Inter1[0], blue_KOVERS)
Mix.add(C_Inter1[1], red_KOVERS)
Mix.add(C_Inter1[2], yellow_KOVERS)
R_Inter1 = Mix.getR()
Inter1 = Observation(light_source, viewer, Mix.getR())
compare_3 = Compare(Inter1, STD)
RMS_Inter1 = compare_3.RMS()
DeltaE_Inter1 = compare_3.delta_E()


# Method 4 Interpolation using R
pr = 10
Dis1 = np.linspace(0*C_First[0], 2*C_First[0], pr)
Dis2 = np.linspace(0*C_First[1], 2*C_First[1], pr)
Dis3 = np.linspace(0*C_First[2], 2*C_First[2], pr)

R_Lookup = []
C_Lookup = []

Mix = Mixture(R_sub)
for x in range(pr):
    for y in range(pr):
        for z in range(pr):
            Mix.clear()
            Mix.add(Dis1[x][0], blue_KOVERS)
            Mix.add(Dis2[y][0], red_KOVERS)
            Mix.add(Dis3[z][0], yellow_KOVERS)
            R_Lookup.append(Mix.getR().T[0])
            C_Lookup.append([Dis1[x][0], Dis2[y][0], Dis3[z][0]])

R_Lookup = R_Lookup[20:53]
R_Lookup = np.array(R_Lookup)
C_Lookup = np.array(C_Lookup)
calc = MyDelaunay(R_Lookup, 'QJ')
res = calc.getResult(R_std.T[0], C_Lookup)
C_Inter2 = res[0]
print(C_Inter2)
Mix.clear()
Mix.add(C_Inter2[0], blue_KOVERS)
Mix.add(C_Inter2[1], red_KOVERS)
Mix.add(C_Inter2[2], yellow_KOVERS)
R_Inter2 = Mix.getR()
Inter2 = Observation(light_source, viewer, Mix.getR())
compare_4 = Compare(Inter2, STD)
RMS_Inter2 = compare_4.RMS()
DeltaE_Inter2 = compare_4.delta_E()
"""

Showing Results

"""
# Draw K / S For 3 Dyes
"""
fig, axs = plt.subplots(3)
axs[0].plot(wave_length, blue_KOVERS)
axs[0].set_title('Blue')
axs[1].plot(wave_length, red_KOVERS)
axs[1].set_title('Red')
axs[2].plot(wave_length, yellow_KOVERS)
axs[2].set_title('Yellow')
plt.gcf().canvas.set_window_title('K / S For Dyes شکل1')
for ax in axs.flat:
    ax.set(xlabel='Wave Length', ylabel=' K / S ')
plt.tight_layout()
fig.set_size_inches(5, 8)
plt.show()
"""

# Draw all R
p1, = plt.plot(wave_length, R_std, color='green', label="R STD")
p2, = plt.plot(wave_length, R_First, color='red', label="R First Method")
p3, = plt.plot(wave_length, R_Last, color='blue', label="R Second Method")
p4, = plt.plot(wave_length, R_Inter1, color='purple',
               label="R Interpolation using XYZ")
# p5, = plt.plot(wave_length, R_Inter2, color='purple', label = "R Interpolation using R")
lines = [p1, p2, p3, p4]
plt.legend(lines, [l.get_label() for l in lines])
plt.gcf().canvas.set_window_title('Bigger Comparison شکل2-2')
plt.xlabel('Wave Length')
plt.ylabel('R')
plt.gcf().set_size_inches(8, 8)
plt.show()
