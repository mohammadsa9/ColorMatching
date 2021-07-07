#!/usr/bin/env python
# coding: utf-8

# # Libraries

# In[1]:


import numpy as np
import numpy.linalg as linalg
import pandas as pd

import matplotlib.pyplot as plt

from scipy.spatial import Delaunay
from pyhull.delaunay import DelaunayTri

from tabulate import tabulate
from prettytable import PrettyTable

from colour.plotting import *


# # My Libraries
# 

# In[2]:


from libraries.Spectrum import *
from libraries.MyPlot import *
import libraries.MyMath as mm


# # Initial Data

# In[3]:


start_wave = 400
end_wave = 700
data_size = 31
distance = (end_wave - start_wave) / (data_size - 1)

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
data = pd.read_excel("data/data.xls")

extract_data = pd.DataFrame(data, columns=["c"])
c = extract_data.to_numpy()
c = mm.cleanNaN(c)

R_red = []
for i in range(1, red_sample_num + 1):
    extract_data = pd.DataFrame(data, columns=["r" + str(i)])
    R_red.append(extract_data.to_numpy())

R_yellow = []
for i in range(1, yellow_sample_num + 1):
    extract_data = pd.DataFrame(data, columns=["y" + str(i)])
    R_yellow.append(extract_data.to_numpy())

R_blue = []
for i in range(1, blue_sample_num + 1):
    extract_data = pd.DataFrame(data, columns=["b" + str(i)])
    R_blue.append(extract_data.to_numpy())


extract_data = pd.DataFrame(data, columns=["Rsub"])
R_sub = extract_data.to_numpy()
k_sub = mm.applyFunction(R_sub, find_KOVERS)

extract_data = pd.DataFrame(data, columns=["Rstd"])
R_std = extract_data.to_numpy()
k_std = mm.applyFunction(R_std, find_KOVERS)

extract_data = pd.DataFrame(data, columns=["xbar"])
xbar = extract_data.to_numpy()

extract_data = pd.DataFrame(data, columns=["ybar"])
ybar = extract_data.to_numpy()

extract_data = pd.DataFrame(data, columns=["zbar"])
zbar = extract_data.to_numpy()

viewer = Viewer(xbar, ybar, zbar)

extract_data = pd.DataFrame(data, columns=["D65"])
E_D65 = extract_data.to_numpy()
light_source = LightSource(E_D65)

# Principal Component Analysis
munsell_size = 1269
munsell_R = []
data = pd.read_excel("data/Munsell400_10_700.xlsx")
for i in range(munsell_size):
    newdata = data[data.columns[i]]
    newdata = newdata.to_numpy()
    munsell_R.append(newdata)
munsell_R = np.array(munsell_R).T

munsell_A = np.cov(munsell_R)
eigenValues, eigenVectors = linalg.eig(munsell_A)
idx = eigenValues.argsort()[::-1]
eigenValues = eigenValues[idx]
eigenVectors = eigenVectors[:, idx]

R_mean = np.array([[sum(row) for row in munsell_R]])
R_mean = R_mean.T / munsell_size


# # K OVER S for dyes

# In[4]:


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


# # Analyze current data

# In[5]:


print(
    tabulate(
        {"λ (nm)": wave_length, "R mean": R_mean}, headers="keys", tablefmt="fancy_grid"
    )
)

# print(eigenVectors)
print("\n3 EigenVectors with Highest EigenValues")
print(eigenVectors[0:3])
print("\n3 Highest EigenValues")
print(eigenValues[0:3])
print()

(p1,) = plt.plot(wave_length, BBB.getR(), color="blue", label="R Blue Dye")
(p2,) = plt.plot(wave_length, RRR.getR(), color="red", label="R Red Dye")
(p3,) = plt.plot(wave_length, YYY.getR(), color="yellow", label="R Yellow Dye")
lines = [p1, p2, p3]
draw_R_style1(lines)

print()

(p1,) = plt.plot(wave_length, BBB.getKOVERS(), color="blue", label="K/S Blue Dye")
(p2,) = plt.plot(wave_length, RRR.getKOVERS(), color="red", label="K/S Red Dye")
(p3,) = plt.plot(wave_length, YYY.getKOVERS(), color="yellow", label="K/S Yellow Dye")
lines = [p1, p2, p3]
draw_KoverS_style1(lines)

print()

# Plotting the *CIE 1931 Chromaticity Diagram*.
plot_chromaticity_diagram_CIE1931(standalone=False)

# D65
OBS = Observation(light_source, viewer, 1)
xy_D65 = OBS.getxy()
xy = xy_D65
x, y = xy
plt.plot(x, y, "o-", color="black")
plt.annotate(
    "D65",
    xy=xy,
    xytext=(-50, 30),
    textcoords="offset points",
    arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=-0.2"),
)
# Blue
OBS = Observation(light_source, viewer, BBB.getR())
xy_D65 = OBS.getxy()
xy = xy_D65
x, y = xy
plt.plot(x, y, "o-", color="black")
plt.annotate(
    "Blue dye",
    xy=xy,
    xytext=(-50, 30),
    textcoords="offset points",
    arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=-0.2"),
)

# Red
OBS = Observation(light_source, viewer, RRR.getR())
xy_D65 = OBS.getxy()
xy = xy_D65
x, y = xy
plt.plot(x, y, "o-", color="black")
plt.annotate(
    "Red dye",
    xy=xy,
    xytext=(-50, 30),
    textcoords="offset points",
    arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=-0.2"),
)

# Yellow
OBS = Observation(light_source, viewer, YYY.getR())
xy_D65 = OBS.getxy()
xy = xy_D65
x, y = xy
plt.plot(x, y, "o-", color="black")
plt.annotate(
    "Yellow dye",
    xy=xy,
    xytext=(-50, 30),
    textcoords="offset points",
    arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=-0.2"),
)
# Displaying the plot.
render(standalone=True, limits=(-0.1, 0.9, -0.1, 0.9), x_tighten=True, y_tighten=True)


# # SET R Substrate for New Surface = 1

# In[6]:


R_sub = np.array([mm.array_repeat(1, data_size)]).T


# # Creating Look up Table

# In[7]:


pr = 11
Dis1 = np.linspace(0, 1, pr)
Dis2 = np.linspace(0, 1, pr)
Dis3 = np.linspace(0, 1, pr)


R_Lookup = []
XYZ_Lookup = []
C_Lookup = []

Mix = Mixture(R_sub)
for x in range(pr):
    for y in range(pr):
        for z in range(pr):
            Mix.clear()
            Mix.add(Dis1[x], blue_KOVERS)
            Mix.add(Dis2[y], red_KOVERS)
            Mix.add(Dis3[z], yellow_KOVERS)
            Temp = Observation(light_source, viewer, Mix.getR())
            XYZ_Lookup.append([Temp.getX(), Temp.getY(), Temp.getZ()])
            C_Lookup.append([Dis1[x], Dis2[y], Dis3[z]])
            R_Lookup.append(Mix.getR().T[0])

R_Lookup = np.array(R_Lookup)
XYZ_Lookup = np.array(XYZ_Lookup)
C_Lookup = np.array(C_Lookup)

print(Dis1)


# # Creating Samples to Evaluate Methods

# In[8]:


pr = 4
Dis1 = np.array([0.25, 0.45, 0.65, 0.85])
Dis2 = np.array([0.25, 0.45, 0.65, 0.85])
Dis3 = np.array([0.25, 0.45, 0.65, 0.85])

R_Samples = []
XYZ_Samples = []
C_Samples = []

Mix = Mixture(R_sub)
for x in range(pr):
    for y in range(pr):
        for z in range(pr):
            Mix.clear()
            Mix.add(Dis1[x], blue_KOVERS)
            Mix.add(Dis2[y], red_KOVERS)
            Mix.add(Dis3[z], yellow_KOVERS)
            Temp = Observation(light_source, viewer, Mix.getR())
            XYZ_Samples.append([Temp.getX(), Temp.getY(), Temp.getZ()])
            C_Samples.append([Dis1[x], Dis2[y], Dis3[z]])
            R_Samples.append(Mix.getR().T[0])

count = 1
for R_Sample in R_Samples:
    (p1,) = plt.plot(
        wave_length, R_Sample, color="green", label="R Sample " + str(count)
    )
    lines = [p1]
    count = count + 1
    # print("GFC: ",GFC(munsell_R,R_Sample))

draw_R_style1(lines)
print(Dis1)


# # Method 1 Interpolation using R with reduced dimension Based on PCA (3D) => Project purpose

# In[9]:


dim = 3
R_Lookup = mm.array_PC(R_Lookup, dim, eigenVectors, R_mean)
calc = MyDelaunay(R_Lookup)

""" Another method to calculate with delaunay
    from scipy.interpolate import NearestNDInterpolator
    why = LinearNDInterpolator(R_Lookup,C_Lookup)
    why = NearestNDInterpolator(R_Lookup,C_Lookup)
    print("lib",why(R_Find))
"""

# Result
M_R_RMS = 0
M_R_DeltaE = 0
M_R_DeltaC = 0

M_R_minRMS = 0
M_R_maxRMS = 0

M_R_minE = 0
M_R_maxE = 0

for i in range(len(R_Samples)):
    # print(R_Samples[i])
    R_std = R_Samples[i]
    R_Find = R_Samples[i]
    R_Find = mm.PC(R_Find, dim, eigenVectors, R_mean)

    res = calc.getResult(R_Find, C_Lookup)
    C_Inter2 = res[0]

    Mix.clear()
    Mix.add(C_Inter2[0], blue_KOVERS)
    Mix.add(C_Inter2[1], red_KOVERS)
    Mix.add(C_Inter2[2], yellow_KOVERS)
    R_Inter2 = Mix.getR()

    Inter2 = Observation(light_source, viewer, R_Inter2)
    STD = Observation(light_source, viewer, R_std)

    compare_4 = Compare(STD, Inter2)
    RMS_Inter2 = compare_4.RMS()
    DeltaE_Inter2 = compare_4.delta_E()
    GFC = 0

    text_R = "R: " + str(R_std)
    text_RMS = "RMS: " + str(RMS_Inter2)
    text_DeltaE = "ΔE: " + str(DeltaE_Inter2)
    text_GFC = "GFC: " + str(GFC)
    text_C_Real = "Real C: " + str(C_Samples[i])
    text_C_Cal = "Interpolated C: " + str(C_Inter2)

    text_all = (
        text_RMS
        + "\n"
        + text_DeltaE
        + "\n"
        + text_GFC
        + "\n"
        + text_C_Real
        + "\n"
        + text_C_Cal
        + "\n\n"
        + text_R
    )

    # Result
    M_R_RMS = M_R_RMS + RMS_Inter2
    M_R_DeltaE = M_R_DeltaE + DeltaE_Inter2
    M_R_DeltaC = M_R_DeltaC + mm.RMS(C_Samples[i], C_Inter2)

    if i == 0:
        M_R_minRMS = RMS_Inter2
        M_R_maxRMS = RMS_Inter2

        M_R_minE = DeltaE_Inter2
        M_R_maxE = DeltaE_Inter2

    if RMS_Inter2 < M_R_minRMS:
        M_R_minRMS = RMS_Inter2

    if RMS_Inter2 > M_R_maxRMS:
        M_R_maxRMS = RMS_Inter2

    if DeltaE_Inter2 < M_R_minE:
        M_R_minE = DeltaE_Inter2

    if DeltaE_Inter2 > M_R_maxE:
        M_R_maxE = DeltaE_Inter2

    (p1,) = plt.plot(wave_length, R_std, color="green", label="R Sample " + str(i + 1))
    (p2,) = plt.plot(wave_length, R_Inter2, color="black", label="R Interpolated (PCA)")
    lines = [p1, p2]
    draw_R_style1(lines, comment=text_all)

# Result
M_R_RMS = M_R_RMS / len(R_Samples)
M_R_DeltaE = M_R_DeltaE / len(R_Samples)
M_R_DeltaC = M_R_DeltaC / len(R_Samples)

print("mean RMS: ", M_R_RMS)
print("mean ΔE: ", M_R_DeltaE)
print("mean ΔC: ", M_R_DeltaC)
print("Min RMS: ", M_R_minRMS)
print("Max RMS: ", M_R_maxRMS)
print("Min ΔE: ", M_R_minE)
print("Max ΔE: ", M_R_maxE)


# # Method 2 Interpolation using XYZ => For Comparison

# In[10]:


calc = MyDelaunay(XYZ_Lookup)

# Result
M_XYZ_RMS = 0
M_XYZ_DeltaE = 0
M_XYZ_DeltaC = 0

M_XYZ_minRMS = 0
M_XYZ_maxRMS = 0

M_XYZ_minE = 0
M_XYZ_maxE = 0

for i in range(len(R_Samples)):
    # print(R_Samples[i])
    R_std = R_Samples[i]
    Find = XYZ_Samples[i]
    Temp = Observation(light_source, viewer, R_std)

    res = calc.getResult(Find, C_Lookup)
    C_Inter1 = res[0]

    Mix.clear()
    Mix.add(C_Inter1[0], blue_KOVERS)
    Mix.add(C_Inter1[1], red_KOVERS)
    Mix.add(C_Inter1[2], yellow_KOVERS)
    R_Inter1 = Mix.getR()
    Inter1 = Observation(light_source, viewer, R_Inter1)
    STD = Observation(light_source, viewer, R_std)
    compare_3 = Compare(Inter1, STD)
    RMS_Inter1 = compare_3.RMS()
    DeltaE_Inter1 = compare_3.delta_E()

    text_R = "R: " + str(R_std)
    text_RMS = "RMS: " + str(RMS_Inter1)
    text_DeltaE = "ΔE: " + str(DeltaE_Inter1)
    text_GFC = "GFC: " + str("?")
    text_C_Real = "Real C: " + str(C_Samples[i])
    text_C_Cal = "Interpolated C: " + str(C_Inter1)

    text_all = (
        text_RMS
        + "\n"
        + text_DeltaE
        + "\n"
        + text_GFC
        + "\n"
        + text_C_Real
        + "\n"
        + text_C_Cal
        + "\n\n"
        + text_R
    )

    # Result
    M_XYZ_RMS = M_XYZ_RMS + RMS_Inter1
    M_XYZ_DeltaE = M_XYZ_DeltaE + DeltaE_Inter1
    M_XYZ_DeltaC = M_XYZ_DeltaC + mm.RMS(C_Samples[i], C_Inter1)

    if i == 0:
        M_XYZ_minRMS = RMS_Inter1
        M_XYZ_maxRMS = RMS_Inter1

        M_XYZ_minE = DeltaE_Inter1
        M_XYZ_maxE = DeltaE_Inter1

    if RMS_Inter1 < M_XYZ_minRMS:
        M_XYZ_minRMS = RMS_Inter1

    if RMS_Inter1 > M_XYZ_maxRMS:
        M_XYZ_maxRMS = RMS_Inter1

    if DeltaE_Inter1 < M_XYZ_minE:
        M_XYZ_minE = DeltaE_Inter1

    if DeltaE_Inter1 > M_XYZ_maxE:
        M_XYZ_maxE = DeltaE_Inter1

    (p1,) = plt.plot(wave_length, R_std, color="green", label="R Sample " + str(i + 1))
    (p2,) = plt.plot(wave_length, R_Inter1, color="black", label="R Interpolated (XYZ)")
    lines = [p1, p2]
    draw_R_style1(lines, comment=text_all)

# Result
M_XYZ_RMS = M_XYZ_RMS / len(R_Samples)
M_XYZ_DeltaE = M_XYZ_DeltaE / len(R_Samples)
M_XYZ_DeltaC = M_XYZ_DeltaC / len(R_Samples)

print("mean RMS: ", M_XYZ_RMS)
print("mean ΔE: ", M_XYZ_DeltaE)
print("mean ΔC: ", M_XYZ_DeltaC)
print("Min RMS: ", M_XYZ_minRMS)
print("Max RMS: ", M_XYZ_maxRMS)
print("Min ΔE: ", M_XYZ_minE)
print("Max ΔE: ", M_XYZ_maxE)


# # Showing Results

# In[19]:


x = PrettyTable()

x.field_names = [
    "Method name",
    "Mean RMS",
    "Mean ΔE",
    "Mean ΔC",
    "Min RMS",
    "Max RMS",
    "Min ΔE",
    "Max ΔE",
]

M_R_RMS = round(M_R_RMS, 5)
M_R_DeltaE = round(M_R_DeltaE, 5)
M_R_DeltaC = round(M_R_DeltaC, 5)
M_R_minRMS = round(M_R_minRMS, 5)
M_R_maxRMS = round(M_R_maxRMS, 5)
M_R_minE = round(M_R_minE, 5)
M_R_maxE = round(M_R_maxE, 5)

M_XYZ_RMS = round(M_XYZ_RMS, 5)
M_XYZ_DeltaE = round(M_XYZ_DeltaE, 5)
M_XYZ_DeltaC = round(M_XYZ_DeltaC, 5)
M_XYZ_minRMS = round(M_XYZ_minRMS, 5)
M_XYZ_maxRMS = round(M_XYZ_maxRMS, 5)
M_XYZ_minE = round(M_XYZ_minE, 5)
M_XYZ_maxE = round(M_XYZ_maxE, 5)

x.add_row(
    [
        "Principal Component Coordinates",
        M_R_RMS,
        M_R_DeltaE,
        M_R_DeltaC,
        M_R_minRMS,
        M_R_maxRMS,
        M_R_minE,
        M_R_maxE,
    ]
)

x.add_row(
    [
        "XYZ",
        M_XYZ_RMS,
        M_XYZ_DeltaE,
        M_XYZ_DeltaC,
        M_XYZ_minRMS,
        M_XYZ_maxRMS,
        M_XYZ_minE,
        M_XYZ_maxE,
    ]
)

print(x)


# In[ ]:




