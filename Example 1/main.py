import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Spectrum import *
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
wave_length = mm.array_distance(400, distance, 700)
example_data = mm.array_repeat(1, data_size)

"""
Start Getting Data From Excel file
"""
data = pd.read_excel('data.xls')

extract_data = pd.DataFrame(data, columns=['c'])
c = extract_data.to_numpy()
c = mm.cleanNaN(c)

R_red = []
for i in range(1, 8):
    extract_data = pd.DataFrame(data, columns=['r'+str(i)])
    R_red.append(extract_data.to_numpy())

R_yellow = []
for i in range(1, 8):
    extract_data = pd.DataFrame(data, columns=['y'+str(i)])
    R_yellow.append(extract_data.to_numpy())

R_blue = []
for i in range(1, 8):
    extract_data = pd.DataFrame(data, columns=['b'+str(i)])
    R_blue.append(extract_data.to_numpy())


extract_data = pd.DataFrame(data, columns=['Rsub'])
r_sub = extract_data.to_numpy()
k_sub = mm.applyFunction(r_sub, find_KOVERS)

extract_data = pd.DataFrame(data, columns=['Rstd'])
r_std = extract_data.to_numpy()
k_std = mm.applyFunction(r_std, find_KOVERS)

extract_data = pd.DataFrame(data, columns=['xbar'])
xbar = extract_data.to_numpy()

extract_data = pd.DataFrame(data, columns=['ybar'])
ybar = extract_data.to_numpy()

extract_data = pd.DataFrame(data, columns=['zbar'])
zbar = extract_data.to_numpy()

viewer = Viewer(xbar, ybar, zbar)

extract_data = pd.DataFrame(data, columns=['D65'])
E_D65 = extract_data.to_numpy()
light_source = LightSource(E_D65, viewer)

"""
Getting Data is Finished
"""

# initial object to find K OVER S for Blue Dye
BBB = Dye(blue_sample_num, data_size)
BBB.setR(r_sub)
BBB.setC(c)
BBB.setSub(r_sub)
BBB.setR(R_blue)


# initial object to find K OVER S for Red Dye
RRR = Dye(red_sample_num, data_size)
RRR.setR(r_sub)
RRR.setC(c)
RRR.setSub(r_sub)
RRR.setR(R_red)


# initial object to find K OVER S for Yellow Dye
YYY = Dye(yellow_sample_num, data_size)
YYY.setR(r_sub)
YYY.setC(c)
YYY.setSub(r_sub)
YYY.setR(R_yellow)

blue_KOVERS = BBB.getKOVERS()
red_KOVERS = RRR.getKOVERS()
yellow_KOVERS = YYY.getKOVERS()


all_KOVERS = np.hstack((blue_KOVERS, red_KOVERS, yellow_KOVERS))
delta_KOVERS = mm.sum([k_std, -1*k_sub])
C_First = findC1(all_KOVERS, delta_KOVERS)


KOVERS_First = mm.sum([C_First[0]*blue_KOVERS, C_First[1] *
                       red_KOVERS, C_First[2]*yellow_KOVERS, k_sub])
R_First = mm.applyFunction(KOVERS_First, find_r)


EST = Spectrum(R_First, light_source)
STD = Spectrum(r_std, light_source)
RMS_First = RMS(R_First, r_std)
DeltaE_First = delta_E(EST, STD)

Data3 = findC2(light_source, C_First, all_KOVERS, r_std, r_sub, maxRMS)
C_Last = Data3[0]
all_E = Data3[1]
num_tried = Data3[2]

KOVERS_Last = mm.sum([C_Last[0]*blue_KOVERS, C_Last[1] *
                      red_KOVERS, C_Last[2]*yellow_KOVERS, k_sub])
R_Last = mm.applyFunction(KOVERS_Last, find_r)
ESTN = Spectrum(R_Last, light_source)

RMS_Last = RMS(R_Last, r_std)
DeltaE_Last = delta_E(ESTN, STD)

"""

Showing Results

"""
# Draw K / S For 3 Dyes
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

# Draw R For Diffrent Methodes
fig, axs = plt.subplots(3)
axs[0].set_title('Method 1')
axs[0].plot(wave_length, r_std, color='green')
axs[0].plot(wave_length, R_First, color='red')

axs[1].set_title('Method 2')
axs[1].plot(wave_length, r_std, color='green')
axs[1].plot(wave_length, R_First, color='blue')

axs[2].set_title('All in One')
axs[2].plot(wave_length, r_std, color='green')
axs[2].plot(wave_length, R_First, color='red')
axs[2].plot(wave_length, R_Last, color='blue')

plt.gcf().canvas.set_window_title('Compare شکل2-1')
for ax in axs.flat:
    ax.set(xlabel='Wave Length', ylabel='R')
plt.tight_layout()
fig.set_size_inches(4, 8)
plt.show()

# Draw R Better Way
plt.plot(wave_length, r_std, color='green')
plt.plot(wave_length, R_First, color='red')
plt.plot(wave_length, R_Last, color='blue')
plt.gcf().canvas.set_window_title('Bigger Comparison شکل2-2')
plt.xlabel('Wave Length')
plt.ylabel('R')
plt.gcf().set_size_inches(8, 8)
plt.show()

# Draw Delta E in Loop
plt.plot(mm.array_distance(1, 1, 81), all_E, color='red')
plt.gcf().canvas.set_window_title('Delta E شکل3')
plt.xlabel('Number of tries')
plt.ylabel('Δ E')
plt.show()

# Draw Delta E & RMS & C For 3 Dyes in two methods
m1_means, m1_std = (RMS_First, DeltaE_First,
                    C_First[0][0], C_First[1][0], C_First[2][0]), (0, 0, 0, 0, 0)
m2_means, m2_std = (RMS_Last, DeltaE_Last,
                    C_Last[0][0], C_Last[1][0], C_Last[2][0]), (0, 0, 0, 0, 0)

ind = np.arange(len(m1_means))  # the x locations for the groups
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, m1_means, width, yerr=m1_std,
                label='Method 1')
rects2 = ax.bar(ind + width/2, m2_means, width, yerr=m2_std,
                label='Method 2')

# Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Scores')
#ax.set_title('Scores by STH')
ax.set_xticks(ind)
ax.set_xticklabels(('RMS', 'Δ E', "C Blue", "C Red", "C Yellow"))
ax.legend()
fig.tight_layout()
autolabel(rects1, ax, "center", precise)
autolabel(rects2, ax, "center", precise)
plt.gcf().set_size_inches(10, 10)
plt.gcf().canvas.set_window_title('جدول 1')
plt.show()

#print("C Blue:", C_First[0], "C Red", C_First[1], "C Yellow", C_First[2])
# print(C_Last)
