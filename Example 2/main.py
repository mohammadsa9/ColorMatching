import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

# My Libraries
from Functions import *
from MyPlot import *
import MyMath as mm

"""
Initial Data
"""
start_wave = 400
end_wave = 700
distance = 10
data_size = int((end_wave - start_wave)/distance)+1
munsell_size = 1269
munsell_size2 = 32
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
munsell_R = []
data = pd.read_excel('Munsell400_10_700.xlsx')
for i in range(munsell_size):
    newdata = data[data.columns[i]]
    newdata = newdata.to_numpy()
    #newdata = np.array([newdata])
    munsell_R.append(newdata)


data = pd.read_excel('sample-data.xlsx')
extract_data = pd.DataFrame(data, columns=['xbar'])
xbar = extract_data.to_numpy()
xbar = mm.cleanNaN(xbar)

extract_data = pd.DataFrame(data, columns=['ybar'])
ybar = extract_data.to_numpy()
ybar = mm.cleanNaN(ybar)

extract_data = pd.DataFrame(data, columns=['zbar'])
zbar = extract_data.to_numpy()
zbar = mm.cleanNaN(zbar)

viewer = Viewer(xbar, ybar, zbar)

extract_data = pd.DataFrame(data, columns=['D65'])
E_D65 = extract_data.to_numpy()
E_D65 = mm.cleanNaN(E_D65)
light_source = LightSource(E_D65)

samples = []
for i in range(1, 15):
    extract_data = pd.DataFrame(data, columns=['s'+str(i)])
    samples.append(mm.cleanNaN(extract_data.to_numpy()))
"""
Getting Data is Finished
"""

# Calculate XYZ for Munsell R values
munsell_XYZ = []
temp = Observation(light_source, viewer, 1)

for i in range(munsell_size):
    temp = Observation(light_source, viewer, munsell_R[i])
    temp_arr = [temp.getX(), temp.getY(), temp.getZ()]
    munsell_XYZ.append(temp_arr)


# Draw Munsell Data for a better sense
for i in range(munsell_size):
    plt.plot(wave_length, munsell_R[i])
plt.gcf().canvas.set_window_title('All Munsell Data')
plt.xlabel('Wave Length')
plt.ylabel('R')
plt.gcf().set_size_inches(8, 8)
plt.show()


# Part 1 - Interpolating R from XYZ
# Delaunay Calculation
points = np.array(munsell_XYZ)
calc = MyDelaunay(points)
for i in range(14):
    OBS = Observation(light_source, viewer, samples[i])
    R_calc, check = Interploration(calc, OBS, munsell_R)
    # R_calc2 = myInterploration(OBS, munsell_XYZ, munsell_R)
    CAL = Observation(light_source, viewer, R_calc)
    compare = Compare(OBS, CAL)
    RMS = compare.RMS()
    DeltaE = compare.delta_E()
    gamut = 'in-gamut data'
    if check == False:
        gamut = 'out-of-gamut data'
    text = 'Sample'+str(i+1)+': '+gamut+' , ΔE = ' + \
        str(DeltaE)+', RMS = '+str(round(RMS, 6))
    R_sample = OBS.R[0]
    # Showing Data
    p1, = plt.plot(wave_length, R_sample, color='green', label="R STD")
    p2, = plt.plot(wave_length, R_calc, color='red', label="R Interpolation")
    # p3, = plt.plot(wave_length, R_calc2[0], color='black', label="R 3")
    lines = [p1, p2]
    plt.legend(lines, [l.get_label() for l in lines])
    plt.gcf().canvas.set_window_title('Delaunay Calculation - Sample: '+str(i+1))
    plt.xlabel('Wave Length')
    plt.ylabel('R')
    plt.gcf().set_size_inches(8, 8)
    plt.ylim([0, 2])
    plt.text(400, 1.5, r''+text)
    plt.show()


# Part 2 - Interpolating XYZ from R
test = np.array([munsell_R[1]]).T
munsell_R = munsell_R[0:munsell_size2]
points = np.array(munsell_R)
calc = MyDelaunay(points)
for i in range(14):
    OBS = Observation(light_source, viewer, samples[i])
    # OBS = Observation(light_source, viewer, test)
    XYZ_calc, check = revInterploration(calc, OBS, munsell_XYZ)
    gamut = 'in-gamut data'
    if check is False:
        gamut = 'out-of-gamut data'

    # Show Result
    m1_means, m1_std = (OBS.getX(), OBS.getY(), OBS.getZ()), (0, 0, 0)
    m2_means, m2_std = (XYZ_calc[0], XYZ_calc[1],
                        XYZ_calc[2]), (0, 0, 0)

    ind = np.arange(len(m1_means))  # the x locations for the groups
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(ind - width/2, m1_means, width, yerr=m1_std,
                    label='STD')
    rects2 = ax.bar(ind + width/2, m2_means, width, yerr=m2_std,
                    label='Interpolation')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel(gamut)
    # ax.set_title('Scores by STH')
    ax.set_xticks(ind)
    ax.set_xticklabels(('X', 'Y', "Z"))
    ax.legend()
    fig.tight_layout()
    autolabel(rects1, ax, "center", precise)
    autolabel(rects2, ax, "center", precise)
    plt.gcf().set_size_inches(10, 10)
    plt.gcf().canvas.set_window_title('Sample '+str(i+1))
    plt.show()
