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
    munsell_R.append(newdata.to_numpy())

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


# Delaunay Calculation
points = np.array(munsell_XYZ)
calc = MyDelaunay(points)
for i in range(14):
    OBS = Observation(light_source, viewer, samples[i])
    R_calc, check = Interploration(calc, OBS, munsell_R)
    CAL = Observation(light_source, viewer, R_calc)
    compare = Compare(OBS, CAL)
    RMS = compare.RMS()
    DeltaE = compare.delta_E()
    gamut = 'in-gamut data'
    if check == False:
        gamut = 'out-of-gamut data'
    text = 'Sample'+str(i+1)+': '+gamut+' , ΔE = ' + \
        str(DeltaE)+', RMS = '+str(round(RMS, 6))
    R_sample = OBS.R
    # Showing Data
    p1, = plt.plot(wave_length, R_sample, color='green', label="R STD")
    p2, = plt.plot(wave_length, R_calc, color='red', label="R Interpolation")
    lines = [p1, p2]
    plt.legend(lines, [l.get_label() for l in lines])
    plt.gcf().canvas.set_window_title('Delaunay Calculation')
    plt.xlabel('Wave Length')
    plt.ylabel('R')
    plt.gcf().set_size_inches(8, 8)
    plt.ylim([0, 2])
    plt.text(400, 1.5, r''+text)
    plt.show()
