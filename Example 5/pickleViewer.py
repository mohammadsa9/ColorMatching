import matplotlib.pyplot as plt
import pickle

# Load figure from disk and display
figx = pickle.load(open("FigureObject.fig.pickle", "rb"))

plt.show()
