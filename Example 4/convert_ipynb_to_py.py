import os


desired_path = os.path.dirname(os.path.realpath(__file__))

os.chdir(desired_path)

list_of_directory = os.listdir(desired_path)

for file in list_of_directory:
    if ".ipynb" in str(file):
        os.system("ipython nbconvert --to script " + str(file))
