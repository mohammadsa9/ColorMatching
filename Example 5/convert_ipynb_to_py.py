import os


desired_path = os.path.dirname(os.path.realpath(__file__))

os.chdir(desired_path)

list_of_directory = os.listdir(desired_path)

for file in list_of_directory:
    if ".ipynb" in str(file):
        # reformat .ipynb file with black
        os.system("black-nb " + str(file))
        # convert .ipynb file to .py file
        os.system("jupyter nbconvert --to script " + str(file))