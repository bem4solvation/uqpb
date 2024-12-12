import os
import numpy
import pandas
import matplotlib.pyplot as plt
import seaborn
import sys

# Agregra el directorio raiz (experiments) al path
current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, '../../'))
sys.path.append(root_dir)

from experiments.energies_and_samplers.create_sample_folders import run_sampler, create_folder_name

# Setear molecula de prueba y parametros

pqr_file = "../mobley_test_pqr/1112_tetrachloroethane.pqr"
n_test = 1000
n_workers = 1
# Crear carpetas para guardar los resultados

sampler_list = ["pseudo", "LHS", "Halton", "Hammersley", "Sobol"]
folder = "tests"
for sampler in sampler_list:
    folder_name = os.path.join(folder, create_folder_name(pqr_file,sampler))
    print("Creating folder: ", folder_name)
    os.system("python ..\\..\\sampler.py -nt {} -nw {} -pqr {} -f {} -sp {}".format(
                n_test, n_workers, pqr_file, folder_name, sampler)
                )
    


