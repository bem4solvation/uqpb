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
n_test = 5000
n_workers = 10
n_samples = n_test//n_workers
exec_sampler = os.path.join("..","..","sampler.py")

# 1. Crear carpetas para guardar los resultados

sampler_list = ["pseudo", "LHS", "Halton", "Hammersley", "Sobol"]
folder = "tests"
for sampler in sampler_list:
    folder_name = os.path.join(folder, create_folder_name(pqr_file,sampler))
    print("Creating folder: ", folder_name)
    os.system("python {} -nt {} -nw {} -pqr {} -f {} -sp {}".format(
                exec_sampler,n_test, n_workers, pqr_file, folder_name, sampler)
                )
    
# 2. Leer resultados de coeficiente de agitación de radio, primera columna de los archivos _coeff.txt

i_molecule = 0
# data es una lista con una matriz por sampler
# cada matriz tiene n_workers filas y n_samples columnas
data = [numpy.zeros((n_samples,n_workers)) for _ in range(len(sampler_list))]
for sampler in sampler_list:
    folder_name = os.path.join(folder, create_folder_name(pqr_file,sampler))
    i = sampler_list.index(sampler)

    for job in range(n_workers):
        job_folder = os.path.join(folder_name, f"job_{job:02}")
        list_of_coeff = os.listdir(job_folder)
        list_of_coeff = [x for x in list_of_coeff if x.endswith("_coeff.txt")]
        m = 0
        for coeff_file in list_of_coeff:
            with open(os.path.join(job_folder, coeff_file), "r") as f:
                k = 0
                while k<=i_molecule:
                    line = f.readline()
                    k += 1
            coeff = float(line.strip().split()[0]) # 0 por valor del radio
            data[i][m,job] = coeff
            m+=1

# 3. Se guardan los datos de data en en archivos de texto

for i,sampler in enumerate(sampler_list):
    output_folder = "random_coefficients"
    if output_folder not in os.listdir():
        os.mkdir(output_folder)
    file_name = sampler + "_coeff.txt"
    # Save file
    numpy.savetxt(os.path.join(output_folder,file_name),data[i],delimiter=',',newline='\n')








