import numpy
import os
import pbj
import time

"""Este script contiene funciones para calcular la energía electrostática de una molécula de 36 atomos.
En este caso, argnina.pqr
Corre 1000 pruebas para cada sampler y guarda los resultados en archivos csv.
"""

import sys

# Agregra el directorio raiz (experiments) al path
current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, '../../'))
sys.path.append(root_dir)

from experiments.energies_and_samplers.create_sample_folders import  create_folder_name

# Setear molecula de prueba y parametros

pqr_file = "arginina.pqr"

n_test = 1000
n_workers = 1
digits_workers = len(str(n_workers))
mesh_generator = "nanoshaper"
#n_samples = n_test//n_workers
exec_sampler = os.path.join("..","..","sampler.py")
exec_solver = os.path.join("..","..","MC_Solver.py")
exec_get_moments = os.path.join("..","..","get_moments.py")
# 1. Crear carpetas para guardar los resultados

# Boolean value to create the folder and its subfolders with shaked values
bool_create_folders = True

if bool_create_folders:

    sampler_list = ["pseudo", "LHS", "Halton", "Hammersley", "Sobol"]
    #sampler_list = ["pseudo", "Halton", "Hammersley", "Sobol"]
    folder = "tests"
    for sampler in sampler_list:
        folder_name = os.path.join(folder, create_folder_name(pqr_file,sampler))
        print("Creating folder: ", folder_name)
        os.system("python {} -nt {} -nw {} -pqr {} -f {} -sp {}".format(
                    exec_sampler,n_test, n_workers, pqr_file, folder_name, sampler)
                    )
        
# 2. Correr pruebas
i = 0
ss = ""
for sampler in sampler_list:
    folder_name = os.path.join(folder, create_folder_name(pqr_file,sampler))
    #os.system("cls")
    i+=1
    ss+=f"{i} - Reading and executing folder: "+folder_name+'\n'
    print(ss)
    os.system("python {} -f {} -mg {} ".format(
                exec_solver,folder_name,mesh_generator)
                )
    ss+=f"{i} - Getting moments for folder: "+ folder_name,'\n'
    #os.system('cls')
    print(ss)
    os.system("python {} -f {}".format(
                folder_name)
                )
    ss+="Finished folder: "+ folder_name,'\n'

    time.sleep(5)


