


# Comienza programa para calcular energias de solvatacion de molecula informada en input_config.prm del directorio

#Librerias importantes
import os 
import numpy as np
import pbj
from time import time
from utilities import masas
from thermal_functions import nombre_atomo, randomX, new_name, leer_archivo_entrada, shake_file
import gc
from scipy.sparse import csr_matrix, SparseEfficiencyWarning
import warnings
warnings.simplefilter("ignore",SparseEfficiencyWarning)

#Se carga diccionario con datos de input_config.prm
parametros=leer_archivo_entrada('input_config.prm')

# Se comienza escritura de archivo de salida
output_file = parametros['OUTPUT_FILE']
csv_data=open(output_file,'w')
csv_data.write('ITERATION,SOLV. ENERGY,GMRES Iterations,Elapsed time,Number of elements\n')
csv_data.close()

# Average Thermal Length: definicion de funciones y creacion de diccionario
ATL_dic = {}
kB = 1.3806e-23 # Constante de Stefan-Boltzmann
t  = parametros['TIEMPO'] # Tiempo característico
T  = 298 # Temperatura en Kelvin
vt = lambda m,T: np.sqrt(kB*t/m) # thermal velocity en funcion de masa y temperatura (298K por defecto)
L  = lambda m,T: vt(m,T)*t*1e10  #1e10 convierte la unidad a Angstroms
for atomo,masa in masas.items():
    ATL_dic[atomo] = L(masa,T)

# Se agitan las moleculas utilizando los parametros iniciales
path = parametros['TESTPATH']

n_tests         = parametros['N_TESTS']
mainfile        = parametros['PQRFILE']
nombre_molecula = parametros['NAME']

# Se cargan datos de archivo pqr
remarks=""""""
atoms=""""""
end=""""""
arch=open(mainfile)
post_atoms=False
for line in arch:
    data=line.split()
    if data[0]!='ATOM' and not post_atoms:
        remarks+=line
    elif data[0]=='ATOM':
        post_atoms=True
        atoms+=line
    else:
        end+=line
arch.close()


for i in range(n_tests):
    # Se crea nuevo archivo pqr
    test_file=os.path.join(path,new_name(nombre_molecula+'.pqr',i,n_tests))
    shake_file(test_file,i,path,remarks,end,nombre_molecula,atoms,n_tests,ATL_dic)
    
    try:
        start_time=time()
        
        # Carga de archivo pqr 
        molecule=pbj.Solute(test_file,mesh_generator='msms',mesh_density=parametros['DENSIDAD'])
        
        # Parametros del objeto molecule
        molecule.gmres_tolerance        = 1e-5
        molecule.gmres_max_iterations   = 400
        molecule.ep_in                  = parametros['EP_IN']
        molecule.kappa                  = parametros['KAPPA']

        #Calculo de energia de solvatacion electrostatica
        molecule.calculate_solvation_energy()
        
        #Impresion por pantalla
        ET=time()-start_time
        print('INFO: {0}\ti:{1},\t {2}\t(kcal/mol), {3} it.,\t {4} [s],{5} elem.'.format(output_file.split('.')[0],i, 
        molecule.solvation_energy,
        molecule.solver_iteration_count, 
        round(ET,3),
        molecule.mesh.number_of_elements))

        # Escritura en archivo de salida
        csv_data=open(output_file,'a')
        csv_data.write('{0},{1},{2},{3},{4}\n'.format(
            i,
            molecule.solvation_energy,
            molecule.solver_iteration_count,
            round(ET,3),
            molecule.mesh.number_of_elements))
        csv_data.close()
    
    # Pausa o interrupcion de usuario
    except KeyboardInterrupt:
        print('Interrupcion')
        error_file=open('error_file.txt','w')
        error_file.close()
        break

    except OSError:
        #Algunas veces ocurren errores y mesh_temp queda ocupado, con esto se evita arrastrar el error a las moleculas siguientes
        if 'mesh_temp' in os.listdir():
            os.chdir('mesh_temp')
            for arch in os.listdir():
                os.remove(arch)
            os.chdir('..')
            os.rmdir('mesh_temp')
    except Exception:
        import sys
        exc_type,value,traceback=sys.exc_info()
        print(exc_type.__name__)
        print(traceback)

        #Algunas veces ocurren errores y mesh_temp queda ocupado, con esto se evita arrastrar el error a las moleculas siguientes
        if 'mesh_temp' in os.listdir():
            os.chdir('mesh_temp')
            for arch in os.listdir():
                os.remove(arch)
            os.chdir('..')
            os.rmdir('mesh_temp')
