# Librerias importantes
import os
import numpy 
import pbj
import sys
import glob
import argparse
from time import time
import gc
from scipy.sparse import csr_matrix, SparseEfficiencyWarning
import warnings

def check_parser(argv):
    """ 
    Reads the inputs from the command line
    """

    parser = argparse.ArgumentParser(description="generate random variables for montecarlo")
    parser.add_argument('-f','--folder', dest='folder', type=str, default=None, help='Folder with pqr samples')
    parser.add_argument('-k','--kappa', dest='kappa', type=float, default=0.125, help='Inverse of Debye length, defaults to 0.125 angs^-1')
    parser.add_argument('-e1','--epsilon_in', dest='epsilon_in', type=float, default=4., help='Dielectric constant in molecule region. Defaults to 4.')
    parser.add_argument('-d','--mesh_density', dest='mesh_density', type=float, default=4., help='Vertices per square angs of surface mesh. Defaults to 4.')
    parser.add_argument('-sub','--n_subset', dest='n_subset', type=int, default=None, help='Number of cases if only a subset will be run. Takes first n_subset cases')

    args = parser.parse_args(argv)

    if args.folder[-1] != "/":
        args.folder += "/"

    return args.folder, args.kappa, args.epsilon_in, args.mesh_density, args.n_subset


def run_mc(folder, kappa=0.125, epsilon_in=4., mesh_density=4, n_subset=None):

    warnings.simplefilter("ignore", SparseEfficiencyWarning)

    output_file = folder + "output.csv"

    csv_data = open(output_file, "w")
    csv_data.write("ITERATION,SOLV. ENERGY,Elapsed time,Number of elements\n")

    pqr_files = glob.glob(folder + "*.pqr")

    if n_subset == None:
        n_cases = len(pqr_files)
    else:
        n_cases = n_subset


    for i in range(n_cases):

        try:
            start_time = time()
            
            test_file = pqr_files[i]
            
            # Se crea y limpia objeto simulation
            simulation = pbj.electrostatics.Simulation()

            # Carga de archivo pqr
            molecule = pbj.Solute(
                test_file, mesh_generator="msms", mesh_density=mesh_density
            )

            # Se agrega soluto
            simulation.add_solute(molecule)

            # Parametros de la simulacion
            simulation.gmres_tolerance = 1e-5
            simulation.gmres_max_iterations = 400
            simulation.kappa = kappa

            simulation.solutes[0].ep_in = epsilon_in

            # Calculo de energia de solvatacion electrostatica
            simulation.calculate_solvation_energy()

            # Impresion por pantalla
            ET = time() - start_time
            print(
                "INFO: {0}\ti:{1},\t {2}\t(kcal/mol)\t {3} [s],{4} elem.".format(
                    output_file.split(".")[0],
                    i,
                    molecule.results["solvation_energy"],
                    round(ET, 3),
                    molecule.mesh.number_of_elements,
                )
            )

            # Escritura en archivo de salida
            csv_data.write(
                "{0},{1},{2},{3}\n".format(
                    i,
                    molecule.results["solvation_energy"],
                    round(ET, 3),
                    molecule.mesh.number_of_elements,
                )
            )

        # Pausa o interrupcion de usuario
        except KeyboardInterrupt:
            print("Interrupcion")
            error_file = open("error_file.txt", "w")
            error_file.close()
            break

        except OSError:
            # Algunas veces ocurren errores y mesh_temp queda ocupado, con esto se evita arrastrar el error a las moleculas siguientes
            if "mesh_temp" in os.listdir():
                os.chdir("mesh_temp")
                for arch in os.listdir():
                    os.remove(arch)
                os.chdir("..")
                os.rmdir("mesh_temp")

        except Exception:
            import sys

            exc_type, value, traceback = sys.exc_info()
            print(exc_type.__name__)
            print(traceback)

            # Algunas veces ocurren errores y mesh_temp queda ocupado, con esto se evita arrastrar el error a las moleculas siguientes
            if "mesh_temp" in os.listdir():
                os.chdir("mesh_temp")
                for arch in os.listdir():
                    os.remove(arch)
                os.chdir("..")
                os.rmdir("mesh_temp")

    csv_data.close()

if __name__ == "__main__":

    folder, kappa, epsilon_in, mesh_density, n_subset = check_parser(sys.argv[1:])
    run_mc(folder, kappa, epsilon_in, mesh_density, n_subset)
