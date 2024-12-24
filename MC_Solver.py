# Librerias importantes
import os
import numpy 
import pbj
import sys
import glob
import argparse
from time import time
from datetime import datetime
import gc
from scipy.sparse import csr_matrix, SparseEfficiencyWarning
import warnings

def check_parser(argv):
    """ 
    Reads the inputs from the command line
    """

    parser = argparse.ArgumentParser(description="Run PBJ calculations for all pqr files in a folder")
    parser.add_argument('-f','--folder', dest='folder', type=str, default=None, help='Folder with pqr samples')
    parser.add_argument('-of','--output_file_name', dest='output_file_name', type=str, default=None, help='Output file name')
    parser.add_argument('-k','--kappa', dest='kappa', type=float, default=0.125, help='Inverse of Debye length, defaults to 0.125 angs^-1')
    parser.add_argument('-e1','--epsilon_in', dest='epsilon_in', type=float, default=4., help='Dielectric constant in molecule region. Defaults to 4.')
    parser.add_argument('-d','--mesh_density', dest='mesh_density', type=float, default=4., help='Vertices per square angs of surface mesh. Defaults to 4.')
    parser.add_argument('-mg','--mesh_generator', dest='mesh_generator', type=str, default="nanoshaper", help='Mesh generator, msms or nanoshaper.')
    parser.add_argument('-sub','--n_subset', dest='n_subset', type=int, default=None, help='Number of cases if only a subset will be run. Takes first n_subset cases')

    args = parser.parse_args(argv)

    if args.folder[-1] != "/":
        args.folder += "/"

    return args.folder, args.output_file_name, args.kappa, args.epsilon_in, args.mesh_density, args.n_subset

def generate_unique_file_name(file_name):

    if not os.path.exists(file_name):
        return file_name

    # Extract name and extension
    file_clean, file_extension = os.path.splitext(file_name)

    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

    new_file_name = f"{file_clean}_{timestamp}{file_extension}"

    return new_file_name

def run_mc(folder, output_file_name=None, kappa=0.125, epsilon_in=4., mesh_density=4, mesh_generator="nanoshaper", n_subset=None):

    warnings.simplefilter("ignore", SparseEfficiencyWarning)


    if output_file_name==None:
        output_file = folder + "output.csv"
    else:
        output_file = folder + output_file_name

    output_file = generate_unique_file_name(output_file)

    csv_data = open(output_file, "w")
    csv_data.write("iter,pqr_file,elec_energy,cav_energy,disp_energy,nonpolar_energy,solv_energy,elapsed_time,number_of_elements\n")

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
            simulation = pbj.implicit_solvent.Simulation()

            # Carga de archivo pqr
            molecule = pbj.Solute(
                test_file, mesh_generator=mesh_generator, mesh_density=mesh_density
            )

            # Se agrega soluto
            simulation.add_solute(molecule)

            # Parametros de la simulacion
            simulation.gmres_tolerance = 1e-5
            simulation.gmres_max_iterations = 400
            simulation.kappa = kappa

            simulation.solutes[0].ep_in = epsilon_in

            # Calculo de energia de solvatacion
            simulation.calculate_solvation_energy(electrostatic_energy=True, nonpolar_energy=True)

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
                "{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(
                    i,
                    test_file,
                    simulation.solutes[0].results["electrostatic_solvation_energy"],
                    simulation.solutes[0].results["cavity_energy"],
                    simulation.solutes[0].results["dispersion_energy"],
                    simulation.solutes[0].results["nonpolar_solvation_energy"],
                    simulation.solutes[0].results["solvation_energy"],
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

    folder, output_file_name, kappa, epsilon_in, mesh_density, mesh_generator, n_subset = check_parser(sys.argv[1:])
    run_mc(folder, output_file_name, kappa, epsilon_in, mesh_density, mesh_generator, n_subset)
