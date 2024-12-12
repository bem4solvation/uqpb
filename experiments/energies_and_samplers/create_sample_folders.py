"""
This script contains auxiliary functions and execution to run the script sampler.py from the console

python sampler.py -nt N_TESTS -nw N_WORKERS -pqr BASE_PQR_FILE -f FOLDER

"""

# To run the script sampler.py from the console and handle directories
import os

# To set a timer between samples and time the execution
import time


def run_sampler(n_test, n_workers, pqr_file, folder, sampler):
    """Runs the sampler.py script from the console with the given parameters"""
    if sampler != "":
        os.system(
            "python sampler.py -nt {} -nw {} -pqr {} -f {} -sp {}".format(
                n_test, n_workers, pqr_file, folder, sampler
            )
        )
    else: 
        os.system(
            "python sampler.py -nt {} -nw {} -pqr {} -f {}".format(
                n_test, n_workers, pqr_file, folder
            )
        )

def create_folder_name(pqr_file,sampler):
    """Creates a folder name from the pqr file name and sampler"""
    if sampler=="":
        sampler = "pseudo"
    folder_name = "output_{pqr_file}_{sampler}".format(pqr_file=pqr_file.split("/")[-1].split(".")[0],sampler=sampler)
    return folder_name

if __name__ == "__main__":

    sampler_list = ["", "LHS", "Halton", "Hammersley", "Sobol"]
    n_test = 1000
    n_workers = 2
    pqr_file = "mobley_test_pqr/1112_tetrachloroethane.pqr"
    folder = "tests"
    for sampler in sampler_list:
        folder_name = folder + "/" + create_folder_name(pqr_file,sampler)
        run_sampler(n_test, n_workers, pqr_file, folder_name, sampler)
        print("Creating folder: ", folder_name)
        time.sleep(10)