import argparse
import numpy
import sys
import os
import pbj
import numpy as np
import skopt

from utilities import masas, amino_acids


def check_parser(argv):
    """
    Reads the inputs from the command line
    """

    parser = argparse.ArgumentParser(
        description="generate random variables for montecarlo"
    )

    parser.add_argument(
        "-nt", "--n_test", dest="n_test", type=int, default=4, help="Number of MC tests"
    )
    parser.add_argument(
        "-nw",
        "--n_workers",
        dest="n_workers",
        type=int,
        default=1,
        help="Number of folders to store configurations",
    )
    parser.add_argument(
        "-pqr",
        "--pqr_file",
        dest="pqr_file",
        type=str,
        default="mobley_test_pqr/1112_tetrachloroethane.pqr",
        help="Base pqr file name",
    )
    parser.add_argument(
        "-f",
        "--folder",
        dest="folder",
        type=str,
        default="tests",
        help="Working folder, defaults to name of pqr file",
    )
    parser.add_argument(
        "-sp",
        "--sampler",
        dest="sampler",
        type=str,
        default="pseudo",
        help="Sampler used. Can be pseudorandom or quasirandom samples.",
    )
    parser.add_argument(
        "-pd",
        "--prob_dist",
        dest="prob_dist",
        type=str,
        default="uniform",
        help="Probability distribution in shake. Uniform (default) or normal. Normal is the half normal distribution for the thermal radius.",
    )
    parser.add_argument(
        "-sr",
        "--shake_radius",
        dest="shake_radius",
        type=float,
        default=None,
        help="Shake radius. Defaults to None, which uses the thermal radius.",
    )
    parser.add_argument(
        "-tt",
        "--t_thermal",
        dest="t_thermal",
        type=float,
        default=1e-8,
        help="Characteristic thermal time for thermal radius. Defaults to 1e-8.",
    )

    args = parser.parse_args(argv)

    if args.pqr_file[-4:] != ".pqr":
        args.pqr_file += ".pqr"

    index_folder = args.pqr_file.rfind("/")
    pqr_file_name = args.pqr_file[index_folder:]
    pqr_dir = args.pqr_file[:index_folder]

    if args.folder is None:
        args.folder = pqr_file_name[:-4]

    if args.folder[-1] != "/":
        args.folder += "/"

    return (
        pqr_dir,
        pqr_file_name,
        args.folder,
        args.n_test,
        args.n_workers,
        args.sampler,
        args.prob_dist,
        args.shake_radius,
        args.t_thermal,
    )


def count_pqr_atoms(pqr_file):
    """
    Counts the number of atoms in pqr_file
    """

    pqr = open(pqr_file, "r")
    pqr_data = pqr.readlines()
    pqr.close()

    n_atom = 0
    for line in pqr_data:
        if line[:4] == "ATOM":
            n_atom += 1

    return n_atom


def atom_name_fix(atom, is_biomolecule=False):
    """
    Tomando un string correspondiente al átomo de un archivo pdb (formato C1, C2, etc)
    se retorna sólo el elemento químico sin el número.
    En caso de que sea biomolecula, puede que algunos atomos se nombren en formato CA (carbono alpha)
    y derivados, en ese caso se retorna solo CHONSP"""

    final = ""
    for caracter in atom:
        if caracter not in "0123456789":
            final += caracter
    if is_biomolecule:
        for elem in "CHONSP":
            if elem in final:
                return elem
    return final


def sample(n_samples, dimension, sampler="pseudo"):
    """Generate pseudorandom or quasirandom samples in [0, 1]^dimension.
    Reference: DeepXDE https://github.com/lululxvi/deepxde.git

    Args:
        n_samples (int): The number of samples.
        dimension (int): Space dimension.
        sampler (string): One of the following: "pseudo" (pseudorandom), "LHS" (Latin
            hypercube sampling), "Halton" (Halton sequence), "Hammersley" (Hammersley
            sequence), or "Sobol" (Sobol sequence).
    """
    if sampler == "pseudo":
        return pseudorandom(n_samples, dimension)
    if sampler in ["LHS", "Halton", "Hammersley", "Sobol"]:
        return quasirandom(n_samples, dimension, sampler)
    raise ValueError("f{sampler} sampling is not available.")


def pseudorandom(n_samples, dimension):
    """Pseudo random."""
    # If random seed is set, then the rng based code always returns the same random
    # number, which may not be what we expect.
    # rng = np.random.default_rng(config.random_seed)
    # return rng.random(size=(n_samples, dimension), dtype=config.real(np))
    return np.random.random(size=(n_samples, dimension))


def quasirandom(n_samples, dimension, sampler):
    # Certain points should be removed:
    # - Boundary points such as [..., 0, ...]
    # - Special points [0, 0, 0, ...] and [0.5, 0.5, 0.5, ...], which cause error in
    #   Hypersphere.random_points() and Hypersphere.random_boundary_points()
    skip = 0
    if sampler == "LHS":
        sampler = skopt.sampler.Lhs()
    elif sampler == "Halton":
        # 1st point: [0, 0, ...]
        sampler = skopt.sampler.Halton(min_skip=1, max_skip=1)
    elif sampler == "Hammersley":
        # 1st point: [0, 0, ...]
        if dimension == 1:
            sampler = skopt.sampler.Hammersly(min_skip=1, max_skip=1)
        else:
            sampler = skopt.sampler.Hammersly()
            skip = 1
    elif sampler == "Sobol":
        # 1st point: [0, 0, ...], 2nd point: [0.5, 0.5, ...]
        sampler = skopt.sampler.Sobol(randomize=False)
        if dimension < 3:
            skip = 1
        else:
            skip = 2
    space = [(0.0, 1.0)] * dimension
    return np.asarray(sampler.generate(space, n_samples + skip)[skip:])


def average_thermal_length(atom_name, res_name, t_thermal):
    """
    Calculate the average thermal length for each atom
    """
    # Average Thermal Length: definicion de funciones y creacion de diccionario
    ATL_dic = {}
    kB = 1.3806e-23  # Constante de Stefan-Boltzmann
    T = 298  # Temperatura en Kelvin
    vt = lambda m, T: numpy.sqrt(
        kB * t_thermal / m
    )  # thermal velocity en funcion de masa y temperatura (298K por defecto)
    L = lambda m, T: vt(m, T) * t_thermal * 1e10  # 1e10 convierte la unidad a Angstroms
    for atom, mass in masas.items():
        ATL_dic[atom] = L(mass, T)

    r_thermal = numpy.zeros(len(atom_name))
    for i, atom_raw in enumerate(atom_name):
        if res_name[i] in amino_acids:
            is_biomolecule = True
        else:
            is_biomolecule = False

        atom = atom_name_fix(atom_raw, is_biomolecule)

        if atom not in ATL_dic:
            print("Atom " + atom + " not found")
            return
        else:
            r_thermal[i] = ATL_dic[atom]

    return r_thermal


def write_new_pqr(pqr_file, x_new, out_pqr_name):
    """
    Writes pqr file with new (shaken) positions
    """
    old_pqr = open(pqr_file, "r")
    old_pqr_data = old_pqr.readlines()
    old_pqr.close()

    new_pqr = open(out_pqr_name, "w")

    i = 0
    for line in old_pqr_data:
        if line[:4] != "ATOM":
            new_pqr.write(line)

        else:
            line_space = line.replace("-", " -")

            words = line_space.split()

            if len(words) == 11:
                x_old = words[6]
                y_old = words[7]
                z_old = words[8]

            elif len(words) == 10:
                x_old = words[5]
                y_old = words[6]
                z_old = words[7]

            line = line.replace(x_old, "%1.4f" % (x_new[i, 0]))
            line = line.replace(y_old, "%1.4f" % (x_new[i, 1]))
            line = line.replace(z_old, "%1.4f" % (x_new[i, 2]))

            new_pqr.write(line)

            i += 1

    new_pqr.close()


def read_pqr(pqr_file):
    """
    Read in teh x,y,z positions in the pqr
    """
    pqr = open(pqr_file, "r")
    pqr_data = pqr.readlines()
    pqr.close()

    n_atom = 0
    for line in pqr_data:
        if line.startswith("ATOM"):
            n_atom += 1

    i = 0
    x_atom = numpy.zeros((n_atom, 3), dtype=float)
    atom_name = numpy.empty((n_atom,), dtype=object)
    res_name = numpy.empty((n_atom,), dtype=object)
    for line in pqr_data:
        if line.startswith("ATOM"):
            line_space = line.replace("-", " -")

            words = line_space.split()

            if len(words) == 11:
                x = words[6]
                y = words[7]
                z = words[8]

            elif len(words) == 10:
                x = words[5]
                y = words[6]
                z = words[7]

            x_atom[i, 0] = float(x)
            x_atom[i, 1] = float(y)
            x_atom[i, 2] = float(z)

            atom_name[i] = words[2]
            res_name[i] = words[3]

            i += 1

    return x_atom, atom_name, res_name


def generate_random_samples(
    n_test,
    n_workers,
    n_atom,
    folder,
    pqr_dir,
    pqr_file_name,
    sampler="pseudo",
    prob_dist="uniform",
    shake_radius=None,
    t_thermal=1e-8,
):
    """
    Generates random coefficients and modified pqr's
    """

    digits = len(str(n_test))
    if not os.path.exists(folder):
        os.makedirs(folder)

    if n_workers > n_test:
        print("n_workers has to be smaller than n_test, defaulting to n_workers=1")
        n_workers = 1

    jobs_per_worker = n_test // n_workers
    digits_workers = len(str(n_workers))
    for i in range(n_workers):
        worker_folder = "job_" + str(i).zfill(digits_workers) + "/"
        if not os.path.exists(folder + worker_folder):
            os.makedirs(folder + worker_folder)

    pqr_file = pqr_dir + pqr_file_name

    x_base, atom_name, res_name = read_pqr(pqr_file)

    if shake_radius == None:
        r_thermal = average_thermal_length(atom_name, res_name, t_thermal)

    if prob_dist == "uniform":
        coeffs = sample(n_test, n_atom * 3, sampler=sampler)

    for i in range(n_test):
        if shake_radius == None:
            if prob_dist == "uniform":
                shake = coeffs[i].reshape((n_atom, 3))
                # shake = numpy.random.uniform(low = 0.0, high = 1.0, size = (n_atom,3))
            elif prob_dist == "normal":
                shake_rad = numpy.absolute(
                    numpy.random.normal(loc=0.0, scale=1.0, size=(n_atom))
                )
                shake_angle = numpy.random.uniform(low=0.0, high=1.0, size=(n_atom, 2))

                shake = numpy.zeros((n_atom, 3))
                shake[:, 0] = shake_rad
                shake[:, 1:] = shake_angle[:, :]

            else:
                print("Probability distribution not understood")
                return

        else:
            r_thermal = shake_radius

        worker_number = i // jobs_per_worker
        if worker_number >= n_workers:
            worker_number = i % n_workers

        worker_folder = "job_" + str(worker_number).zfill(digits_workers) + "/"

        out_file_name = (
            folder
            + worker_folder
            + pqr_file_name[:-4]
            + "_"
            + str(i).zfill(digits)
            + "_coeff.txt"
        )
        out_pqr_name = (
            folder
            + worker_folder
            + pqr_file_name[:-4]
            + "_"
            + str(i).zfill(digits)
            + ".pqr"
        )

        r = shake[:, 0] * r_thermal
        theta = shake[:, 1] * 2 * numpy.pi
        phi = shake[:, 2] * numpy.pi

        x_new = numpy.zeros((n_atom, 3))
        x_new[:, 0] = x_base[:, 0] + r * numpy.cos(theta) * numpy.sin(phi)
        x_new[:, 1] = x_base[:, 1] + r * numpy.sin(theta) * numpy.sin(phi)
        x_new[:, 2] = x_base[:, 2] + r * numpy.cos(phi)

        write_new_pqr(pqr_file, x_new, out_pqr_name)

        numpy.savetxt(out_file_name, shake)


if __name__ == "__main__":
    (
        pqr_dir,
        pqr_file_name,
        folder,
        n_test,
        n_workers,
        sampler,
        prob_dist,
        shake_radius,
        t_thermal,
    ) = check_parser(sys.argv[1:])
    n_atom = count_pqr_atoms(pqr_dir + pqr_file_name)
    generate_random_samples(
        n_test,
        n_workers,
        n_atom,
        folder,
        pqr_dir,
        pqr_file_name,
        sampler,
        prob_dist,
        shake_radius,
        t_thermal,
    )
