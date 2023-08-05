import random
import numpy as np
import os


# Thermal functions
def new_name(main, i, n):
    """Crea un nombre en formato string para un archivo nuevo correspondiente
    a las N pruebas del análisis de Montecarlo"""
    i = str(i)
    num = i
    if len(i) < len(str(n)):
        dif = len(str(n)) - len(i)
        num = "0" * dif + i
    return main.split(".")[0] + num + "." + main.split(".")[1]


def randomX(x, R):  # R indica el radio atómico en Angstroms
    """Teniendo un vector x, se agita una esfera que lo contiene entregando una
    nueva posición"""
    r = random.uniform(0, 1) * R
    theta = random.uniform(0, 1) * 2 * np.pi
    phi = random.uniform(0, 1) * np.pi
    new_x = [0, 0, 0]
    new_x[0] = x[0] + r * np.cos(theta) * np.sin(phi)
    new_x[1] = x[1] + r * np.sin(theta) * np.sin(phi)
    new_x[2] = x[2] + r * np.cos(phi)
    return new_x


def dec3(x):
    """Redondeo de un número float a 3 decimales para correcta escritura en .pqr"""
    valor = str(round(x, 3))
    if "." in valor:
        decimal = valor.split(".")[-1]
        if len(decimal) < 3:
            decimal = decimal + "000"
            decimal = decimal[:3]
            return valor.split(".")[0] + "." + decimal
        else:
            return valor
    else:
        return valor + ".000"


def nombre_atomo(atomo, es_biomolecula=False):
    """Tomando un string correspondiente al átomo de un archivo pdb (formato C1, C2, etc)
    se retorna sólo el elemento químico sin el número.
    En caso de que sea biomolecula, puede que algunos atomos se nombren en formato CA (carbono alpha)
    y derivados, en ese caso se retorna solo CHONSP"""
    final = ""
    for caracter in atomo:
        if caracter not in "0123456789":
            final += caracter
    if es_biomolecula:
        for elem in "CHONSP":
            if elem in final:
                return elem
    return final


def leer_archivo_entrada(prm_file):
    """Esta funcion lee un archivo de parametros de entrada para el solver y retorna un diccionario, de llaves predefinidas.
    El formato de cada linea del archivo de paramatros de entrada `input_config.prm` es PARAMETRO=VALOR
    """
    # Diccionario de parametros por defecto
    parametros = {
        "NAME": "",
        "TESTPATH": "",
        "PQRFILE": "",
        "N_TESTS": 100,
        "TIEMPO": 0,
        "SIGMA_MOBLEY": 0,
        "EP_IN": 1,
        "KAPPA": 0,
        "DENSIDAD": 3,
        "OUTPUT_FILE": "",
    }

    # A continuacion se lee el archivo de configuracion de entrada para actualizar diccionario
    arch = open(prm_file)
    for line in arch:
        if "=" in line:
            parametro, valor = line.strip().split("=")
            # Se separa transformacion de valores dependiendo del tipo de dato
            if parametro in ["TIEMPO", "EP_IN", "KAPPA", "DENSIDAD"]:
                parametros[parametro] = float(valor)
            elif parametro in ["N_TESTS"]:
                parametros[parametro] = int(valor)
            else:
                # Para valores de tipo strings
                parametros[parametro] = valor
    arch.close()
    return parametros


def shake_file(
    mainfile, i, path, remarks, end, nombre_molecula, atoms, n_tests, ATL_dic
):
    # Si no existe el path se crea
    if not os.path.exists(path):
        os.makedirs(path)
    file = open(os.path.join(path, new_name(nombre_molecula + ".pqr", i, n_tests)), "w")
    file.write(remarks)
    list_atoms = atoms.split("\n")
    for i in range(len(list_atoms)):
        linea = list_atoms[i].split()
        if len(linea) == 0:
            continue
        atomo = nombre_atomo(linea[2], True)  # 2 para pqr, -1 para pdb
        _x = list(map(float, linea[5:8]))  # 5:8 pqr, 6:9 pdb
        _x = randomX(_x, ATL_dic[atomo])
        x, y, z = list(map(str, _x))
        lineareal = list(list_atoms[i])
        extension = mainfile.split(".")[-1]
        # Si es pdb
        if extension == "pdb":
            c = 0
            for j in z[::-1] + " ":
                lineareal[53 - c] = j
                c += 1
            c = 0
            for j in y[::-1] + " ":
                lineareal[45 - c] = j
                c += 1
            c = 0
            for j in x[::-1] + " ":
                lineareal[37 - c] = j
                c += 1
            file.write("".join(lineareal) + "\n")
        # Si es pqr
        elif extension == "pqr":
            linea[5:8] = list(map(str, _x))
            file.write("\t".join(linea) + "\n")
    file.write(end)
    file.close()
    return None


shake_file(
    "tests/tetrafluoromethane/tetrafluoromethane0.pqr",
    0,
    "tests\\tetrafluoromethane",
    "",
    "",
    "tetrafluoromethane",
    "",
    5,
    {},
)
