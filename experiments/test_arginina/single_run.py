import pbj

# Se crea y limpia objeto simulation
simulation = pbj.implicit_solvent.Simulation()

# Carga de archivo pqr
molecule = pbj.Solute(
    test_file, mesh_generator="msms", mesh_density=mesh_density
)

# Se agrega soluto
simulation.add_solute(molecule)

# Parametros de la simulacion
simulation.gmres_tolerance = 1e-5
simulation.gmres_max_iterations = 400
simulation.kappa = 0

simulation.solutes[0].ep_in = 4

# Calculo de energia de solvatacion
simulation.calculate_solvation_energy(electrostatic_energy=True, nonpolar_energy=True)

# Impresion por pantalla

print(
    "INFO: {0}\ti:{1},\t {2}\t(kcal/mol)\t {3} [s],{4} elem.".format(
        output_file.split(".")[0],
        i,
        molecule.results["solvation_energy"],
        round(ET, 3),
        molecule.mesh.number_of_elements,
    )
)