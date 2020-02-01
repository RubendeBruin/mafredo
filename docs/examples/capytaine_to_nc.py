"""
This is a basic example of how to generate a an output file from capytaine.

It is copy-pasted together from the rao and animation examples from the capytaine package and is provided here
for convenience.

"""


import logging

import numpy as np
from numpy import pi

import capytaine as cpt

logging.basicConfig(level=logging.INFO, format='%(levelname)-8s: %(message)s')
bem_solver = cpt.BEMSolver()


def make_database(body, omegas, wave_directions):
    # SOLVE BEM PROBLEMS
    problems = []
    for wave_direction in wave_directions:
        for omega in omegas:
            problems += [cpt.RadiationProblem(omega=omega, body=body, radiating_dof=dof) for dof in body.dofs]
            problems += [cpt.DiffractionProblem(omega=omega, body=body, wave_direction=wave_direction)]
    results = [bem_solver.solve(problem) for problem in problems]
    *radiation_results, diffraction_result = results
    dataset = cpt.assemble_dataset(results)

    dataset['diffraction_result'] = diffraction_result

    return dataset


def generate_boat() -> cpt.FloatingBody:
    boat = cpt.RectangularParallelepiped(size=(100,30,8),resolution=(40, 12, 8))
    boat.add_all_rigid_body_dofs()
    boat.keep_immersed_part()

    return boat

def generate_mesh() -> cpt.FloatingBody:
    boat = cpt.FloatingBody.from_file(r"C:\data\docker\output.stl", file_format="stl", name="cheetah")

    # boat.clip(Plane(normal=(0, 0, 1), point=(0, 0, -0.01)))  # <-- optional, depending on the input mesh
    boat.add_all_rigid_body_dofs()
    boat.keep_immersed_part()

    return boat



if __name__ == '__main__':

    omega = [0.01, 0.02, 0.04,0.06,0.08,0.1,0.15, 0.2,0.25, 0.3,0.35, 0.4,0.45, 0.5,0.55, 0.6,0.65, 0.7, 0.8,0.9, 1.0,1.1,1.2,1.4,1.6,1.8,2.0,4.0]

    body = generate_mesh()
    # body = generate_boat()

    body.show()

    dataset = make_database(body=body, omegas=omega, wave_directions=np.linspace(0,pi,9))
    print(dataset)

    from capytaine.io.xarray import separate_complex_values
    sep = separate_complex_values(dataset)
    sep.to_netcdf(r"C:\data\python\rao\docs\examples\cheetah.nc",
                  encoding={'radiating_dof': {'dtype': 'U'},
                            'influenced_dof': {'dtype': 'U'},
                            'diffraction_result': {'dtype': 'U'}})
