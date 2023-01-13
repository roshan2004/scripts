import MDAnalysis as mda
import numpy as np
import argparse
import math
import subprocess


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_gro", type = str, help = 'input gro file')
parser.add_argument("-ci", "--water_gro", type = str, help = 'input water gro file')
parser.add_argument("-o", "--output_gro", type = str, help = "output gro file")

args = parser.parse_args()

input_gro = args.input_gro

output_gro = args.output_gro

water_gro = args.water_gro

command = "gmx insert-molecules"




u = mda.Universe(input_gro)


volume = u.atoms.ts.volume * 1e-30 # Angstrom^3 ----> m ^3

density = 1000 # kg/m^3

mass_one_water_molecule = (18 * 1e-3)/(6.023e23)

number_of_AA_water_mols = (density * volume) / mass_one_water_molecule

number_of_CG_water_mols = math.ceil(number_of_AA_water_mols/4)

args = str(f"-f {input_gro} -ci {water_gro} -o {output_gro} -nmol {str(number_of_CG_water_mols)}")





subprocess.run(["gmx","insert-molecules", "-f",input_gro,"-ci",water_gro,"-o", output_gro, "-nmol", str(number_of_CG_water_mols)], check = True)


