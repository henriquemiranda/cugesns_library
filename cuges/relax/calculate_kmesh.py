from pymatgen.io.vasp.inputs import Incar, Poscar, Kpoints, Potcar
from pymatgen.io.vasp.outputs import Incar, Poscar, Kpoints, Potcar
from pymatgen.io.vasp.sets import *

vasprun = Vasprun('vasprun.xml')
kpoints = Kpoints.automatic_density(vasprun.final_structure, 3000)
kpoints.write_file('KPOINTS')
poscar = Poscar(vasprun.final_structure)
poscar.write_file('POSCAR')
