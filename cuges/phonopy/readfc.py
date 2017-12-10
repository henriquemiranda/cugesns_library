#
# 1 February 2016
#
# Obtain phonondb structures  from
# http://phonondb.mtl.kyoto-u.ac.jp
# and plot the pohnon dispersion on the path obtained from the materials project API
#
from phonopy.units import *
from phonopy import Phonopy
from phonopy.interface.phonopy_yaml import *
from phonopy.structure.atoms import Atoms
import phonopy.file_IO as file_IO

#read input
ph_yaml = PhonopyYaml(calculator='vasp')
ph_yaml.read('phonopy_disp.yaml')
atoms            = ph_yaml._unitcell
supercell_matrix = ph_yaml._data['supercell_matrix']
phonon = Phonopy(atoms,supercell_matrix)

if 1:
    #get forces
    force_sets = file_IO.parse_FORCE_SETS(filename='FORCE_SETS')
    phonon.set_displacement_dataset(force_sets)

    #get force constants
    phonon.produce_force_constants()
    phonon.symmetrize_force_constants_by_space_group()

    #write force constants
    force_constants = phonon._force_constants
    file_IO.write_FORCE_CONSTANTS(force_constants,filename='FORCE_CONSTANTS')

#get forces
force_sets = file_IO.parse_FORCE_CONSTANTS(filename='FORCE_CONSTANTS')
phonon.set_force_constants(force_sets)

#get frequencies at gamma
for f in phonon.get_frequencies((0,0,0))*THzToCm:
    print "%12.8lf" % f
 
