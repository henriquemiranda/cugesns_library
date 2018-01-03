import sys
from phonopy.units import *
from phonopy import Phonopy
from phonopy.interface.phonopy_yaml import *
from phonopy.structure.atoms import Atoms
import phonopy.file_IO as file_IO
import matplotlib.pyplot as plt

options = {'all':    {'forces':True, 'masses':True },
           'masses': {'forces':False,'masses':True },
           'forces': {'forces':True, 'masses':False}}
emin=50
emax=450
energies_cm1 = np.linspace(emin,emax,500)

def get_dos(x,freqs,sigma=3):
    y = np.zeros(len(x))

    #lorentzian stuff
    s2 = (.5*sigma)**2
    c = (.5*sigma)

    for e in freqs:
        x1 = (x-e)**2
        y += c/(x1+s2)
    return x,y

def read_ingredients(path):
    """
    read the ingredients to calculate the phonons from the file
    """

    #read input
    ph_yaml = PhonopyYaml(calculator='vasp')
    ph_yaml.read('%s/phonopy_disp.yaml'%path)
    atoms            = ph_yaml._unitcell
    supercell_matrix = ph_yaml._data['supercell_matrix']

    #get forces
    force_constants = file_IO.parse_FORCE_CONSTANTS(filename='%s/FORCE_CONSTANTS'%path)

    return atoms, supercell_matrix, force_constants

def init_phonon(atoms, supercell_matrix, force_constants):
    """
    initialize phononpy class from input data
    """
    
    #create phonon structure
    phonon = Phonopy(atoms,supercell_matrix)
    phonon.set_force_constants(force_constants)

    return phonon

def read_phonon(path):
    return init_phonon(*read_ingredients(path))

def mix_atoms(at1,at2,ratio1):
    """
    return a structure that is a mix between a1 and a2
    at = at1*ratio1 + at2*(1-ratio1)
    """
    #at1
    cell1 = at1.get_cell() 
    pos1  = at1.get_scaled_positions()
 
    #at2
    cell2 = at2.get_cell() 
    pos2  = at2.get_scaled_positions()
  
    #mix
    cell = cell1*(1-ratio1) + cell1*ratio1
    pos  = pos1*(1-ratio1)  + pos2*ratio1

    numbers = at1.get_atomic_numbers()

    #at
    at = Atoms(numbers=numbers,scaled_positions=pos,cell=cell)
    return at
     
def mix_phonon(struct1,struct2,ratio1,forces=True,cell=True,masses=True):
    """
    mix the phonons from two structures according to 0 < ratio1 < 1
    the final result is struct1*ratio1 + struct2*(1-ratio1)

    a structure is: atoms, supercell_matrix, force sets
    """
    #unpack 
    at1, sc1, fc1 = struct1
    at2, sc2, fc2 = struct2

    #consistency check
    if sc1 != sc2:
        raise ValueError("Incompatible calculations")
    else:
        supercell_matrix = sc1

    #mix forces   
    if forces:
        force_constants = fc1*(1-ratio1) + fc2*ratio1
    else:
        force_constants = fc1  

    #mix structure 
    if cell:
        atoms = mix_atoms(at1,at2,ratio1)
    else:
        atoms = at1 

    #mix masses
    if masses:
        mass1 = at1.get_masses()
        mass2 = at2.get_masses()
        mass = mass1*(1-ratio1) + mass2*ratio1
        atoms.set_masses(mass)

    #create phonopy structure
    phonon = Phonopy(atoms,supercell_matrix)
    phonon.set_force_constants(force_constants)

    return phonon


cuges = read_ingredients('cuges')
cusns = read_ingredients('cusns')

steps = np.linspace(0,1,10)

def calculate_freqs(tag, s1=cusns, s2=cuges):
    phfreqs = []
    for x in steps:
        phonon = mix_phonon(s1, s2, x, **options[tag])
        freqs = phonon.get_frequencies((0,0,0))*THzToCm
        phfreqs.append(freqs)
        print("%5.3lf %12.8lf" % (x, freqs[-1]))
    return phfreqs

def plot_spectra(phfreqs,ax=plt,save=False):
    phfreqs = np.array(phfreqs)
    freq0 = phfreqs[0,-3]
    freqx = phfreqs[-1,-3]

    wmin = phfreqs[-1,-3]
    wmax = phfreqs[ 0,-3]
    ax.axvline(wmin,ls='--',lw=1)
    ax.axvline(wmax,ls='--',lw=1)
    ax.text(wmax,-.05,"%6.2lf"%wmax)
    ax.text(wmin,1.10,"%6.2lf"%wmin)
    print(wmax-wmin)

    for i,s in enumerate(steps):
        e,dos = get_dos(energies_cm1,phfreqs[i])
        ax.plot(e,dos*0.1+s)
    ax.set_xlim(emin,emax)
    ax.set_xlabel('frequencies cm$^{-1}$')
    if save: plt.savefig('cuges_hse06_%s_spectra.pdf'%tag)
    return ax

def plot_frequencies(phfreqs,ax=plt,save=False):
    phfreqs = np.array(phfreqs)
    freq0 = phfreqs[0,-3]
    freqx = phfreqs[-1,-3]

    wmin = phfreqs[-1,-3]
    wmax = phfreqs[ 0,-3]
    ax.axhline(wmin,ls='--',lw=1,c='r')
    ax.axhline(wmax,ls='--',lw=1,c='r')
    ax.text(.01,wmax+2,"%6.2lf"%wmax)
    ax.text(.78,wmin+2,"%6.2lf"%wmin)
    print(wmax-wmin)

    for phfreq in phfreqs:
        ax.plot(steps,phfreqs,c='C0',alpha=0.3)

    ax.set_xlim(min(steps),max(steps))
    ax.set_ylim(200,450)
    ax.set_xlabel('x')
    if save: plt.savefig('cuges_hse06_%s_freqs.pdf'%tag)
    return ax

size=3   
fig = plt.figure(figsize=(size*3,size+2))
#fig.suptitle("CuGeS to Cu(Sn$_{1-x}$Ge$_x$)S HSE06")

ax = plt.subplot(1,3,1)
ax.set_title('a) interpolation of forces\n and masses')
ax.set_ylabel('frequencies cm$^{-1}$')
phfreqs = calculate_freqs('all')
plot_frequencies(phfreqs,ax)

ax = plt.subplot(1,3,2)
ax.set_title('b) interpolation of forces only\n using the masses of CTS')
#ax.set_ylabel('frequencies cm$^{-1}$')
#ax.get_yaxis().set_visible(False)
phfreqs = calculate_freqs('forces')
plot_frequencies(phfreqs,ax)

ax = plt.subplot(1,3,3)
ax.set_title('c) interpolation of masses only\n using the forces of CGS')
#ax.set_ylabel('frequencies cm$^{-1}$')
#ax.get_yaxis().set_visible(False)
phfreqs = calculate_freqs('masses', s1=cuges, s2=cusns)
phfreqs.reverse()
plot_frequencies(phfreqs,ax)
plt.subplots_adjust(left=0.1, right=0.95,wspace=0.2,bottom=0.15,top=0.8)
plt.savefig('freqs_3colum.pdf')
plt.show()
