import sys
sys.path.insert(1,'/global/project/projectdirs/metatlas/anaconda/lib/python2.7/site-packages' )
from metatlas import metatlas_objects as metob

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Descriptors

def ppm_error(mass, theoretical_mass):
    """
    returns the ppm error of a given observed mass, and theoretical mass
    """
    theoretical_mass = Descriptors.ExactMolWt(mol)
    ppm = (mass - theoretical_mass) / theoretical_mass * 1e6
    return ppm

def ppm_window(mass, ppm=5, result='bounds'):
    """
    given a mass and a ppm error, returns lower and upper bounds corresponding to that ppm window
    result can be "bounds" or "error"
    bounds is a list of lower and upper bounds
    error is an amu value of the error, so you can add or subtract the error as you wish
    """
    error = ppm/1e6 * mass
    lower_bound = mass - error
    upper_bound = mass + error
    if result.lower() == 'bounds':
        return [lower_bound, upper_bound]
    elif result.lower() == 'error':
        return error
    else:
        raise RuntimeError('%s is not a valid result argument' %(result))
