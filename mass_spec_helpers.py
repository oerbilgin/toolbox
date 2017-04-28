import sys
sys.path.insert(1,
    '/global/project/projectdirs/metatlas/anaconda/lib/python2.7/site-packages'
    )
from metatlas import metatlas_objects as metob

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Descriptors

def ppm_error(mass, theoretical_mass):
    """
    Returns the ppm error of a given observed mass, and theoretical mass
    """
    if not isinstance(theoretical_mass, float):
        if isinstance(theoretical_mass, int):
            theoretical_mass = float(theoretical_mass)
        else:
            theoretical_mass = Descriptors.ExactMolWt(theoretical_mass)
    ppm = (mass - theoretical_mass) / theoretical_mass * 1e6
    return abs(ppm)

def ppm_window(mass, ppm=5, result='bounds'):
    """
    Given a mass and a ppm error, returns lower and upper bounds 
    corresponding to that ppm window. 

    Inputs
    ------
    mass:   monoisotopic mass
    ppm:    the ppm error to construct the window
    result: "bounds" returns a list of lower and upper bounds
            "error" returns the amu value of the error, if you wanted to
            add and/or subtract the error as you wish

    Outputs
    -------
    Either a list of lower and upper mass bounds, or a single value
    corresponding to the amu error
    """
    error = ppm/1e6 * mass
    lower_bound = mass - error
    upper_bound = mass + error
    if result.lower() == 'bounds':
        return [lower_bound, upper_bound]
    elif result.lower() == 'error':
        return error
    else:
        raise RuntimeError(
            '%s is not a valid result argument' %(result)
            )

def accurate_mass_match(mass, compound_df=None, ppm=5, extract='inchi_key'):
    """
    Accurate mass searching against a compound database.

    Inputs
    ------
    mass:           An accurate monoisotopic mass
    compound_df:    A dataframe of compounds. Must have a column named
                    "mono_isotopic_molecular_weight"
    ppm:            ppm error to allow the mass matcher
    extract:        What compound information to return. Must correspond
                    to a valid column in compound_df
    
    Outputs
    -------
    cpd:            List of compounds that were matched
    """

    err = ppm_window(mass, ppm=ppm, result='error')

    potential_compounds = compound_df[
                            abs(compound_df['mono_isotopic_molecular_weight'] \
                            - mass) <= err]
    cpds = potential_compounds[extract].values.tolist()
    if len(cpds) == 0:
        cpds = None
    return cpds

def find_mass(target_weight, ppm, df, 
    mz_column_name='mz', rt_column_name='rt_peak', rt=None, rt_bound=0.05):
    """
    Given a target mass and ppm window, return the rows of a dataframe 
    corresponding to the mass. 
    Optionally enter a retention time and bound
    
    Inputs
    ------
    target_weight:  the m/z you are looking for
    ppm:            the ppm error you are willing to have
    df:             a dataframe with rows being individual features
    mz_column_name: the column name that has the feature's m/z, metatlas
                    default: "mz"
    rt_column_name: the column name that has the reature's retention 
                    time peak, metatlas default: "rt_peak
    rt:             the retention time for the feature
    rt_bound:       the upper and lower bound you're willing to search 
                    for the retention time, default +/- 0.05 minutes
    
    Outputs
    -------
    A slice of the input dataframe that only contains rows that matched 
    the search parameters
    """
    w = ppm_window(target_weight, ppm)
    a = df[(df[mz_column_name] >= w[0]) & (df[mz_column_name] <= w[1])]
    if rt is not None:
        a = a[(a[rt_column_name] > rt - rt_bound) & \
              (a[rt_column_name] < rt + rt_bound)]
    return a

def data_convert(row, data_pattern='Peak height', search='endswith'):
    """
    Converts an input Series into a dataframe with 'x' and 'y' columns, 
    with 'x' corresponding to run order, and 'y' corresponding to a 
    measurement of interest.

    Run number is determined by splitting by the word "Run" and then by 
    a period. For example: "experimentinfo_Run1.mzml"

    If your filenames don't look like that, then the run number finder 
    probably won't work
    
    Inputs
    ------
    row:            a Pandas Series object, where the colums contain run
                    information
    data_pattern:   is the text pattern of the column that contains data
                    you want to plot, default is "Peak height"
    search:         the Pandas string methods search type:
                    - if search is "endswith" or "startswith" then it 
                    will be a str.endswith or str.startswith search; 
                    multiple data_patterns to search should be in a 
                    TUPLE
                    - if search is "contains" then it will be a 
                    str.contains search; multiple data_patterns to 
                    search should be in one string separated by '|'
        
    Outputs
    -------
    A dataframe with columns 'x' and 'y' that is sorted in ascending 
    order of 'x'
    """
    if not isinstance(row, pd.Series):
        raise RuntimeError('input was not a Series datatype!')
    if search.lower()=='endswith':
        c = row.keys()[row.keys().str.endswith(data_pattern)]
    elif search.lower()=='startswith':
        c = row.keys()[row.keys().str.startswith(data_pattern)]
    elif search.lower()=='contains':
        c = row.keys()[row.keys().str.contains(data_pattern)]
    else:
        raise RuntimeError(
            'Did not recognize %s as a valid Pandas search pattern, or it is \
            not encompassed yet' %(search))
        
    y = row[c].values.tolist()
    x = []
    for cc in c:
        x.append(int(cc.split('Run')[1].split('.')[0]))
    t = pd.DataFrame({'x':x, 'y':y})
    t.sort_values('x', inplace=True)
    
    return t

def ctrl_expt_compare(medium, experimental, fold_change=4):
    """
    Removes features that are within a window of detection for both 
    medium control and an experimental condition.
    Input dataframes are data_convert() outputs! (i.e. 'x' and 'y' 
    columns corresponding to run number and ion intensity, respectively)
    Compares an experimental data point with the medium controls that 
    were run within a 20-run window (+/- 10 runs) of the experimental
    
    Keeps features that: 
    1. Are much higher in experimental than medium control
    2. Are much lower in experimental than medium control
    
    Inputs
    ------
    medium:         data_convert() output corresponding to medium 
                    control runs for one feature
    experimental:   data_convert() output corresponding to experimental
                    condition runs for one feature
    fold_change:    boundary to toss a feature if found to be 
                    fold_change higher in experimental over medium, or 
                    1/fold_change lower in experimental over medium; 
                    default 4
    
    Outputs
    -------
    boolean True or False, where True means that the feature is has a 
    very different intensity (as judged by fold_change) in experimental
    compared to the nearest control run
    """
    # for each medium control point, look for experimental points within
    # 10 runs of a control point
    for idx, ctrl in medium.iterrows():
        expt = experimental[(experimental['x'] < (ctrl['x']+10)) & \
                            (experimental['x'] > (ctrl['x']-10))]   
        # if there is an experimental point within 20 of a control
        if len(expt) > 0:
            for y in expt['y']:
                # if the control is zero, don't chuck it regardless of 
                # what experimental is
                if ctrl['y'] == 0.:
                    continue
                # otherwise, look at the fold difference between 
                # experimental and control
                else:
                    fold = y / ctrl['y']
                    # if the experimental is greater than 4x the control
                    # or less than .25x of the control, dont chuck it
                    if fold > fold_change or fold < (1 / fold_change):
                        continue
                    # if the experimental is ever within a 0.25 - 4x 
                    #range of the control, toss it immediately
                    else:
                        return False
    # if the experimental is always less than 1/4x of control or greater
    # than 4x of control across all control data points, then keep it
    return True

def cpd_from_inchikey(inchikey, blocks=1):
    '''
    Get a metabolite atlas compound from inchikey.
    
    Inputs
    ------
    inchikey:   a standard inchikey, e.g. RYYVLZVUVIJVGH-UHFFFAOYSA-N
    blocks:     which blocks of the inchikey to use as a search pattern 
                to find compound(s) in metatlas, example:
                1: RYYVLZVUVIJVGH
                2: RYYVLZVUVIJVGH-UHFFFAOYSA
                3: RYYVLZVUVIJVGH-UHFFFAOYSA-N

    '''
    if blocks not in [1,2,3]:
        raise RuntimeError('block argument must be 1, 2, or 3')

    if (inchikey != '') and (inchikey is not None):
        search_pattern = '-'.join(inchikey.split('-')[:blocks])
        cpds = metob.retrieve('Compounds', inchi_key='%s%%' \
                                %(inchikey.split('-')[0]), username='*')
        if len(cpds) != 0:
            return cpds[0].name
        else:
            return None
    else:
        return inchikey


