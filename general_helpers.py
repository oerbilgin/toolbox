from collections import defaultdict
import math
import pandas as pd

def load_dataframe(fname, filetype=None, key=None):
    """
    Uses the appropriate pandas function to load a file based on the
    file extension:
    .pkl or .pickle: pickle file
    .xls or .xlsx: excel file
    .csv: comma separated file
    .txt, .tab, or .tsv: tab separated file
    .h5 or .hdf5: HDF5 formatted files. To load these, you must also
                  pass a key argument!

    Inputs
    ------
    fname: path to the file to be loaded as a pandas dataframe
    filetype: Used to override the autodetection of a filetype based on
              its file extension. Accepts file extensions listed above,
              but without the preceding "." (e.g. pass it "pkl")
    key: They key for the table to be loaded in the HDF5 file. Only
         required/used when loading HDF5 tables. 

    Outputs
    -------
    df: the loaded pandas dataframe

    WARNING: for binary file types (pickle and hdf5), be wary of
    saving and loading from different versions of pandas, as this will
    very likely break the loader.
    """

    if filetype is None:
        file_ext = fname.split('.')[-1]
    else:
        file_ext = filetype
    if file_ext in ['pkl', 'pickle']:
        df = pd.read_pickle(fname)
    elif file_ext in ['xls', 'xlsx']:
        df = pd.read_excel(fname)
    elif file_ext in ['csv']:
        df = pd.read_csv(fname)
    elif file_ext in ['txt', 'tab', 'tsv']:
        df = pd.read_csv(fname, sep='\t')
    elif file_ext in ['h5', 'hdf5']:
        if key is None:
            raise IOError('"key" argument must be used when loading\
                an HDF5 file')
        else:
            df = pd.read_hdf(fname, key)
    else:
        raise IOError('could not infer what type of file %s is... please \
            make it a csv, tab, or pickle file')

    # remove rows that are empty
    df = df[pd.notnull(df).all(axis=1)]
    return df

def partition_indexes(totalsize, numberofpartitions, offset=0):
	"""
	Used to split up an iterable into multiple chunks
	This function returns the indices to create chunks of roughly the
	same size

	From:
	http://stackoverflow.com/questions/33878939/
		get-indices-of-roughly-equal-sized-chunks

	Inputs
	------
	totalsize: length of the iterable to partition
	numberofpartitions: number of partitions to split the iterable into
	offset: number to add to the partition indices. Used when you're
			partitioning a slice of something, for example if you want
			to parititon array[100:200] into 4 partitions, this function
			would return [[0, 25], [25, 50], [50, 75], [75, 100]] if
			offset were kept default. So you would set offset to 100,
			which would return
			[[100, 125], [125, 150], [150, 175], [175, 200]]
	"""

	chunk_indices = []
	# Compute the chunk size (integer division; i.e. assuming Python 2.7)
	chunksize = totalsize / numberofpartitions
	# How many chunks need an extra 1 added to the size?
	remainder = totalsize - chunksize * numberofpartitions
	a = 0
	for i in xrange(numberofpartitions):
		b = a + chunksize + (i < remainder)
		# Yield the inclusive-inclusive range
#         chunk_indices.append([a, b - 1])
		chunk_indices.append([a, b])
		a = b

	offset_chunk_indices = []
	for indices in chunk_indices:
		i = indices[0] + offset
		j = indices[1] + offset
		offset_chunk_indices.append([i, j])

	return offset_chunk_indices
def upper_tri_index_unravel(i, j, n, printer=False):
    '''
    Formula to get linear index of an upper triangle of a square matrix (nxn) while ignoring the main diagonal
    
    Input:
    i, j: coordinates in the upper triangle (first one is (0,1); (0,0) doesn't exist)
    n: shape of the square (e.g. for 4x4, n=4)
    
    Output: single value of the linear index
    
    Formula from:
    http://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
    '''
    
    k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
    
    if printer:
        print 'Coordinates:', (i,j)
        print 'Linear index:', k
    
    return k
def upper_tri_index_ravel(k, n):
    """
    Formulae to get coordinates from a linear index of a square matrix while ignoring the main diagonal

    k is the linear index, n is the length of one side of matrix

    Formula from:
    http://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
    """
    i = int(n - 2 - math.floor(math.sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5))
    j = int(k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2)
    return (i, j)

def etree_to_dict(t):
	"""
	from https://stackoverflow.com/questions/2148119/how-to-convert-an-xml-string-to-a-dictionary-in-python

	input needs to be made like this:
		from xml.etree import cElementTree as ET
		xml_text = <XML STRING>
		xml_code = ET.XML(xml_text)
		xml_dict = etree_to_dict(xml_code)

	output is a dictionary
	"""
	d = {t.tag: {} if t.attrib else None}
	children = list(t)
	if children:
		dd = defaultdict(list)
		for dc in map(etree_to_dict, children):
			for k, v in dc.items():
				dd[k].append(v)
		d = {t.tag: {k:v[0] if len(v) == 1 else v for k, v in dd.items()}}
	if t.attrib:
		d[t.tag].update(('@' + k, v) for k, v in t.attrib.items())
	if t.text:
		text = t.text.strip()
		if children or t.attrib:
			if text:
				d[t.tag]['#text'] = text
		else:
			d[t.tag] = text
	return d
