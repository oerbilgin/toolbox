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
