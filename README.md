# toolbox
Helpful python functions aggregated into one place

## Mass Spectrometry helpers
* MetAtlas is an import, but is not required for most functions
* rdkit is an import, but is currently not required for any functions
* Most of these functions are designed to work with MetAtlas, Pactolus, or MAGI outputs

### `ppm_error()`
Returns ppm error between two masses

### `ppm_window()`
Returns mass-ppm and mass+ppm

### `find_mass()`
Given a mass and ppm error, and optionally a retention time and window, gives you dataframe row(s) that have features that match your search

### `data_convert()`
Converts a Pandas.Series object into dataframe for easier plotting

### `ctrl_expt_compare()`
Removes peaks that were found in two samples

### `cpd_from_inchikey()`
Uses MetAtlas database to get a compound object from an InChI-key

## Matplotlib helpers
### `cmap_color()`
Color maps one liner
[cmap example Notebook](https://github.com/oerbilgin/toolbox/blob/master/example_notebooks/cmap.ipynb)

![cmap example](https://github.com/oerbilgin/toolbox/blob/master/images/cmap.png)

### `subplot_row_col()`
Makes plotting lots of subplots easy

### `calc_perp()`
Calculates a point on a line perpendicular to an arbitrary point

### `twod_dist()`
Calculates distance between two 2D points

### `elbow_point()`
Finds the elbow point of an elbow plot

[elbowpoint example Notebook](https://github.com/oerbilgin/toolbox/blob/master/example_notebooks/elbowfinder.ipynb)

![elbowpoint example](https://github.com/oerbilgin/toolbox/blob/master/images/elbow.png)

### `trendline()`
Calculates a polynomial fit, optionally gives you confidence intervals for the polynomial

[trendline example Notebook](https://github.com/oerbilgin/toolbox/blob/master/example_notebooks/trendlines.ipynb)

![trendline example](https://github.com/oerbilgin/toolbox/blob/master/images/trend.png)

### `stddev_from_trendline()`
Calculates the linear regression of a 2D dataset, and gives you the vertical distance from the trendline that encompasses n-standard deviations of the data

[trendline example Notebook](https://github.com/oerbilgin/toolbox/blob/master/example_notebooks/trendlines.ipynb)

![trendline example](https://github.com/oerbilgin/toolbox/blob/master/images/linreg_std.png)

### `multivar_ellipse()`
Draws an ellipse around multivariate data that encompasses n-sigma of the data

[ellipse example Notebook](https://github.com/oerbilgin/toolbox/blob/master/example_notebooks/ellipse.ipynb)

![ellipse example](https://github.com/oerbilgin/toolbox/blob/master/images/ellipse.png)

