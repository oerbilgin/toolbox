import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy.stats import linregress

def cmap_color(curve_idx, n_curves, color_map='viridis'):
    """
    returns the suitable color from the given color map for the plot number
    
    Inputs
    ------
    curve_idx: index position of the curve you are getting the plot for
    n_curves: total number of curves you want to plot with the color map
    color_map: a valid matplotlib color map. Available choices:
    
    cmaps = [('Perceptually Uniform Sequential',
                            ['viridis', 'inferno', 'plasma', 'magma']),
         ('Sequential',     ['Blues', 'BuGn', 'BuPu',
                             'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
                             'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
                             'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']),
         ('Sequential (2)', ['afmhot', 'autumn', 'bone', 'cool',
                             'copper', 'gist_heat', 'gray', 'hot',
                             'pink', 'spring', 'summer', 'winter']),
         ('Diverging',      ['BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
                             'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
                             'seismic']),
         ('Qualitative',    ['Accent', 'Dark2', 'Paired', 'Pastel1',
                             'Pastel2', 'Set1', 'Set2', 'Set3', 'Vega10',
                             'Vega20', 'Vega20b', 'Vega20c']),
         ('Miscellaneous',  ['gist_earth', 'terrain', 'ocean', 'gist_stern',
                             'brg', 'CMRmap', 'cubehelix',
                             'gnuplot', 'gnuplot2', 'gist_ncar',
                             'nipy_spectral', 'jet', 'rainbow',
                             'gist_rainbow', 'hsv', 'flag', 'prism'])]

    """
    cmap = plt.get_cmap('viridis')
    c_norm = colors.Normalize(vmin=0, vmax=range(n_curves)[-1])
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=cmap)
    color = scalar_map.to_rgba(curve_idx)
    return color

def subplot_row_col(n_plots, n_cols=5):
    """
    determines number of rows  for subplots based on the total number of plots and number of columns.

    For use with matplotlib.subplots():
        n_rows, n_cols = subplot_row_col(7)
        fig, axs = plt.subplots(n_rows, n_cols)
    
    Inputs
    ------
    n_plots: total number of subplots
    n_cols: how many columns you want
    
    Outputs
    ------
    n_rows: integer value of number of rows to encompass the total number of plots
    """
    a = np.floor(n_plots / float(n_cols))
    b = n_plots % float(n_cols)
    if b != 0.:
        n_rows = int(a + 1)
    else: n_rows = int(a)
    return n_rows
def trendline(x, y, poly=1, err=False, sigma=2):
    """
    Returns trendline coordinates and optionally error
    
    Adapted from:
    http://stackoverflow.com/questions/28505008/numpy-polyfit-how-to-get-1-sigma-uncertainty-around-the-estimated-curve
    
    Inputs
    ------
    x: x coordinate data of sample
    y: y coordinate data of sample
    poly: degree of polynomial fit
    err: False to just get x_fit and y_fit coordinates; True to get error bound
    sigma: how many standard deviations of error to return
    
    Outputs
    -------
    x: sorted x values from sample data
    y_fit: y values for the best fit
    sigma_y_fit: error of y_fit; returned when err==True
    
    Example Usage (linear regression)
    ---------------------------------
    xfit, yfit, err = trendline(x,y, err=True)
    fig, ax = plt.subplots()
    ax.plot(x, y, color='coral', marker='o', linestyle='') # data points
    ax.plot(xfit, yfit, color='k') # linear regression
    ax.plot(xfit, yfit+err, color='k', linestyle='--') # + 2sigma
    ax.plot(xfit, yfit-err, color='k', linestyle='--') # - 2sigma
    """
    p, cov = np.polyfit(x, y, poly, cov=True)
    sorted_idx = np.argsort(x)
    x = x[sorted_idx]
    
    ## this is copied almost directly from the SO post
    TT = np.vstack([x**(poly-i) for i in range(n+1)]).T
    y_fit = np.dot(TT, p)  # matrix multiplication calculates the polynomial values
    cov_y_fit = np.dot(TT, np.dot(cov, TT.T)) # C_y = TT*C_z*TT.T
    sigma_y_fit = np.sqrt(np.diag(cov_y_fit)) * sigma # Standard deviations are sqrt of diagonal
    ##
    
    if not err:
        return x, y_fit
    else:
        return x, y_fit, sigma_y_fit
