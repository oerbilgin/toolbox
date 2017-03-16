import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy.stats import linregress
import math

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
    determines number of rows  for subplots based on the total number of plots 
    and number of columns.

    For use with matplotlib.subplots():
        n_rows, n_cols = subplot_row_col(7)
        fig, axs = plt.subplots(n_rows, n_cols)
    
    Inputs
    ------
    n_plots: total number of subplots
    n_cols: how many columns you want
    
    Outputs
    ------
    n_rows: integer value of number of rows to encompass the total number of 
            plots
    """
    a = np.floor(n_plots / float(n_cols))
    b = n_plots % float(n_cols)
    if b != 0.:
        n_rows = int(a + 1)
    else: n_rows = int(a)
    return n_rows
def trendline(x, y, poly=1, CI=False, PI=False, sigma=2):
    """
    Returns trendline coordinates and optionally Confidence Interval (CI) 
    and/or Prediction Interval (PI) and points within and out of CI or PI

    From Wikipedia:
    prediction intervals predict the distribution of individual future points, 
    whereas confidence intervals and credible intervals of parameters predict 
    the distribution of estimates of the true population mean or other quantity
    of interest that cannot be observed
    
    Adapted from:
    http://stackoverflow.com/questions/28505008/numpy-polyfit-how-to-get-1-sigma-uncertainty-around-the-estimated-curve
    
    Inputs
    ------
    x: x coordinate data of sample
    y: y coordinate data of sample
    poly: degree of polynomial fit
    err: False to just get x_fit and y_fit coordinates; True to get error bound
    sigma: corresponds to confidence interval:
        2= 95% CI
        3= 99% CI
    
    Outputs
    -------
    x: sorted x values from sample data
    trend_fn: numpy.poly1d function of the trendline
    sigma_y_fit: error of y_fit; returned when err==True
    within_bounds: 2D array of [x] and [y] arrays corresponding to data that 
                    fall within the error bounds
    outof_bounds: 2D array of [x] and [y] arrays corresponding to data that 
                    fall outside of the error bounds
    
    Example Usage (linear regression)
    ---------------------------------
    x = np.array([4.0,2.5,3.2,5.8,7.4,4.4,8.3,8.5])
    y = np.array([2.1,4.0,1.5,6.3,5.0,5.8,8.1,7.1])
    xfit, trend_fn, err, within, without = trendline(x,y, err=True)
    fig, ax = plt.subplots()
    # data points within error
    ax.plot(within[0], within[1], color='lightskyblue', marker='o', 
            linestyle='') 
    # data points within error
    ax.plot(without[0], without[1], color='coral', marker='o', linestyle='') 
    # linear regression
    ax.plot(xfit, trend_fn(xfit), color='k') 
    # + 2sigma
    ax.plot(xfit, trend_fn(xfit)+err, color='k', linestyle='--') 
    # - 2sigma
    ax.plot(xfit, trend_fn(xfit)-err, color='k', linestyle='--') 
    """
    # sort the data by x dimension
    sorted_idx = np.argsort(x)
    x = x[sorted_idx]
    y = y[sorted_idx]
    
    # calculate the trendline
    p, cov = np.polyfit(x, y, poly, cov=True)
    trend_fn = np.poly1d(p)
    
    if not err:
        return x, trend_fn
    else:
        ## this is copied almost directly from the SO post
        TT = np.vstack([x**(poly-i) for i in range(poly+1)]).T
        # matrix multiplication calculates the polynomial values
        y_fit = np.dot(TT, p)
        # C_y = TT*C_z*TT.T
        cov_y_fit = np.dot(TT, np.dot(cov, TT.T)) 
        # Standard deviations are sqrt of diagonal
        sigma_y_fit = np.sqrt(np.diag(cov_y_fit)) * sigma 
        ##
        
        # determine which data points are in or out of bounds
        within_bounds = []
        outof_bounds = []
        for i, xval in enumerate(x):
            if (y[i] > trend_fn(xval) - sigma_y_fit[i]) and \
                y[i] < (trend_fn(xval) + sigma_y_fit[i]):
                within_bounds.append([xval, y[i]])
            else:
                outof_bounds.append([xval, y[i]])
                
        # transform the arrays
        within_bounds = np.asarray(within_bounds).T
        outof_bounds = np.asarray(outof_bounds).T
        return x, trend_fn, sigma_y_fit, within_bounds, outof_bounds

def calc_perp(A, B, C):
    """
    Given a line segment connecting points A-B and a third point C, 
        calculate where point D is on the line that makes a perpendicular 
        line to connect C to A-B
    
    Inputs
    ------
    A, B, and C are (x,y) coordinates for points A, B, and C

    Outputs
    -------
    Dx, Dy: x and y coordinates of point D that forms a perpendicular line from
            AB to C
    
    from:
    http://stackoverflow.com/questions/10301001/perpendicular-on-a-line-segment-from-a-given-point
    http://stackoverflow.com/questions/1811549/perpendicular-on-a-line-from-a-given-point
    """  
    Ax,Ay = A
    Bx,By = B
    Cx,Cy = C
    t= ((Cx-Ax)*(Bx-Ax) + (Cy-Ay)*(By-Ay)) / ((Bx-Ax)**2 + (By-Ay)**2)
    Dx = Ax + t*(Bx-Ax)
    Dy = Ay + t*(By-Ay)
    return [Dx, Dy]

def twod_dist(a, b):
    """
    Calculates the distance between two points in 2-dimensional space

    Inputs
    ------
    a, b: [x,y] coordinates for the two points; can be list or tuple

    Outputs
    -------
    dist: distance
    """
    try:
        if len(a) != 2 or len(b) != 2:
            raise RuntimeError('Too many values in one or both inputs')
    except TypeError:
        raise RuntimeError('The inputs must be [x,y] coordinates')
    
    a = [float(x) for x in a]
    b = [float(x) for x in b]
    
    dist = np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)
    
    return dist


def elbow_point(x_vals, y_vals):
    """
    calculates the elbow point of a plot by finding the point furthest from the
        line created between first and last values

    Inputs
    ------
    x_vals, y_vals: x and y values of the plot that forms the elbow

    Outputs
    -------
    elbow_idx: The index of x_vals and y_vals that corresponds to the data
                point that is the elbow

    from:
    http://stackoverflow.com/questions/4033821/using-a-smoother-with-the-l-method-to-determine-the-number-of-k-means-clusters
    """
    # get the formula for the hypotenuse
    A = [x_vals[0], x_vals[-1]]
    B = [y_vals[0], y_vals[-1]]
    z = np.polyfit(A, B, 1)
    p = np.poly1d(z)
    maxdist=0
    for i, y in enumerate(y_vals):
        C = [x_vals[i], y]
        D = calc_perp(A, B, C)
        dist = twod_dist(C, D)
        if dist > maxdist:
            maxdist = dist
            elbow_idx = i
    return elbow_idx
