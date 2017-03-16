import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy.stats import t

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

def calculate_stddev_lines(x, y, sigma=2):
    """
    Calculates the vertical distance from a linear regression that encompasses a certain number of standard deviations of the sample set being modeling
    
    Note: Only works for linear regression!
    
    Key portion adapted from http://stackoverflow.com/questions/133897/how-do-you-find-a-point-at-a-given-perpendicular-distance-from-a-line
    
    Inputs
    ------
    x, y: x and y arrays for the data being modeled
    sigma: number of standard deviations
    
    Outputs
    -------
    std_dy: an array of same dimensions as x and y that is the vertical offset from the regression line
    trend_fn: np.poly1d() object representing the linear regression of the data
    
    Example Usage
    -------------
    std_dy = calculate_stddev_lines(x, y)
    plt.scatter(x, y) # data
    plt.plot(x, trend_fn(y)) # trendline
    plt.plot(x, trend_fn(y)+std_dy, linestyle='--') # std line upper bound
    plt.plot(x, trend_fn(y)-std_dy, linestyle='--') # std line lower bound 
    plt.show()
    """
    p = np.polyfit(x, y, 1)
    trend_fn = np.poly1d(p)
    
    minidx = np.where(x == x.min())
    maxidx = np.where(x == x.max())
    minx = x[minidx][0]
    maxx = x[maxidx][0]
    a = [minx, trend_fn(minx)] # trendline low-left
    b = [maxx, trend_fn(maxx)] # trendline upper-right

    distances=[]
    for i, xval in enumerate(x):
        c = [xval, y[i]]
        d = calc_perp(a, b, c)
        distances.append(twod_dist(c, d))
    distances=np.asarray(distances)
    davg = distances.mean() # mean distance from points to trendline
    dstd = distances.std() # standard deviation
    perp_dxdy = sigma * dstd # perpendicular distance from trendline

    ## Adapted from http://stackoverflow.com/questions/133897/how-do-you-find-a-point-at-a-given-perpendicular-distance-from-a-line
    # calculate dy for std
    # Calculate unit vector perpendicular to trendline
    tdist = twod_dist(a, b)
    dx = (a[0] - b[0]) / tdist
    dy = (a[1] - b[1]) / tdist
    # a point that is the appropriate perpendicular distance away from the beginning of the trendline
    q = [(a[0] + perp_dxdy*dy), (a[1] - perp_dxdy*dx)]

    # find the vertical distance from point q to the trendline
    std_dy = abs(q[1] - trend_fn(q[0]))
    
    return std_dy, trend_fn

def trendline(x, y, poly=1, confint=False, conf=0.95, sigline=False, sigma=2):
    """
    Returns trendline coordinates and optionally confidence interval and/or lines that encompass standard deviation(s) of the data, and points within and 
    out of these standard deviation bounds. The standard deviation lines currently only work for linear regressions.
    
    Confidence intervals show the confidence of the **fit**, whereas the standard deviation lines show the range around the trendline encompassing a certain number of standard deviations around the average data point distance from the trendline
    
    Adapted from:
    http://stackoverflow.com/questions/28505008/numpy-polyfit-how-to-get-1-sigma-uncertainty-around-the-estimated-curve
    
    Inputs
    ------
    x: x coordinate data of sample
    y: y coordinate data of sample
    poly: degree of polynomial fit
    confint: True to get bounds for confidence interval **of the trendline**
    conf: what confidence interval you want 
    sigline: True to get bounds for lines that encompass standard deviation(s) **of the data**
    sigma: how many standard deviations of error to return
    
    Outputs
    -------
    x: sorted x values from sample data
    trend_fn: numpy.poly1d function of the trendline
    sigma_y_fit: error of y_fit; returned when err==True
    within_bounds: 2D array of [x] and [y] arrays corresponding to data that fall 
                    within the error bounds
    outof_bounds: 2D array of [x] and [y] arrays corresponding to data that fall 
                    outside of the error bounds
    
    Example Usage (linear regression with all bells and whistles)
    -------------
    xfit, trend_fn, conf_dy, sig_dy, within, without = trendline(x,y, ci=True, pi=True, sigma=2)
    # plot trendline
    plt.plot(xfit, trend_fn(xfit), color='k') 

    # Plot confidence intervals of fit
    plt.plot(xfit, trend_fn(xfit)+conf_dy, color='k', linestyle='--') 
    plt.plot(xfit, trend_fn(xfit)-conf_dy, color='k', linestyle='--') 

    # Plot lines that encompass ~ 95% of data
    plt.plot(xfit, trend_fn(xfit)-sig_dy, color='r', linestyle='--')
    plt.plot(xfit, trend_fn(xfit)+sig_dy, color='r', linestyle='--') 
    
    # plot data points that fall within or out of the range
    plt.plot(within[0], within[1], color='lightskyblue', marker='o', linestyle='') # data points within the range
    plt.plot(without[0], without[1], color='coral', marker='o', linestyle='') # data points out of the range

    ToDo:
    [ ] Add ability to calculate stddev lines for higher order polynomial fits
    """
    # sort the data by x dimension
    sorted_idx = np.argsort(x)
    x = x[sorted_idx]
    y = y[sorted_idx]
    
    # calculate the trendline
    p = np.polyfit(x, y, poly)
    trend_fn = np.poly1d(p)
    
    if not confint and not sigline:
        return x, trend_fn
    else:
        ## copied almost directly from 
        ## https://github.com/KirstieJane/STATISTICS/blob/master/CIs_LinearRegression.py
        # Calculate error of fit
        y_err = y - trend_fn(x)
        
        # Calculate the confidence interval
        mean_x = np.mean(x) # mean of x
        n = len(x) # number of samples in original fit
        tstat = t.ppf(conf, n-1) # find appropriate t value
        s_err = np.sum(np.power(y_err,2)) # sum of the squares of the residuals

        conf_dy = tstat * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((x-mean_x),2)/
                ((np.sum(np.power(x,2)))-n*(np.power(mean_x,2))))))
    if not sigline:
        return x, trend_fn, conf_dy
    elif poly == 1:
        # Calculate dy from curve that encompasses n*sigma of data points
        # Calculate distance of each datapoint from the fit curve
        sig_dy, _ = calculate_stddev_lines(x, y, sigma=sigma)
            
        # Determine which data points are in or out of bounds
        within_bounds = []
        outof_bounds = []
        for i, xval in enumerate(x):
            if (y[i] > trend_fn(xval) - sig_dy) and \
                y[i] < (trend_fn(xval) + sig_dy):
                within_bounds.append([xval, y[i]])
            else:
                outof_bounds.append([xval, y[i]])
                
        # transform the arrays
        within_bounds = np.asarray(within_bounds).T
        outof_bounds = np.asarray(outof_bounds).T
        
        if confint and sigline:
            return x, trend_fn, conf_dy, sig_dy, within_bounds, outof_bounds
        else:
            return x, trend_fn, sig_dy, within_bounds, outof_bounds
    else:
        raise RuntimeError('Currently cannot calculate population stdev lines for anything other than linear regression, sorry')