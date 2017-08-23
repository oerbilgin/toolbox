import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy.stats import t
from matplotlib.patches import Ellipse
import matplotlib.patches as mpatches
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def hist_silhouette(hist_obj, ignore_zero=True):
    """
    Takes in a matplotlib plt.hist() or numpy.histogram() output 
    and returns x and y lists corresponding to a line connecting the 
    tops of the bins
    
    plt.hist() should return a tuple where:
    0: array of bin heights
    1: array of bin edges
    2: matplotlib patches

    numpy.histogram() returns a tuple where
    0: array of bin heights
    1: array of bin edges

    Inputs
    ------
    hist_obj:       tuple where the first element is a list of bin heights,
                    and the second element is a list of bin edges
    ignore_zeros:   When true, the output x and y lists do not have
                    values for histogram bars with a height of 0

    Outputs
    -------
    x_list: list of x-coordinates to be plotted, corresponding to the 
            middle point between two bin edges; the center of the bin bar
    y)list: list of y-coordinates to be plotted; the height of the bin 
            bar
    """
    counts = hist_obj[0]
    bin_edges = hist_obj[1]
    y_list = []
    x_list = []
    for i, pt in enumerate(counts):
        y = pt
        x = (bin_edges[i+1] + bin_edges[i]) / 2
        if ignore_zero:
            if y != 0:
                y_list.append(y)
                x_list.append(x)
        else:
            y_list.append(y)
            x_list.append(x)
    return x_list, y_list

def cmap_color(curve_idx, n_curves, color_map='viridis'):
    """
    returns the suitable color from the given color map for the plot 
    number. Meant for use when plotting multiple curves or scatters in
    an iterative loop. 

    Example Usage:
    --------------
    data = [list of dicts corresponding to x and y coordinates]
    for i, subdata in enumerate(data):
        plt.plot(subdata['x'], subdata['y'], 
                 color=cmap_color(i, len(data)))
        OR:
        plt.scatter(subdata['x'], subdata['y'], 
                    color=cmap_color(i, len(data)))

    If you instead want to have color map plotting within one data
    series, plt.scatter() has its own color map functionality for this:
    x = [0,1,2,3,4]
    y = [0,1,2,3,4]
    plt.scatter(x, y, c=range(len(x)), cmap='viridis')
    
    Inputs
    ------
    curve_idx:  index position of the curve you are getting the plot for
    n_curves:   total number of curves you want to plot with the color 
                map
    color_map:  a valid matplotlib color map. Available choices:
    
    cmaps = 
    [('Perceptually Uniform Sequential',
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

    Outputs
    -------
    Hex color string.

    """
    cmap = plt.get_cmap(color_map)
    c_norm = colors.Normalize(vmin=0, vmax=range(n_curves)[-1])
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=cmap)
    color = colors.rgb2hex(scalar_map.to_rgba(curve_idx))
    return color

def subplot_row_col(n_plots, n_cols=5):
    """
    determines number of rows  for subplots based on the total number of
    plots and number of columns.

    For use with matplotlib.subplots():
        n_rows, n_cols = subplot_row_col(7)
        fig, axs = plt.subplots(n_rows, n_cols)
    
    Inputs
    ------
    n_plots:    total number of subplots
    n_cols:     how many columns you want
    
    Outputs
    ------
    n_rows: integer value of number of rows to encompass the total 
            number of plots
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
    calculate where point D is on the line that makes a 
    perpendicular line to connect C to A-B
    
    Inputs
    ------
    A, B, and C are (x,y) coordinates for points A, B, and C

    Outputs
    -------
    Dx, Dy: x and y coordinates of point D that forms a perpendicular 
            line from AB to C
    
    from:
    http://stackoverflow.com/questions/10301001/perpendicular-on-a-line-
        segment-from-a-given-point
    http://stackoverflow.com/questions/1811549/perpendicular-on-a-line-
        from-a-given-point
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
    calculates the elbow point of a plot by finding the point furthest 
    from the line created between first and last values

    Inputs
    ------
    x_vals, y_vals: x and y values of the plot that forms the elbow

    Outputs
    -------
    elbow_idx:  The index of x_vals and y_vals that corresponds to the 
                data point that is the elbow

    from:
    http://stackoverflow.com/questions/4033821/using-a-smoother-with-
        the-l-method-to-determine-the-number-of-k-means-clusters
    """
    # get the formula for the hypotenuse
    line_x = [x_vals[0], x_vals[-1]]
    line_y = [y_vals[0], y_vals[-1]]
    
    A = [line_x[0], line_y[0]]
    B = [line_x[1], line_y[1]]
    
    z = np.polyfit(line_x, line_y, 1)
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

def trendline(x, y, poly=1, confint=False, conf=0.95):
    """
    Returns trendline coordinates and optionally confidence interval 

    The standard deviation lines currently only work for linear 
    regressions.
    
    Confidence intervals show the confidence of the **fit**
    
    Adapted from:
    https://github.com/KirstieJane/STATISTICS/blob/master/CIs_LinearRegression.py
    
    Inputs
    ------
    x:          x coordinate data of sample
    y:          y coordinate data of sample
    poly:       degree of polynomial fit
    confint:    True to get bounds for confidence interval **of the 
                trendline**
    conf:       what confidence interval you want 
    
    Outputs
    -------
    x:              sorted x values from sample data
    trend_fn:       numpy.poly1d function of the trendline
    
    Example Usage (linear regression with all bells and whistles)
    -------------
    xfit, trend_fn, conf_dy = trendline(x, y, confint=True, conf=2)
    
    # plot data
    plt.scatter(x, y, c='lightskyblue')

    # plot trendline
    plt.plot(xfit, trend_fn(xfit), color='k') 

    # Plot confidence intervals of fit
    plt.plot(xfit, trend_fn(xfit)+conf_dy, color='k', linestyle='--') 
    plt.plot(xfit, trend_fn(xfit)-conf_dy, color='k', linestyle='--') 
    """
    # sort the data by x dimension
    sorted_idx = np.argsort(x)
    x = x[sorted_idx]
    y = y[sorted_idx]
    
    # calculate the trendline
    p = np.polyfit(x, y, poly)
    trend_fn = np.poly1d(p)
    
    if not confint:
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

        conf_dy = tstat * np.sqrt((s_err/(n-2))*(1.0/n + \
                    (np.power((x-mean_x),2)/
                    ((np.sum(np.power(x,2)))-n*(np.power(mean_x,2))))))
        return x, trend_fn, conf_dy

def multivar_ellipse(x,y, sigma=2, 
                    ellipse_color='none', 
                    ellipse_border_color='black', 
                    ellipse_linestyle='--', 
                    ellipse_hatch=''):
    """
    Given x and y data, returns an ellipse patch and the points in and 
    out of the ellipse
    
    Adapted from http://stackoverflow.com/questions/20126061/creating-a-
    confidence-ellipses-in-a-sccatterplot-using-matplotlib
    
    Inputs
    ------
    x,y: x and y data
    sigma: number of standard deviations of data to encompass
    ellipse_color: face color of ellipse patch
    ellipse_border: border color of ellipse patch
    
    Outputs
    -------
    
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    cov = np.cov(x, y)
    vals, vecs = eigsorted(cov)
    
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    
    w, h = 2 * sigma * np.sqrt(vals)
    
    ell = Ellipse(xy=(np.mean(x), np.mean(y)),
                  width=w, height=h,
                  angle=theta, color=ellipse_border_color,
                  linestyle=ellipse_linestyle, hatch=ellipse_hatch)
    ell.set_facecolor(ellipse_color)

    within_bounds = []
    outof_bounds = []
    for i, xval in enumerate(x):
        if ell.contains_point([xval, y[i]]):
            within_bounds.append([xval, y[i]])
        else:
            outof_bounds.append([xval, y[i]])
    
    # Transform into 2xN array, where each array is the x or y coords
    # For easier plotting
    within_bounds = np.asarray(within_bounds).T
    outof_bounds = np.asarray(outof_bounds).T
    
    return ell, within_bounds, outof_bounds

def stddev_from_trendline(x, y, sigma=4):
    """
    Calculates the distance from a linear regression that 
    encompasses a certain number of standard deviations of the sample 
    set being modeling
    
    First calculates the perpendicular distance of each data point to 
    the linear regression. Then calculates the standard deviation of 
    these distances, and multiplies by the given value sigma. Then
    calculates the vertical offset from the trendline corresponding to a
    line parallel to the trendline that is the standard deviation times 
    sigma perpendicular distance away.
    
    Note: Only works for linear regression!
    
    Key portion adapted from http://stackoverflow.com/questions/133897/
        how-do-you-find-a-point-at-a-given-perpendicular-distance-from-
        a-line
    
    Inputs
    ------
    x, y:   x and y arrays for the data being modeled
    sigma:  number of standard deviations
    
    Outputs
    -------
    std_dy:     the vertical offset from the regression line 
                corresponding to the number of standard deviations away 
                from the linear regression
    trend_fn:   np.poly1d() object representing the linear regression of
                the data
    
    Example Usage
    -------------
    dy, trend_fn = stddev_from_trendline(x, y)
    # data
    plt.scatter(x, y) 
    # trendline
    plt.plot(x, trend_fn(x)) 
    # std line upper bound
    plt.plot(x, trend_fn(x)+dy, linestyle='--') 
    # std line lower bound 
    plt.plot(x, trend_fn(x)-dy, linestyle='--') 

    # find lines that are outside this range
    idx_dwn = np.where(y < trend_fn(x)-dy)[0]
    idx_up = np.where(y > trend_fn(x)+dy)[0]
    idx = np.concatenate((idx_dwn, idx_up))
    outside_x = x[idx]
    outside_y = y[idx]

    #plot them
    plt.scatter(outside_x, outside_y, s=100, color='', edgecolors='k')
    """
    p = np.polyfit(x, y, 1)
    trend_fn = np.poly1d(p)
    
    minidx = np.where(x == x.min())[0]
    maxidx = np.where(x == x.max())[0]
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
    
    # return dy for given sigma for plotting a line
    perp_dxdy = sigma * dstd # perpendicular distance from trendline

    # calculate dy for std
    # Calculate unit vector perpendicular to trendline
    tdist = twod_dist(a, b)
    dx = (a[0] - b[0]) / tdist
    dy = (a[1] - b[1]) / tdist
    # a point that is the appropriate perpendicular distance away from 
    # the beginning of the trendline
    q = [(a[0] + perp_dxdy*dy), (a[1] - perp_dxdy*dx)]

    # find the vertical distance from point q to the trendline
    std_dy = abs(q[1] - trend_fn(q[0]))
    
    return std_dy, trend_fn

def trendline_deprecated(x, y, poly=1, confint=False, conf=0.95, sigline=False,
                         sigma=2):
    """
    Returns trendline coordinates and optionally confidence interval 
    and/or lines that encompass standard deviation(s) of the data, 
    and points within and out of these standard deviation bounds. 

    The standard deviation lines currently only work for linear 
    regressions.
    
    Confidence intervals show the confidence of the **fit**, whereas the
    standard deviation lines show the range around the trendline 
    encompassing a certain number of standard deviations around the 
    average data point distance from the trendline
    
    Adapted from:
    http://stackoverflow.com/questions/28505008/numpy-polyfit-how-to-get
        -1-sigma-uncertainty-around-the-estimated-curve
    
    Inputs
    ------
    x:          x coordinate data of sample
    y:          y coordinate data of sample
    poly:       degree of polynomial fit
    confint:    True to get bounds for confidence interval **of the 
                trendline**
    conf:       what confidence interval you want 
    sigline:    True to get bounds for lines that encompass standard 
                deviation(s) **of the data**
    sigma:      how many standard deviations of error to return
    
    Outputs
    -------
    x:              sorted x values from sample data
    trend_fn:       numpy.poly1d function of the trendline
    sigma_y_fit:    error of y_fit; returned when err==True
    within_bounds:  2D array of [x] and [y] arrays corresponding to data
                    that fall within the error bounds
    outof_bounds:   2D array of [x] and [y] arrays corresponding to data
                    that fall outside of the error bounds
    
    Example Usage (linear regression with all bells and whistles)
    -------------
    xfit, trend_fn, conf_dy = trendline(x, y, confint=True, conf=2)
    
    # plot data
    plt.scatter(x, y, c='lightskyblue')

    # plot trendline
    plt.plot(xfit, trend_fn(xfit), color='k') 

    # Plot confidence intervals of fit
    plt.plot(xfit, trend_fn(xfit)+conf_dy, color='k', linestyle='--') 
    plt.plot(xfit, trend_fn(xfit)-conf_dy, color='k', linestyle='--') 
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

        conf_dy = tstat * np.sqrt((s_err/(n-2))*(1.0/n + \
                (np.power((x-mean_x),2)/
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
