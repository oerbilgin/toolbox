# Python Cheatsheet
## Argument Parsing
```python
import argparse
# parse arguments
parser = argparse.ArgumentParser()
# a couple examples
parser.add_argument('-l', '--level', help='how many levels deep to search the chemical network', type=int, choices=[1,2,3], default=2, required=True)
parser.add_argument('--mute', help='mutes pandas warnings', action='store_true')
# store arguments
args = parser.parse_args()
# use arguments
args.level
args.mute
```
## Multiprocessing
This opens a numver of python processes, which on some systems can cause memory issues because it *appears* that the amount of memory used has increased by the number of processes open.
I'm sure there's a way to work around this.
```python
import multiprocessing as mp
max_cpu = mp.cpu_count()
p = mp.Pool(max_cpu)
out = p.map(func, iterable)
p.close()
p.terminate()
```

## Matplotlib

### First, some of my favorite matplotlib colors to use when you only need to use a few colors:
```
coral, yellowgreen, lightskyblue, gold, gray, black
```
### Basic plotting
```python
import matplotlib.pyplot as plt
%matplotlib inline #only when in an Ipython notebook
# it is always better to use axes instead of plt for more control of 
# multiple axes and/or multiple plots
ax1, fig = plt.subplots() 
ax1.plot(x, y, color=, label=, marker=, linestyle=, linewidth=)
ax1.set_ylabel('label')
ax1.set_xlabel('label')
plt.title('title')
plt.show()
```

### Fun axes functions
```python
ax1.set_yscale('log') # try to capture any y-values that are zero, otherwise they won't be plotted
ax1.set_ylim(low, high)
ax2 = ax1.twinx() # for separate y-axes
ax2.plot(x, y)
```

### Basic legend stuff
This will place the legend box outside the plot, with the upper left corner of the plot just to the right of the upper right corner of the plot
```python
ax1.legend(loc='upper left', bbox_to_anchor=(1.1, 1.))
```
The `loc` argument says how to orient the `bbox_to_anchor`, and bbox_to_anchor is in units of x and y axes positions: (0,0) is bottom left, (1,1) it top right

### Setting figure size
```python
fig.set_size_inches(12,6) # (width, height)
```

### Subplots
It is really useful to ravel subplots, so you can act on them individually
```python
fig, axs = plt.subplots(rows, cols, sharex=True, sharey=True)
axs = axs.ravel()
for i, coordinates in list_of_coordinates:
	x = coordinates[0]
	y = coordinates[1]
	axs[i].plot(x, y)
	axs[i].set_title('individual subplot title')

### Adjust the spacing of subplots
fig.subplots_adjust(wspace=, **kwargs) # I find this to be way better than tight_layout(), and usually use wspace, but there are other kwargs that give you more control
fig.suptitle('title of overall plot', fontsize=24)
```

### Fun error options
Scatterplot with error bars
```python
ax.errorbar(x, y, yerr=, xerr=, fmt=data_marker, color=fmt_color, ecolor=error_bar_color, *kwargs) # many cool kwargs to customize error bars
```

Showing error via shading (good for time series)
```python
ax.plot(x, y) # main data
ax.fill_between(x, y-error, y+error, color=, alpha=0.2)
```
