# Python Cheatsheet
## Pandas
### NumPy `islcose()`
* Incredibly useful for things like matching mz, rt info between mzmine, metatlas, etc.
```python
mz = 123.4567
# find all detected m/z within .001 daltons
df[np.isclose(df['mz'], mz, atol=.001)]
# could alternately calculate atol with my ppm functions for better matching
```

### String methods reminders
* Useful options are `case=False` and `regex=False`. Turning regex off is sometimes critical to avoid errors, but then you can't search for multiple patterns using `contains()`.

* **.str.contains()** 
	* uses '**|**' operator in between patterns as an `OR` statement, only works if `regex=True`
	* I do not know of a way to do an `AND` statement within the function
* **.str.startswith()** and **.str.endswith()**
	* Both of them use **tuples** to separate patterns with an OR statement

## Decorators
From https://pabloariasal.github.io/python-decorators-from-the-ground-up/

Example decorator for timing a function

```python
import time

# this is the "decorating" function
def timer(func):
    def inner_func(*args, **kwargs):
        start = time.time() # set up timing
        func_out = func(*args, **kwargs) # do the funciton
        print time.time() - start # print the log
        return func_out # return the function's output
    return inner_func # return the function's output

# decorate a function
@timer
def test_func(a, b, c='jedi'):
    return a + b + c

# use the function as normal
>>> test_func('the', 'last')
1.19209289551e-06
'thelastjedi'
```

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
This opens a number of python processes, which on some systems can cause memory issues because it *appears* that the amount of memory used has increased by the number of processes open.
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
plt.rcParams['pdf.fonttype'] = 42 # makes editing text in Illustrator much easier
%matplotlib inline # only when in an Ipython notebook
# it is always better to use axes instead of plt for more control of 
# multiple axes and/or multiple plots
ax1, fig = plt.subplots() 
ax1.plot(x, y, color=, label=, marker=, linestyle=, linewidth=)
ax1.set_ylabel('label')
ax1.set_xlabel('label')
plt.title('title')
plt.show()
```

### color bar
From https://stackoverflow.com/questions/8342549/matplotlib-add-colorbar-to-a-sequence-of-line-plots
```python
fig, ax = plt.subplots()
... # whatever plot you want to make
sm = plt.cm.ScalarMappable(
	cmap='viridis', # any color map
	norm=plt.Normalize(vmin=0, vmax=1) # vmin and vmax are your map's min and max
	)
sm._A = []
cbar = fig.colorbar(sm, orientation='vertical')
cbar.ax.set_ylabel('tanimoto score')
cba.ax.set_yticks() # can do normal axes functions this way too
```

### Fun axes functions
```python
ax1.set_yscale('log') # try to capture any y-values that are zero, otherwise they won't be plotted
ax1.set_ylim(low, high)
ax2 = ax1.twinx() # for separate y-axes
ax2.plot(x, y)
```
### Tick labels
```python
# using plot object
plt.xticks(tick_positions, labels)

# using axes object
ax.set_xticks(tick_positions)
ax.set_xticklabels(labels)

# rotation (kwargs are the same for plt or ax)
# set 'ha' to 'right' to align right when rotating between 0 and 90
plt.xticks(tick_positions, labels, rotation=45, ha='right')
```

### Basic legend stuff
* This will place the legend box outside the plot, with the upper left corner of the plot just to the right of the upper right corner of the plot
```python
ax1.legend(loc='upper left', bbox_to_anchor=(1.1, 1.))
```
The `loc` argument says how to orient the `bbox_to_anchor`, and bbox_to_anchor is in units of x and y axes positions: (0,0) is bottom left, (1,1) it top right

* This will make a legend most useful for bar charts
```python
import matplotlib.patches as mpatches
handle_list = []
sample_colors = ['red', 'green', 'blue']
for i, label in enumerate(['control', 'sample_1', 'sample_2']):
    patch = mpatches.Patch(color=sample_colors[i], label=label)
    handle_list.append(patch)
ax.legend(handles=handle_list, loc='lower left', bbox_to_anchor=(-1., 0))
```

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
