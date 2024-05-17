import numpy as np
import seaborn as sns

### Basic Abstraction Functions ###

def axplot(figure, axes_list, in_data_ob, pltfun, gs, figure_ind, title):
    axes_list.append(figure.add_subplot(gs[figure_ind])) # may want to automate the index creation process
    pltfun(in_data_ob, axes_list[-1])
    axes_list[-1].set_title(title)
    return len(axes_list)

def axtitle(axes_list, x_title, y_title):
    axes_list[-1].set_xlabel(x_title)
    axes_list[-1].set_ylabel(y_title)

def violin(data, axes):
    sns.violinplot(y=data, color="0.8", ax=axes)
    sns.stripplot(y=data, size=1, jitter=True, zorder=1, ax=axes)

def scatter(data, axes):
    axes.scatter(data[0], data[1], s=2)

def bar(data, axes):
    axes.bar(data[0], data[1])
    
def hist(data, axes):
    axes.hist(data[0], bins=data[1], range=data[2], log=True)
    
def hist_bin(X):
    return np.histogram(X, bins=10, range=(0.0,1.0))
