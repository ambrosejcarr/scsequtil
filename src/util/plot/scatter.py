import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import gaussian_kde
from collections import Iterable


def density_2d(x, y):
    """return x and y and their density z, sorted by their density (smallest to largest)

    :param np.ndarray x: coordinate data
    :param np.ndarray y: coordinate data
    :return (np.ndarray,): sorted x, y, and density
    """
    xy = np.vstack([np.ravel(x), np.ravel(y)])
    z = gaussian_kde(xy)(xy)
    return np.ravel(x), np.ravel(y), np.arcsinh(z)


def map_categorical_to_cmap(data, cmap=plt.get_cmap()):
    """
    create a discrete colormap from cmap appropriate for data

    :param np.ndarray data: categorical vector to map to colormap
    :param str|matplotlib.colors.ListedColormap cmap: cmap to discretize, or 'random'
    :return np.ndarray, dict: vector of colors matching input data, dictionary of labels
      to their respective colors
    """
    categories = np.unique(data)
    n = len(categories)
    if isinstance(cmap, str) and 'random' in cmap:
        colors = np.random.rand(n, 3)
    else:
        colors = cmap(np.linspace(0, 1, n))
    category_to_color = dict(zip(categories, colors))
    return np.array([category_to_color[i] for i in data]), category_to_color


def add_legend_to_categorical_vector(
        colors, labels, ax, loc='best',
        markerscale=0.75, **kwargs):
    """
    Add a legend to a plot where the color scale was set by discretizing a colormap.

    :param colors: np.ndarray, output of map_categorical_vector_to_cmap()
    :param labels: np.ndarray, category labels
    :param ax: axis on which the legend should be plotted
    :param str loc: location of legend
    :param float markerscale: scale for legend markers
    :param kwargs: additional kwargs for legend
    :return: None
    """
    artists = []
    for c in colors:
        artists.append(plt.Line2D((0, 1), (0, 0), color=c, marker='o', linestyle=''))
    ax.legend(
        artists, labels, loc=loc, markerscale=markerscale, **kwargs)


def categorical(
        x, y, c, ax=None, cmap=plt.get_cmap(), legend=True, legend_kwargs=None,
        randomize=True, remove_ticks=False, *args, **kwargs):
    """
    Wrapper for plt.scatter wherein the output is colored by the labels in c.

    Data Inputs:
    :param np.ndarray x: coordinate data to be scattered (required)
    :param np.ndarray y: coordinate data to be scattered (required)
    :param Iterable c: iterable of len == x, y containing data labels (required)

    Optional Plotting Parameters:
    :param mpl.axes._subplots.AxesSubplot ax: axis on which to scatter data
    :param str|matplotlib.colors.ListedColormap cmap: color map or 'random' to plot
      random categories
    :param bool legend: if True, plot legend
    :param bool randomize: if True, randomize order of plotting
    :param bool remove_ticks: if True, removes axes ticks and labels
    :param dict legend_kwargs: additional kwargs for legend
    :param list args: additional args for scatter
    :param dict kwargs: additional kwargs for scatter

    :return mpl.axes._subplots.AxesSubplot: axis containing scattered data
    """
    if not ax:
        ax = plt.gca()

    if legend_kwargs is None:
        legend_kwargs = dict()

    color_vector, category_to_color = map_categorical_to_cmap(np.asarray(c), cmap)

    if randomize:
        ind = np.random.permutation(len(x))
    else:
        ind = np.argsort(np.ravel(c))

    ax.scatter(np.ravel(x)[ind], np.ravel(y)[ind], c=color_vector[ind], *args,
               **kwargs)
    if remove_ticks:
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())

    labels, colors = zip(*sorted(category_to_color.items()))
    if legend:
        add_legend_to_categorical_vector(colors, labels, ax, markerscale=2,
                                         **legend_kwargs)
    return ax


def continuous(x, y, c=None, ax=None, colorbar=True, randomize=True,
               remove_ticks=False, *args, **kwargs):
    """
    Wrapper for plt.scatter wherein data is plotted on a continuous color spectrum
    based on values of c.

    Input Data Parameters:
    :param np.ndarray x: x coordinate data (required)
    :param np.ndarray y: y coordinate data (required)
    :param Iterable c:  continuous vector used to color data points

    Opitional Plotting Parameters
    :param mpl.axes._subplots.AxesSubplot ax: axis on which to scatter data
    :param bool colorbar: if True, include a colorbar
    :param bool randomize: if False, sort by C, plotting highest values last
    :param bool remove_ticks: if True, remove ticks from plot
    :param list args: additional arguments to pass to plt.scatter
    :param dict kwargs: additional keyword arguments to pass to plt.scatter

    :return mpl.axes._subplots.AxesSubplot: axis containing scattered data
    """

    if ax is None:
        ax = plt.gca()

    if c is None:  # plot density if no color vector is provided
        x, y, c = density_2d(x, y)

    if randomize:
        ind = np.random.permutation(len(x))
    else:
        ind = np.argsort(c)

    sm = ax.scatter(x[ind], y[ind], c=c[ind], *args, **kwargs)
    if remove_ticks:
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_major_locator(plt.NullLocator())
    if colorbar:
        cb = plt.colorbar(sm)
        cb.ax.xaxis.set_major_locator(plt.NullLocator())
        cb.ax.yaxis.set_major_locator(plt.NullLocator())
    return ax
