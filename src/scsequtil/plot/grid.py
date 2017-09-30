import numpy as np
import matplotlib.pyplot as plt


class AxesGrid:

    def __init__(self, n, nrows=None, ncols=None, *args, **kwargs):
        """
        Generate a grid of figures to plot on.

        :param int n: total number of plots. If smaller than nrow * ncol, additional plots
          are turned off.
        :param int nrows: number of rows of axes
        :param int ncols: number of columns of axes
        :param list args: additional positional arguments to pass to plt.subplots
        :param dict kwargs: additional keyword arguments to pass to plt.subplots
        """
        self._n = n
        if not all([nrows, ncols]):
            nrows = ncols = np.ceil(np.sqrt(n)).astype(int)
        self._nrows = nrows
        self._ncols = ncols
        self._figure = None
        self._axes = None

        self._figure, self._axes = plt.subplots(
            nrows=self._nrows, ncols=self._ncols, *args, **kwargs)

    @property
    def figure(self):
        return self._figure

    @property
    def axes(self):
        return self._axes

    def __iter__(self):
        for ax in np.ravel(self._axes):
            yield ax

    def __getitem__(self, item):
        return self._axes[item]

    def plot_all(self, arguments, plot_function):
        """

        :param Iterable arguments: tuples of arguments to pass to plot_function as *args
        :param plot_function: plotting function for each axis, must take axis as a
          keyword argument (ax=ax)
        """
        for ax, args in zip(self, arguments):
            plot_function(*args, ax=ax)
