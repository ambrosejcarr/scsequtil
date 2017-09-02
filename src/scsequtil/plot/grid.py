import matplotlib.pyplot as plt
import matplotlib as mpl


class AxesGrid:

    def __init__(self, n, nrows, ncols, figure_size):
        """
        Generate a grid of figures to plot on.

        :param int n: total number of plots. If smaller than nrow * ncol, additional plots
          are turned off.
        :param int nrows: number of rows of axes
        :param int ncols: number of columns of axes
        """
        self._n = n
        self._nrows = nrows
        self._ncols = ncols
        self._size = figure_size
        self._figure = None
        self._axes = None

        self._figure, self._axes = plt.subplots(
            nrows=self._nrows, ncols=self._ncols, figsize=self._size)

    def __iter__(self):
        for ax in self._axes:
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
            plot_function(*arguments, ax=ax)
