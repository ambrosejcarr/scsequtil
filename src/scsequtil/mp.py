from multiprocessing import Pool as Pool_
from multiprocessing import cpu_count
from contextlib import closing
from collections import Iterable
from functools import partial
import subprocess
import shlex
import tempfile


class Pool:

    def __init__(self, func, iterable, ncpu=None, **kwargs):
        """

        :param func: a function that takes a single positional argument, which is filled
          by iterable, and any number of keyword arguments, passed as keyword arguments
          to the constructor of Pool
        :param Iterable iterable: the iterable to map over func.
        :param int ncpu:
        :param dict kwargs: keyword arguments to set as defaults for func.
        """
        self._function = func
        self._iterable = iterable
        self._kwargs = kwargs
        if ncpu is None:
            ncpu = cpu_count()
        self._ncpu=ncpu
        self._result = None

    @property
    def result(self):
        if self._result is not None:
            return self._result
        else:
            print('Run map(), imap(), or imap_unordered() first to calculate result '
                  'object.')

    def map(self):
        """Wrapper for Multiprocessing.Pool.map()

        Documentation for map():
        {doc}

        :return:
        """.format(doc=Pool_.map.__doc__)

        mapfunc = partial(self._function, **self._kwargs)
        with closing(Pool_(processes=self._ncpu)) as pool:
            self._result = pool.map(mapfunc, self._iterable)
        return self.result

    def imap_unordered(self):
        """Wrapper for Multiprocessing.Pool.imap_unordered()
        
        Documentation for imap_unordered():
        {doc}

        :return:
        """.format(doc=Pool_.imap_unordered.__doc__)
        mapfunc = partial(self._function, **self._kwargs)
        with closing(Pool_(processes=self._ncpu)) as pool:
            self._result = pool.imap_unordered(mapfunc, self._iterable)
        return self.result

    def imap(self):
        """Wrapper for Multiprocessing.Pool.imap()

        Documentation for imap():
        {doc}

        :return:
        """.format(doc=Pool_.imap.__doc__)
        mapfunc = partial(self._function, **self._kwargs)
        with closing(Pool_(processes=self._ncpu)) as pool:
            self._result = pool.imap(mapfunc, self._iterable)
        return self.result


class Chain:

    def __init__(self, functions, stdin=None):
        """Convenience class to chain together a series of functions through pipes, using
        Popen's communicate() method to avoid deadlocks.

        :param list functions: a list of functions to be chained
        :param str|bytes stdin: Optional. if provided, data will be sent to first function
          in list.
        """

        # make sure functions are properly formatted
        cmds = [f if isinstance(f, list) else shlex.split(f) for f in functions]
        self._cmds = cmds

        if isinstance(stdin, str):
            stdin = stdin.encode()
        elif isinstance(stdin, bytes):
            pass
        else:
            raise TypeError('stdin must be bytes or str')
        self._stdin = stdin
        self._result = None
        self._error = None

    @property
    def result(self):
        return self._result.decode()

    @property
    def error(self):
        return self._error.decode()

    def run(self):
        """

        :return:
        """

        # create a temporary file for each process to store stderr
        error_files = [tempfile.TemporaryFile() for _ in self._cmds]

        # weave functions together
        processes = [
            subprocess.Popen(
                self._cmds[0],
                stdin=self._stdin,
                stdout=subprocess.PIPE,
                stderr=error_files[0])
        ]

        try:
            for i, func in enumerate(self._cmds[1:]):
                last_p = processes[i + 1]  # get previous function
                p = subprocess.Popen(
                    func,
                    stdin=last_p.stdout,
                    stdout=subprocess.PIPE,
                    stderr=error_files[i + 1])
                processes.append(p)

            # get output
            self._result, _ = p.communicate()

            # get error dictionary
            [f.seek(0) for f in error_files]
            self._error = {
                cmd: f.read().decode() for (cmd, f) in zip(self._cmds, error_files)
            }

        finally:  # close all the files and processes
            for f in error_files:
                f.close()
            for p in processes:
                p.terminate()
