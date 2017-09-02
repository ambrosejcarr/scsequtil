import os
import gzip
import bz2
from collections.abc import Iterable


class Reader:
    """
    Basic reader object that seamlessly loops over multiple input files

    Can be subclassed to create readers for specific file types (fastq, gtf, etc.)
    """

    def __init__(self, files_='-', mode='r', header_comment_char=None):
        """

        :param list|str files_: file or list of files to be read. Defaults to sys.stdin
        :param mode: open mode. Default 'r' will return string objects. Change to 'rb' to
          return bytes objects.
        """

        if isinstance(files_, str):
            self._files = [files_]
        elif isinstance(files_, Iterable):  # test items of iterable
            files_ = list(files_)
            if all(isinstance(f, str) for f in files_):
                self._files = files_
            else:
                raise TypeError('all passed files must be type str')
        else:
            raise TypeError('files_ must be a string filename or a list of such names.')

        # set open mode:
        if mode not in {'r', 'rb'}:
            raise ValueError('mode must be one of r, rb')
        self._mode = mode

        if isinstance(header_comment_char, str) and mode == 'rb':
            self._header_comment_char = header_comment_char.encode()
        else:
            self._header_comment_char = header_comment_char

    @property
    def filenames(self):
        return self._files

    def __len__(self):
        """
        return the length of the Reader object. Note that for sys.stdin, this will
        consume the input
        """
        return sum(1 for _ in self)

    def __iter__(self):
        for file_ in self._files:

            # set correct mode
            if self._mode == 'r' and any(file_.endswith(s) for s in ['.gz', '.bz2']):
                mode = 'rt'
            else:
                mode = self._mode

            # open file
            if file_.endswith('.gz'):
                f = gzip.open(file_, mode)
            elif file_.endswith('.bz2'):
                f = bz2.open(file_, mode)
            else:
                f = open(file_, mode)

            # iterate over the file, dropping header lines if requested
            try:
                file_iterator = iter(f)
                if self._header_comment_char is not None:
                    first_record = next(file_iterator)
                    while first_record.startswith(self._header_comment_char):
                        first_record = next(file_iterator)

                    yield first_record  # avoid loss of first non-comment line

                for record in file_iterator:  # now, run to exhaustion
                    yield record
            finally:  # clean up
                f.close()

    @property
    def size(self):
        """return the collective size of all files being read in bytes"""
        return sum(os.stat(f).st_size for f in self._files)

    def estimate_length(self):
        """should estimate length of all provided files."""
        return NotImplemented
