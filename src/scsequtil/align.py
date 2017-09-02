

class STAR:

    def __init__(self, input_=None, output_prefix=None):
        """
        multipurpose STAR wrapper, allows:
        - alignment of input file and exit
        - return command to chain star using util_old.mp.Chain
        - return open subprocess to chain star

        :param input_:
        :param output_prefix:
        """
        self._input = input_
        self._output_prefix = output_prefix

    def align(self):
        """align reads and exit"""
        raise NotImplementedError

    def command(self):
        """return a command formatted for use with util_old.mp.Chain"""
        raise NotImplementedError

    def subprocess(self):
        """return a chainable subprocess instance"""
        raise NotImplementedError

    @staticmethod
    def remove_index(index):
        """remove a sharedmemory index"""
        raise NotImplementedError
