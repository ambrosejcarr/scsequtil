from nose2.tools import params
import os
import unittest
from util import gtf

_data_dir = os.path.split(__file__)[0] + '/data'
_files = ['%s/%s' % (_data_dir, f) for f in ('test.gtf', 'test.gtf.gz', 'test.gtf.bz2')]


class TestGTFReader(unittest.TestCase):

    @params(*_files)
    def test_opens_file_reads_first_line(self, filename):
        rd = gtf.Reader(filename, 'r', header_comment_char='#')
        record = next(iter(rd))

        print(record)


if __name__ == '__main__':
    unittest.main()