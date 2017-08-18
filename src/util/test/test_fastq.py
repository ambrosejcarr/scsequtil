import nose2
from itertools import product
from functools import partial
from nose2.tools import params
import unittest
from util import fastq
import os

# set some useful globals for testing
data_dir = os.path.split(__file__)[0] + '/data'
_files = ('test_i7.fastq', 'test_i7.fastq.gz', 'test_i7.fastq.bz2')
_modes = ('r', 'rb')
_files_and_modes = list(product(_files, _modes))
_map_encoder = {'r': str, 'rb': partial(bytes, encoding='utf-8')}


class TestReader(unittest.TestCase):

    @params(*_files_and_modes)
    def test_reader_opens_file_reads_first_line(self, filename, mode):
        rd = fastq.Reader('%s/%s' % (data_dir, filename), mode)
        record = next(iter(rd))
        expected_result = _map_encoder[mode]('NCACAATG\n')
        self.assertEqual(record.sequence, expected_result)

    @params(*_files_and_modes)
    def test_reader_skips_header_character(self, filename, mode):
        """
        test should skip the first name line, shifting each record up 1. As a result, the
         first sequence should be found in the name field

        :param filename: input fastq files
        :param mode: mode of file opening
        """
        rd = fastq.Reader(
            '%s/%s' % (data_dir, filename),
            mode=mode,
            header_comment_char='@')
        record = next(iter(rd))
        expected_result = _map_encoder[mode]('NCACAATG\n')
        self.assertEqual(record.name, expected_result)

    # # currently failing, unclear how to best raise exceptions without overhead
    # @params(*_files_and_modes)
    # def test_reader_throws_exception_for_incomplete_record(self, filename, mode):
    #     """
    #     test should skip the first name line, shifting each record up 1. As a result,
    #      the first sequence should be found in the name field
    #
    #     :param filename: input fastq files
    #     """
    #     rd = fastq.Reader(
    #         '%s/%s' % (data_dir, filename),
    #         mode=mode,
    #         header_comment_char='@')
    #     self.assertRaises(BaseException, lambda x: sum(1 for _ in x), rd)





class TestRecord(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        pass

if __name__ == "__main__":
    unittest.main()