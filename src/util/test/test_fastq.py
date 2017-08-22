from itertools import product
from functools import partial
from nose2.tools import params
import unittest
from util import fastq
from util import reader
import string
import os


# set some useful globals for testing
data_dir = os.path.split(__file__)[0] + '/data'
_files = ('test_i7.fastq', 'test_i7.fastq.gz', 'test_i7.fastq.bz2')
_modes = ('r', 'rb')
_files_and_modes = list(product(_files, _modes))
_map_encoder = {'r': str, 'rb': partial(bytes, encoding='utf-8')}


class TestFastqReader(unittest.TestCase):

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

    @params(*_files_and_modes)
    def test_reader_reads_correct_number_of_records(self, filename, mode):

        rd = fastq.Reader(
            '%s/%s' % (data_dir, filename),
            mode=mode,
            header_comment_char='#')

        self.assertEqual(len(rd), 100)

        rd = reader.Reader(
            '%s/%s' % (data_dir, filename),
            mode=mode,
            header_comment_char='#')

        self.assertEqual(len(rd), 400)

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

    @params(*_files_and_modes)
    def test_fields_populate_properly(self, filename, mode):
        encoder = _map_encoder[mode]
        rd = fastq.Reader(
            '%s/%s' % (data_dir, filename),
            mode=mode)
        _ = next(iter(rd))
        # test some things about each record
        name_prefix = encoder('@')
        alphabet = set(encoder('ACGTN'))
        name2_string = encoder('+\n')
        ascii_chars = set(i for i in encoder(string.printable))
        for record in rd:
            self.assertTrue(record.name.startswith(name_prefix))
            self.assertTrue(all(i in alphabet for i in record.sequence.strip()))
            self.assertTrue(record.name2 == name2_string)
            self.assertTrue(all(i in ascii_chars for i in record.quality.strip()))


if __name__ == "__main__":
    unittest.main()
