from itertools import product
from functools import partial
from nose2.tools import params
import unittest
from scsequtil import fastq
from scsequtil import reader
import string
import os
import copy
import numpy as np

# set some useful globals for testing
data_dir = os.path.split(__file__)[0] + '/data'
_i7_files = ('test_i7.fastq', 'test_i7.fastq.gz', 'test_i7.fastq.bz2')
_files = ('test_i7.fastq', 'test_r1.fastq', 'test_r2.fastq')
_gz_files = ('test_i7.fastq.gz', 'test_r1.fastq.gz', 'test_r2.fastq.gz')
_bz2_files = ('test_i7.fastq.bz2', 'test_r1.fastq.bz2', 'test_r2.fastq.bz2')

_modes = ('r', 'rb')
_files_and_modes = list(product(_i7_files, _modes))
_multifiles_and_modes = list(product((_files, _gz_files, _bz2_files), _modes))
_map_encoder = {'r': str, 'rb': partial(bytes, encoding='utf-8')}


class TestFastqReader(unittest.TestCase):

    @params(*_files_and_modes)
    def test_reader_opens_file_reads_first_line(self, filename, mode):
        rd = fastq.Reader('%s/%s' % (data_dir, filename), mode)
        record = next(iter(rd))
        expected_result = _map_encoder[mode]('NCACAATG\n')
        self.assertEqual(record.sequence, expected_result)

    def test_reader_correctly_stores_files(self):
        filenames = ['notreal', 'fake']
        rd = fastq.Reader(filenames, 'r')
        self.assertEqual(rd.filenames, filenames)

        rd = fastq.Reader(filenames[0], 'r')
        self.assertEqual(rd.filenames, [filenames[0]])

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

    @params(*_multifiles_and_modes)
    def test_multifile_read(self, filenames, mode):
        """
        loop over multiple files

        :param filenames:
        :return:
        """
        rd = fastq.Reader(
            ['%s/%s' % (data_dir, f) for f in filenames],
            mode=mode,
            header_comment_char='#')

        self.assertEqual(len(rd), 300)

    @params(*product((_files,), _modes))
    def test_mixed_filetype_read(self, filenames, mode):
        """
        loop over multiple files

        :param filenames:
        :return:
        """
        rd = fastq.Reader(
            ['%s/%s' % (data_dir, f) for f in filenames],
            mode=mode,
            header_comment_char='#')

        self.assertEqual(len(rd), 300)

    def test_wrong_type_raises_exception(self):
        """

        :param filenames:
        :param mode:
        :return:
        """
        self.assertRaises(TypeError, fastq.Reader, 10, 'r')
        self.assertRaises(
            TypeError,
            fastq.Reader,
            ('works', 10), 'r'  # arguments
        )
        self.assertRaises(
            ValueError,
            fastq.Reader,
            'works', 'not_acceptable_open_mode'
        )

    def test_fastq_returns_correct_filesize(self):
        rd = fastq.Reader(
            '%s/%s' % (data_dir, _i7_files[0]),
            mode='r',  # mode irrelevant
            header_comment_char='#')
        self.assertEqual(rd.size, 7774)

        rd = fastq.Reader(
            ['%s/%s' % (data_dir, f) for f in _i7_files],
            mode='r',  # mode irrelevant
            header_comment_char='#')
        self.assertEqual(rd.size, 7774 + 853 + 802)

    @params(*zip(_files, [8., 26., 98.]))
    def test_fastq_reader_estimates_correct_sequence_length(self, file_, expected_length):
        rd = fastq.Reader(
            '%s/%s' % (data_dir, file_),
            mode='r',  # mode irrelevant
            header_comment_char='#')
        mean, std, (lengths, counts) = rd.estimate_sequence_length()
        self.assertEqual(mean, expected_length)
        self.assertEqual(std, 0)
        self.assertTrue(np.array_equal(lengths, np.array([expected_length])))
        self.assertTrue(np.array_equal(counts, np.array([100])))

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
        """
        verify that reader correctly generates each field, and that the fields are properly structured.

        :param filename:
        :param mode:
        :return:
        """
        encoder = _map_encoder[mode]
        rd = fastq.Reader(
            '%s/%s' % (data_dir, filename),
            mode=mode)
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

    @params(*_modes)
    def test_mktag(self, mode):
        """
        test that mktag returns a properly formatted tag.

        :param mode:
        :return:
        """
        encoder = _map_encoder[mode]
        record_class = fastq.BytesRecord if mode == 'rb' else fastq.StrRecord
        string_tag = 'TC', encoder('ACGTACGT')  # could be cell barcode
        float_tag = 'TQ', 37.5  # could be an average quality score
        int_tag = 'TN', 10  # could be an additional bit flag
        self.assertEqual(
            record_class._mktag(*string_tag),
            encoder('TC:Z:ACGTACGT')
        )
        self.assertEqual(
            record_class._mktag(*float_tag),
            encoder('TQ:f:37.5')
        )
        self.assertEqual(
            record_class._mktag(*int_tag),
            encoder('TN:i:10')
        )

    @params(*_files_and_modes)
    def test_addition_and_recovery_of_tags(self, filename, mode):
        """
        verify that reader correctly generates each field, and that the fields are properly structured.

        :param filename:
        :param mode:
        :return:
        """
        encoder = _map_encoder[mode]
        rd = fastq.Reader(
            '%s/%s' % (data_dir, filename),
            mode=mode)
        record1 = next(iter(rd))
        record2 = copy.deepcopy(record1)

        string_tag = 'TC', encoder('ACGTACGT')  # could be cell barcode
        float_tag = 'TQ', 37.5  # could be an average quality score
        int_tag = 'TN', 10  # could be an additional bit flag

        all_tags = (string_tag, float_tag, int_tag)

        record1.set_tag(*string_tag)
        record1.set_tag(*float_tag)
        record1.set_tag(*int_tag)

        record2.set_tags(all_tags)

        self.assertEqual(record1.get_tags(), record2.get_tags())
        self.assertEqual(record1.get_tags(), dict(all_tags))



if __name__ == "__main__":
    unittest.main()
