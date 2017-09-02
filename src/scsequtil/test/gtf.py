from nose2.tools import params
import os
import unittest
from scsequtil import gtf
from itertools import chain

_data_dir = os.path.split(__file__)[0] + '/data'
_files = ['%s/%s' % (_data_dir, f) for f in ('test.gtf', 'test.gtf.gz', 'test.gtf.bz2')]


class TestGTFReader(unittest.TestCase):

    @params(*_files)
    def test_opens_file_reads_first_line(self, filename):
        rd = gtf.Reader(filename, 'r', header_comment_char='#')
        _ = next(iter(rd))

    @params(*_files)
    def test_opens_file_populates_fields_properly(self, filename):
        rd = gtf.Reader(filename, 'r', header_comment_char='#')
        record = next(iter(rd))
        self.assertEqual(record.seqname, 'chr19')
        self.assertEqual(record.chromosome, 'chr19')
        self.assertEqual(record.source, 'HAVANA')
        self.assertEqual(record.feature, 'gene')
        self.assertEqual(record.start, 60951)
        self.assertEqual(record.end, 71626)
        self.assertEqual(record.score, '.')
        self.assertEqual(record.strand, '-')
        self.assertEqual(record.frame, '.')

        expected_features = {
            'gene_id': 'ENSG00000282458.1',
            'gene_type': 'transcribed_processed_pseudogene',
            'gene_status': 'KNOWN',
            'gene_name': 'WASH5P',
            'level': '2',
            'havana_gene': 'OTTHUMG00000180466.8',
        }
        self.assertEqual(record._attribute, expected_features)

        self.assertTrue(all(i in str(record) for i in chain(expected_features.keys(), expected_features.values())))

    @params(*_files)
    def test_set_attribute_verify_included_in_output_string(self, filename):
        rd = gtf.Reader(filename, 'r', header_comment_char='#')
        record = next(iter(rd))
        record.set_attribute('test_attr', 'foo')
        self.assertEqual(record.get_attribute('test_attr'), 'foo')

        # verify in output string
        self.assertTrue('foo' in str(record))

    @params(*_files)
    def test_opens_file_parses_size(self, filename):
        rd = gtf.Reader(filename, 'r', header_comment_char='#')
        record = next(iter(rd))
        self.assertEqual(71626 - 60951, record.size)

        # mangle record, make sure error is raised
        record._fields[3:5] = [record.end, record.start]
        self.assertRaises(ValueError, getattr, record, 'size')


if __name__ == '__main__':
    unittest.main()