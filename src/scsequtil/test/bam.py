import unittest
from nose2.tools import params
import os
from scsequtil.bam import SubsetAlignments
import pysam

# test files have 4446 chr 19 and 873 chr 21 alignments
# test files are sorted, so chr 19 alignments come first
data_dir = os.path.split(__file__)[0] + '/data'
_files = (data_dir + '/test.bam', data_dir + '/test.sam')


class TestSubsetAlignments(unittest.TestCase):

    def test_incorrect_name_raises(self):
        filename = 'test.notabam'
        self.assertRaises(ValueError, SubsetAlignments, filename)

    @params(*_files)
    def test_indices_are_properly_ordered(self, filename):

        # get length of file:
        mode = 'r' if filename.endswith('.sam') else 'rb'
        length = sum(1 for _ in pysam.AlignmentFile(filename, mode))

        n_specific = 20
        n_non_specific = 20
        sa = SubsetAlignments(filename)
        ind_specific, ind_nonspecific = sa.indices_by_chromosome(
            n_specific, '19', include_other=n_non_specific)

        # todo put in better error messages

        # should get indices that are all less than the file's length
        self.assertLessEqual(max(ind_specific + ind_nonspecific), length)

        # chromosome 19 indices should all come before chromosome 21
        self.assertLess(max(ind_specific), min(ind_nonspecific), 'spec: %r \n non: %r \n' % (
            ind_specific, ind_nonspecific))

        # should get correct number of indices provided n < 4446
        self.assertEqual(n_specific, len(ind_specific))
        # should get correct number of indices provided n < 873
        self.assertEqual(n_non_specific, len(ind_nonspecific))

        # first chromosome 19 index should come at start of file
        self.assertEqual(min(ind_specific), 0)

        ind_specific, ind_nonspecific = sa.indices_by_chromosome(
            n_specific, '21', include_other=n_non_specific)

        # should get indices that are all less than the file's length
        self.assertLessEqual(max(ind_specific + ind_nonspecific), length)

        # chromosome 19 indices should all come before chromosome 21
        self.assertLess(max(ind_nonspecific), min(ind_specific), 'spec: %r \n non: %r \n' % (
            ind_specific, ind_nonspecific))

        # should get correct number of indices provided n < 873
        self.assertEqual(n_specific, len(ind_specific))

        # should get correct number of indices provided n < 4446
        self.assertEqual(n_non_specific, len(ind_nonspecific))

        # first chromosome 21 index should come after all chromosome 19 indices
        self.assertEqual(min(ind_specific), 4446)
