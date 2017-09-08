import pysam


class SubsetAlignments:

    def __init__(self, alignment_file, open_mode=None):
        """Wrapper for pysam/htslib that allows non-standard filtering of alignments in a bamfile

        :param str alignment_file: sam or bam file.
        :param str open_mode: optional, mode to read file. Will be autodetected by file type if
          the file contains the correct suffix for its type.
        """

        if open_mode is None:
            if alignment_file.endswith('.bam'):
                open_mode = 'rb'
            elif alignment_file.endswith('.sam'):
                open_mode = 'r'
            else:
                raise ValueError('could not autodetect file type for alignment_file %s '
                                 '(detectible suffixes: .sam, .bam)' % alignment_file)
        self._file = alignment_file
        self._open_mode = open_mode

    def indices_by_chromosome(self, n_specific, chromosome, include_other=0):
        """Return the list of first n_specific indices of reads aligned to selected chromosome.

        If desired, will also return non-specific indices in a second list (can serve as negative
        control reads).

        :param int n_specific: number of aligned reads to return indices for
        :param str chromosome: only reads from this chromosome are considered valid
        :param int include_other: optional (default=0), the number of unaligned reads
          to include in the file.
        :return [int]: list of aligned indicies
        :return [int]: list of unaligned indices, only returned if include_other is not zero.
        """

        # check chromosome
        valid_chromosomes = [str(i) for i in range(1, 23)] + ['M', 'X', 'Y']
        if chromosome not in valid_chromosomes:
            raise ValueError('chromsome %s not valid. Must be one of %r' %
                             (chromosome, valid_chromosomes))

        with pysam.AlignmentFile(self._file, self._open_mode) as fin:
            specific, nonspecific = 0, 0  # counters
            chromosome = str(chromosome)
            chromosome_indices = []
            other_indices = []

            for i, record in enumerate(fin):

                if not record.is_unmapped:  # record is mapped
                    if chromosome in record.reference_name and specific < n_specific:
                        chromosome_indices.append(i)
                        specific += 1
                    elif nonspecific < include_other:
                        other_indices.append(i)
                        nonspecific += 1
                elif nonspecific < include_other:  # record is not mapped
                    other_indices.append(i)
                    nonspecific += 1

                # check termination condition (we have the requisite number of reads
                if specific == n_specific and nonspecific == include_other:
                    break

        if specific < n_specific or nonspecific < include_other:
            print('Warning: Only %d unaligned and %d reads aligned to chromosome %s were found in' 
                  '%s' % (nonspecific, specific, chromosome, self._file))

        if include_other != 0:
            return chromosome_indices, other_indices
        else:
            return chromosome_indices
