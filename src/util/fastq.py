import numpy as np
from . import reader


class Record:
    """Fastq record object

    Defines several properties for accessing fastq record information:
    :property name: name field
    :property sequence: sequence field
    :property name2: second name field
    :property quality: quality field

    Also defines several methods for accessing SEQC annotation fields:
    :property annotations: list of annotations
    :property metadata: dictionary of read metadata (if any present)
    :property average_quality: return the mean quality of FastqRecord
    """

    __slots__ = ['_data', '_tags']

    def __init__(self, record):
        """

        :param [str|bytes] record: list of four strings or bytes objects containing a
          single fastq record
        """
        self._data = list(record)

    @property
    def name(self):
        return self._data[0]

    @name.setter
    def name(self, value):
        self._data[0] = value

    @property
    def sequence(self):
        return self._data[1]

    @sequence.setter
    def sequence(self, value: bytes):
        self._data[1] = value

    @property
    def name2(self):
        return self._data[2]

    @name2.setter
    def name2(self, value):
        self._data[2] = value

    @property
    def quality(self):
        return self._data[3]

    @quality.setter
    def quality(self, value):
        self._data[3] = value

    def __bytes__(self):
        try:
            return b''.join(self._data)
        except TypeError:
            return ''.join(self._data).encode()

    def __str__(self):
        try:
            return ''.join(self._data)
        except TypeError:
            return b''.join(self._data).decode()

    def __repr__(self):
        return "Name: %s\nSequence: %s\nName2: %s\nQuality: %s\n" % (
            self._data[0], self._data[1], self._data[2], self._data[3])

    def __len__(self):
        return len(self.sequence)

    def average_quality(self):
        raise NotImplementedError

    def get_tags(self):
        raise NotImplementedError

    def get_tag(self, key):
        raise NotImplementedError

    def set_tag(self, key, value):
        raise NotImplementedError

    def set_tags(self, keys, values):
        raise NotImplementedError


class BytesRecord(Record):
    """Fastq record object

    Defines several properties for accessing fastq record information:
    :property name: name field
    :property sequence: sequence field
    :property name2: second name field
    :property quality: quality field

    Also defines several methods for accessing SEQC annotation fields:
    :property annotations: list of annotations
    :property metadata: dictionary of read metadata (if any present)
    :property average_quality: return the mean quality of FastqRecord
    """

    _type_map = {int: b'i',
                 str: b'Z',
                 float: b'f'}
    _itype_map = {v: k for k, v in _type_map.items()}

    def get_tags(self):
        """
        :return list: annotations present in the fastq header
        """
        return {
            k: self._itype_map[t](v) for k, t, v in
            (v.split(b':') for v in self.name[1:].split(b';')[-1])}

    def get_tag(self, key):
        try:
            return self.get_tags()[key]
        except KeyError:
            return None

    @classmethod
    def _mktag(cls, key, value):
        """
        construct a tag given a key, value pair
        :param key:
        :param value:
        :return str:
        """
        return '%s:%s:%s' % (
            key,
            cls._type_map[type(value)],
            str(value))

    def set_tag(self, key, value):
        """prepends a tag to annotation name
        :param bytes key:
        :param Object value:
        """
        self.name = '@%s;%s' % (self._mktag(key, value), self.name[1:])

    def set_tags(self, keys, values):
        """

        :param keys:
        :param values:
        :return:
        """
        '@%s;%s' % (
            ';'.join(self._mktag(k, v) for k, v in zip(keys, values)),
            self.name[1:]
        )

    def average_quality(self):
        """"""
        return (
            np.mean(np.frombuffer(self.quality, dtype=np.int8, count=len(self)))
            .astype(int) - 33
        )


class StrRecord(Record):
    """Fastq record object

    Defines several properties for accessing fastq record information:
    :property name: name field
    :property sequence: sequence field
    :property name2: second name field
    :property quality: quality field

    Also defines several methods for accessing SEQC annotation fields:
    :property annotations: list of annotations
    :property metadata: dictionary of read metadata (if any present)
    :property average_quality: return the mean quality of FastqRecord
    """

    _type_map = {int: 'i',
                 str: 'Z',
                 float: 'f'}
    _itype_map = {v: k for k, v in _type_map.items()}

    def get_tags(self):
        """
        :return list: annotations present in the fastq header
        """
        return {
            k: self._itype_map[t](v) for k, t, v in
            (v.split(':') for v in self.name[1:].split(';')[-1])}

    def get_tag(self, key):
        try:
            return self.get_tags()[key]
        except KeyError:
            return None

    @classmethod
    def _mktag(cls, key, value):
        """
        construct a tag given a key, value pair
        :param key:
        :param value:
        :return str:
        """
        return '%s:%s:%s' % (
            key,
            cls._type_map[type(value)],
            str(value))

    def set_tag(self, key, value):
        """prepends a tag to annotation name
        :param bytes key:
        :param Object value:
        """
        self.name = '@%s;%s' % (self._mktag(key, value), self.name[1:])

    def set_tags(self, keys, values):
        """

        :param keys:
        :param values:
        :return:
        """
        '@%s;%s' % (
            ';'.join(self._mktag(k, v) for k, v in zip(keys, values)),
            self.name[1:]
        )

    def average_quality(self):
        """"""
        return (
            np.mean(np.fromstring(self.quality, dtype=np.int8, count=len(self)))
            .astype(int) - 33
        )


class Reader(reader.Reader):
    """
    Fastq Reader, defines some special methods for reading and summarizing fastq data:

    :method __iter__: Iterator over fastq Record objects
    :method __len__: return number of records in file
    :method estimate_sequence_length: estimate the length of fastq sequences in file
    """

    @staticmethod
    def record_grouper(iterable):
        args = [iter(iterable)] * 4
        return zip(*args)

    def __iter__(self):
        record_type = StrRecord if self._mode == 'r' else BytesRecord
        for record in self.record_grouper(super().__iter__()):
            yield record_type(record)

    def estimate_sequence_length(self):
        """
        estimate the sequence length of a fastq file from the first 10000 records of
        the file.

        :return: int mean, float standard deviation, (np.ndarray: observed lengths,
          np.ndarray: counts per length)
        """
        i = 0
        records = iter(self)
        data = np.empty(10000, dtype=int)
        while i < 10000:
            try:
                seq = next(records).sequence
            except StopIteration:  # for fastq files shorter than 10000 records
                data = data[:i]
                break
            data[i] = len(seq) - 1  # last character is a newline
            i += 1
        return np.mean(data), np.std(data), np.unique(data, return_counts=True)


