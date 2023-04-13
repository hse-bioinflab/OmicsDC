from pathlib import Path
from typing import Optional, Literal

import numpy as np

from ._parse import parse as _parse
from ..data_source import DataSource

DEFAULT_ONE_HOT_ENCODING = (
    ('A', [0]),
    ('C', [1]),
    ('G', [2]),
    ('T', [3]),
    ('U', [3]),
    ('M', [0, 1]),
    ('R', [0, 2]),
    ('W', [0, 3]),
    ('S', [1, 2]),
    ('Y', [1, 3]),
    ('K', [2, 3]),
    ('V', [0, 1, 2]),
    ('H', [0, 1, 3]),
    ('D', [0, 2, 3]),
    ('B', [1, 2, 3]),
    ('X', [0, 1, 2, 3]),
    ('N', [0, 1, 2, 3]),
)

NucleotidesOHE = dict[str, list[int]]


class OneHotEncoder(DataSource):
    def __init__(self, fasta: Path, mapping: Optional[NucleotidesOHE] = None):
        if not fasta.is_file():
            raise ValueError(f"Fasta doesn't exist: {fasta}")

        self.fasta = fasta
        self.sequences = _parse(fasta)

        # Prepare the mapping matrix
        self.mapping = self._prepare_mapping(mapping)

    def fetch(self, contig: str, strand: Literal['+', '-', '.'], start: int, end: int) -> np.ndarray:
        if strand not in {'+', '.'}:
            raise ValueError('Reverse strand is not supported')

        if contig not in self.sequences:
            raise ValueError(f'Contig {contig} was not present in the {self.fasta} fasta file.')

        sequence = self.sequences[contig]
        seqlen = len(sequence)
        is_within_sequence = 0 <= start <= end <= seqlen
        if not is_within_sequence:
            error = f'Requested intervals is outside of the sequence with len {seqlen}: {contig}:{start}-{end}'
            raise ValueError(error)

        sequence = sequence[start: end]
        sequence = np.fromiter(sequence.encode('ASCII'), dtype=np.uint8, count=(end - start))
        # One-hot-encoding
        encoded = self.mapping[:, sequence]
        return encoded

    def __eq__(self, other):
        return isinstance(other, OneHotEncoder) \
            and np.array_equal(other.mapping, self.mapping) \
            and other.fasta == self.fasta \
            and other.sequences == self.sequences

    __hash__ = None

    def _prepare_mapping(self, mapping: Optional[NucleotidesOHE]) -> np.ndarray:
        if mapping is None:
            mapping = dict(DEFAULT_ONE_HOT_ENCODING)

        ohe = np.zeros((4, 255), dtype=np.float32)
        for letter, indices in mapping.items():
            if not isinstance(letter, str) or len(letter) != 1:
                raise ValueError('Incorrect mapping format')
            upper, lower = ord(letter.upper()), ord(letter.lower())
            weight = 1 / len(indices)
            for ind in indices:
                ohe[ind, upper] = weight
                ohe[ind, lower] = weight
        return ohe
