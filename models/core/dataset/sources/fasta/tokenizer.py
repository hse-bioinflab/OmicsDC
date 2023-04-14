from pathlib import Path
from typing import Literal

import numpy as np

from ._parse import parse as _parse
from ..data_source import DataSource


class Tokenizer(DataSource):
    def __init__(self, fasta: Path, kmers: int, vocab: dict[str, int]):
        if kmers < 1:
            raise ValueError(f"Kmers must be >= 1, found {kmers}")
        self.fasta = fasta
        self.sequences = _parse(fasta)
        self.kmers = kmers
        self.kmer2id = vocab

    @staticmethod
    def parse_vocab(config: Path):
        if not config.is_file():
            raise ValueError(f"Vocab file doesn't exist: {config}")
        with open(config) as stream:
            config = stream.read().upper()
        tokens = config.split("\n")  # [:-1]
        return {token.strip("\\"): ind for ind, token in enumerate(tokens)}

    def pad(self, sequence: np.ndarray, tolen: int, padtok: str = "[PAD]") -> np.ndarray:
        if len(sequence) > tolen:
            raise ValueError("Padding is only possible when current length < target")
        if len(sequence) < tolen:
            pad = self.kmer2id[padtok]
            sequence = np.pad(sequence, (0, len(sequence) - tolen), "constant", constant_values=(pad, pad))
        assert len(sequence) == tolen  # noqa: S101
        return sequence

    def fetch(self, contig: str, strand: Literal["+", "-", "."], start: int, end: int) -> np.ndarray:
        if strand not in {"+", "."}:
            raise ValueError("Reverse strand is not supported")

        sequence = self.sequences[contig]
        seqlen = len(sequence)
        is_within_sequence = 0 <= start <= end <= seqlen
        if not is_within_sequence:
            coords = f"{contig}:{start}-{end}"
            raise ValueError(
                f"Requested intervals is outside of the sequence with len {seqlen}: {coords}",
            )
        sequence = sequence[start: end]

        # split sequence into kmers & tokenize
        sequence = [sequence[ind: ind + self.kmers] for ind in range(len(sequence) - self.kmers + 1)]
        sequence = [self.kmer2id[kmer] for kmer in sequence]

        return np.asarray(sequence, dtype=np.int32)
