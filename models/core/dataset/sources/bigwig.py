from pathlib import Path
from typing import Literal, Optional

import numpy as np
import pyBigWig

from .data_source import DataSource


class BigWig(DataSource):
    def __init__(
        self, *, naval: float = 0,
        fwd: Optional[Path] = None,
        rev: Optional[Path] = None,
        unstranded: Optional[Path] = None,
    ):
        self.bws = {}
        self.paths = {}
        self.naval = naval

        match (fwd, rev, unstranded):
            case (None, None, _) if unstranded is not None:
                self.paths = {"+": unstranded, "-": unstranded, ".": unstranded}
            case (_, _, None) if fwd is not None and rev is not None:
                self.paths = {"+": fwd, "-": rev}
            case _:
                raise ValueError("Data must be either unstranded (fwd == rev == None) or stranded (unstranded == None)")

    def fetch(self, contig: str, strand: Literal["+", "-", "."], start: int, end: int) -> np.ndarray:
        if strand not in self.paths:
            raise ValueError(f"Trying to fetch data for unavailable strand: {strand}")
        if strand not in self.bws:
            self.bws[strand] = pyBigWig.open(self.paths[strand].as_posix())

        values = self.bws[strand].values(contig, start, end, numpy=True)
        values[np.isnan(values)] = self.naval
        return values

    def __getstate__(self):
        return self.paths, self.naval

    def __setstate__(self, state):
        self.paths = state[0]
        self.naval = state[1]
        self.bws = {}

    def __del__(self):
        for bw in self.bws.values():
            bw.close()
        self.bws.clear()

    def __eq__(self, other):
        return isinstance(other, BigWig) \
            and other.naval == self.naval \
            and other.paths == self.paths

    __hash__ = None