from collections import defaultdict
from pathlib import Path
from typing import Optional, Literal

import intervaltree as it
import numpy as np
from pybedtools import BedTool

from .data_source import DataSource


class BED(DataSource):
    def __init__(self, bed: Path, usecol: Optional[int] = None):
        if not bed.is_file():
            raise ValueError(f"BED file doesn't exist: {bed}")
        self.bed = bed
        self.col = usecol

        self.index = defaultdict(lambda *args: it.IntervalTree())
        for win in BedTool(bed):
            if usecol is None:
                score = 1
            else:
                score = float(win.fields[usecol])
            self.index[(win.chrom, win.strand)].addi(win.start, win.end, score)

        self.index = dict(self.index)

    def fetch(self, contig: str, strand: Literal['+', '-', '.'], start: int, end: int) -> np.ndarray:
        if start > end:
            raise ValueError(f'Start must be <= end, got {start} > {end}')

        result = np.zeros(end - start, dtype=np.float32)
        if (contig, strand) not in self.index:
            return result

        index = self.index[(contig, strand)]
        for hit in index.overlap(start, end):
            hstart, hend = max(hit.begin, start), min(hit.end, end)
            assert start <= hstart <= hend <= end  # noqa: S101
            result[hstart - start: hend - start] = hit.data
        return result

    def __eq__(self, other):
        return isinstance(other, BED) \
            and self.bed == other.bed \
            and self.col == other.col \
            and other.index == self.index

    __hash__ = None
