import gzip
from itertools import chain
from pathlib import Path
from typing import Union, Optional

from Bio import SeqIO
from joblib import Parallel, delayed
from pybedtools import Interval, BedTool, chromsizes as fetch_chromsizes

from .typing import BedLike


def genomic_windows(
    chromsizes: Union[str, dict[str, int]], winsize: int,
    step: Optional[int] = None, exclude: Optional[BedLike] = None,
    dropshort: bool = True,
) -> BedTool:
    if isinstance(chromsizes, str):
        chromsizes = fetch_chromsizes(chromsizes)
        chromsizes = {contig: end for contig, (_, end) in chromsizes.items()}

    assembly = BedTool([Interval(contig, 0, length) for contig, length in chromsizes.items()]).sort()

    # Exclude some regions if needed
    if exclude is not None:
        if isinstance(exclude, list):
            exclude = BedTool(exclude)
        exclude = exclude.sort()
        assembly = assembly.subtract(exclude)

    step = winsize if step is None else step
    windows = assembly.window_maker(b=assembly, w=winsize, s=step)

    if dropshort:
        windows = windows.filter(lambda wd: wd.length == winsize)

    return windows.sort()


def ambiguous_sites(
    fasta: Path,
    allowednuc: tuple[str, ...] = ('A', 'a', 'C', 'c', 'G', 'g', 'T', 't'),
    n_jobs: int = -1,
) -> list[Interval]:
    def job(contig, seq, allowed) -> list[Interval]:  # noqa: WPS430
        ambcursor = None
        ambiguous = []

        for ind, letter in enumerate(seq):
            if letter in allowed:
                if ambcursor is not None:
                    # Save currently tracked ambiguous regions
                    ambiguous.append(Interval(contig, ambcursor, ind))
                    ambcursor = None
            elif ambcursor is None:
                # Start tracking ambiguous regions
                ambcursor = ind

        if ambcursor is not None:
            ambiguous.append(Interval(contig, ambcursor, len(seq)))
        return ambiguous

    if not fasta.is_file():
        raise ValueError(f"Fasta {fasta} doesn't exist.")

    allowednuc = set(allowednuc)

    if fasta.name.endswith('.gz'):
        stream = gzip.open(fasta, 'rt')
    else:
        stream = open(fasta, 'r')   # noqa: WPS515

    fasta = SeqIO.parse(stream, 'fasta')
    intervals = Parallel(n_jobs=n_jobs)(delayed(job)(seq.id, seq.seq, allowednuc) for seq in fasta)
    stream.close()

    return BedTool(list(chain(*intervals))).sort()
