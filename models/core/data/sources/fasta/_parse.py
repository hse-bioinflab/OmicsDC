import gzip
from pathlib import Path

from Bio import SeqIO


def parse(fasta: Path) -> dict[str, str]:
    sequences = {}
    if fasta.name.endswith('gz'):
        stream = gzip.open(fasta, 'rt')
    else:
        stream = open(fasta, 'rt')  # noqa: WPS515

    for contig in SeqIO.parse(stream, 'fasta'):
        sequences[contig.id] = str(contig.seq).upper()

    stream.close()
    return sequences
