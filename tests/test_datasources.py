import pickle

import numpy as np
import pytest

import data
from models.core.data import sources


def test_fasta_one_hot_encoder():
    for path in data.FASTA, data.FASTA_GZ:
        source = sources.fasta.OneHotEncoder(path)

        # Pickling support
        pickled = pickle.loads(pickle.dumps(source))
        assert source == pickled

        # Unknown contig
        with pytest.raises(ValueError):
            source.fetch('MT', '+', 10, 20)

        # Reverse strand is not supported
        with pytest.raises(ValueError):
            source.fetch('18', '-', 10, 20)

        # Outside the range
        with pytest.raises(ValueError):
            source.fetch('18', '+', 100000, 200000)
        with pytest.raises(ValueError):
            source.fetch('18', '+', -10, 0)

        # empty
        ohe = source.fetch('19', '+', 10, 10)
        assert ohe.shape == (4, 0)

        # 'good' paths with default mapping
        ohe = source.fetch('18', '+', 5, 15)
        # TTTTTTGAGA
        assert np.array_equal(ohe.T, [
            [0, 0, 0, 1], [0, 0, 0, 1],
            [0, 0, 0, 1], [0, 0, 0, 1],
            [0, 0, 0, 1], [0, 0, 0, 1],
            [0, 0, 1, 0], [1, 0, 0, 0],
            [0, 0, 1, 0], [1, 0, 0, 0],
        ])

        ohe = source.fetch('__all_upac__', '+', 0, 16)
        assert np.allclose(ohe.T, [
            [1, 0, 0, 0],  # A
            [0, 1, 0, 0],  # C
            [0, 0, 1, 0],  # G
            [0, 0, 0, 1],  # T
            [0, 0, 0, 1],  # U
            [0.5, 0, 0, 0.5],  # W = A|T
            [0, 0.5, 0.5, 0],  # S = C|G
            [0.5, 0.5, 0, 0],  # M = A|C
            [0, 0, 0.5, 0.5],  # K = G|T
            [0.5, 0, 0.5, 0],  # R = A|G
            [0, 0.5, 0, 0.5],  # Y = C|T
            [0, 1 / 3, 1 / 3, 1 / 3],  # B = not A
            [1 / 3, 0, 1 / 3, 1 / 3],  # D = not C
            [1 / 3, 1 / 3, 0, 1 / 3],  # H = not G
            [1 / 3, 1 / 3, 1 / 3, 0],  # V = not T
            [1 / 4, 1 / 4, 1 / 4, 1 / 4],  # N = any
        ])


def test_bed():
    for bed in data.BED, data.BED_GZ:
        source = sources.BED(bed)

        # Pickling support
        pickled = pickle.loads(pickle.dumps(source))
        assert source == pickled

        # Unknown contig / outside the range
        assert np.array_equiv(source.fetch('unknown', '.', -10, 100), 0)
        assert np.array_equiv(source.fetch('18', '.', -10, 15), 0)
        assert np.array_equiv(source.fetch('18', '.', 0, 50), 0)

        # 'good' path
        assert np.array_equiv(source.fetch('18', '+', 0, 15), 1)
        assert np.array_equiv(source.fetch('18', '+', 40, 60), 1)
        assert np.array_equiv(source.fetch('18', '+', 371, 500), 0)
        assert np.array_equiv(source.fetch('18', '+', 370, 372), [1, 0])

        expected = [0] * 8 + [1] * 782 + [0] * 63 + [1] * 583 + [0] * 224  # noqa: WPS221
        assert np.array_equiv(source.fetch('19', '-', 6140, 7800), expected)


def test_bigwig():
    # Entries:
    # bigwig.addHeader([('1', 100), ('MT', 20)])
    # bigwig.addEntries(['1', '1', '1'], [0, 35, 50], ends=[15, 50, 100], values=[0.0, 1.0, 200.0])
    # bigwig.addEntries(['MT', 'MT'], [1, 5], ends=[5, 10], values=[-2.0, 10.0])

    # Single stranded data must be passed only as unstranded
    for kwargs in {'fwd': data.BIGWIG}, {'rev': data.BIGWIG}:
        with pytest.raises(Exception):
            sources.BigWig(**kwargs)

    unstranded = sources.BigWig(unstranded=data.BIGWIG)
    stranded = sources.BigWig(fwd=data.BIGWIG, rev=data.BIGWIG)

    # Unavailable positions must be returned as zeros
    # Why? Because life is short.
    for contig in '2', '3', '12', 'MT':
        for start, end in (100, 200), (20, 25), (-5, 10):
            for strand in '+', '-', '.':
                for bw in stranded, unstranded:
                    with pytest.raises(Exception):
                        bw.fetch(contig, strand, start, end)

    for bw in stranded, unstranded:
        for strand in '+', '-':
            assert np.array_equiv(bw.fetch('MT', strand, 15, 20), 0)
            assert np.array_equiv(
                bw.fetch('MT', strand, 0, 20),
                [0.0] + [-2.0] * 4 + [10.0] * 5 + [0.0] * 10,
            )
            assert np.array_equiv(
                bw.fetch('1', strand, 10, 60),
                [0.0] * 5 + [0.0] * 20 + [1.0] * 15 + [200.0] * 10,
            )
            assert np.array_equiv(bw.fetch('1', strand, 0, 5), 0)

    # Pickling support
    for bw in stranded, unstranded:
        pickled = pickle.loads(pickle.dumps(bw))
        assert bw == pickled


def test_vocab():
    vocab = sources.fasta.Tokenizer.parse_vocab(data.VOCAB)
    assert len(vocab) == 65
    assert vocab['GCT'] == 58

    for fasta in data.FASTA, data.FASTA_GZ:
        tokenizer = sources.fasta.Tokenizer(fasta, 3, vocab)
        # Non ACGT sequences
        with pytest.raises(Exception):
            tokenizer.fetch('__all_upac__', '+', 0, 10)

        # Outside the fasta / negative strand
        workload = [('18', '-', 0, 10), ('20', '+', 500, 600), ('19', '+', -10, 23)]
        for contig, strand, start, end in workload:
            with pytest.raises(Exception):
                tokenizer.fetch(contig, strand, start, end)

        # 'good' path
        assert np.array_equiv(tokenizer.fetch('18', '+', 0, 11), [vocab['TTT']] * 9)
        # GCTGGGA
        assert np.array_equiv(tokenizer.fetch('19', '+', 33, 40), [58, 40, 32, 64, 61])
