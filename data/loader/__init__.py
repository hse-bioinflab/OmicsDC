from pathlib import Path
from loader.omicDC import omics


def assembly(tag: str, saveto: Path, *_, force: bool = False):
    supported = {"GRCh38", "GRCm39"}
    assert tag in supported, f"Requested assembly ({tag}) is not among supported: {','.join(supported)}"

    # TODO: clearer error message
    assert saveto.name.endswith(".gz"), "Loaded assemblies must be saved as indexed & bgzip-ed files"

    # Use GENCODE:
    # After downloading, unzip -> bgzip -> samtools faidx files
    # GRCh38 - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
    # GRCm39 - https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.primary_assembly.genome.fa.gz
    raise NotImplementedError()


__all__ = ["omics", "assembly"]
