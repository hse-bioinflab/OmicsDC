from pathlib import Path


def omics(expid: str, assembly: str, antigen_class: str, antigen: str, cell_type, cell: str, storage: Path):
    # TODO: think about parameters naming, order, etc & write docs for each parameter
    raise NotImplementedError()


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
