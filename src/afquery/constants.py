# PAR1 and PAR2 for GRCh37 and GRCh38
PAR = {
    "GRCh38": {
        "chrX": [(10_001, 2_781_479), (155_701_383, 156_030_895)],
        "chrY": [(10_001, 2_781_479), (56_887_903, 57_217_415)],
    },
    "GRCh37": {
        "chrX": [(60_001, 2_699_520), (154_931_044, 155_260_560)],
        "chrY": [(10_001, 2_649_520), (59_034_050, 59_363_566)],
    },
}

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
SEX_CHROMS = ["chrX", "chrY"]
MITO_CHROM = "chrM"
ALL_CHROMS = AUTOSOMES + SEX_CHROMS + [MITO_CHROM]
CHROM_ORDER: dict[str, int] = {c: i for i, c in enumerate(ALL_CHROMS)}
VALID_GENOME_BUILDS = frozenset(["GRCh37", "GRCh38"])
VALID_SEX = frozenset(["male", "female"])


def normalize_chrom(chrom: str) -> str:
    """Normalize chromosome name to 'chr'-prefixed canonical form.

    '1' -> 'chr1', 'chr1' -> 'chr1', 'X' -> 'chrX', 'MT' -> 'chrM'
    """
    chrom = chrom.strip()
    if chrom.upper() == "MT":
        return "chrM"
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return chrom


def is_autosome(chrom: str) -> bool:
    return normalize_chrom(chrom) in AUTOSOMES


def is_sex_chrom(chrom: str) -> bool:
    return normalize_chrom(chrom) in SEX_CHROMS


def is_mito(chrom: str) -> bool:
    return normalize_chrom(chrom) == MITO_CHROM
