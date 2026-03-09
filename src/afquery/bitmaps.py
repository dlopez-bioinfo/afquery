from pyroaring import BitMap
from .models import Sample


def serialize(bm: BitMap) -> bytes:
    return bm.serialize()


def deserialize(data: bytes) -> BitMap:
    return BitMap.deserialize(data)


def build_sex_bitmaps(samples: list[Sample]) -> dict[str, BitMap]:
    """Returns {'male': BitMap([...]), 'female': BitMap([...])}."""
    male_ids = [s.sample_id for s in samples if s.sex == "male"]
    female_ids = [s.sample_id for s in samples if s.sex == "female"]
    return {
        "male": BitMap(male_ids),
        "female": BitMap(female_ids),
    }


def build_icd10_bitmaps(
    sample_icd10: list[tuple[int, str]]  # list of (sample_id, icd10_code)
) -> dict[str, BitMap]:
    """Returns {icd10_code: BitMap([sample_ids...])}."""
    result: dict[str, list[int]] = {}
    for sample_id, code in sample_icd10:
        result.setdefault(code, []).append(sample_id)
    return {code: BitMap(ids) for code, ids in result.items()}


def build_tech_bitmaps(samples: list[Sample]) -> dict[int, BitMap]:
    """Returns {tech_id: BitMap([sample_ids...])}."""
    result: dict[int, list[int]] = {}
    for s in samples:
        result.setdefault(s.tech_id, []).append(s.sample_id)
    return {tech_id: BitMap(ids) for tech_id, ids in result.items()}
