import os
from dataclasses import dataclass, field


@dataclass
class ParsedSample:
    sample_name: str
    sex: str               # 'male' | 'female'
    tech_name: str
    vcf_path: str          # absolute path
    icd10_codes: list[str]


@dataclass
class ParsedTechnology:
    tech_name: str
    bed_path: str | None   # None = WGS


class ManifestError(ValueError):
    pass


def parse_manifest(
    manifest_path: str,
    bed_dir: str | None = None,
) -> tuple[list[ParsedSample], list[ParsedTechnology]]:
    manifest_dir = os.path.dirname(os.path.abspath(manifest_path))

    with open(manifest_path) as f:
        lines = f.read().splitlines()

    if not lines:
        raise ManifestError("Manifest is empty")

    header = lines[0].split("\t")
    required_cols = {"sample_name", "sex", "tech_name", "vcf_path", "icd10_codes"}
    missing = required_cols - set(header)
    if missing:
        raise ManifestError(f"Missing required columns: {', '.join(sorted(missing))}")

    col_idx = {col: i for i, col in enumerate(header)}

    samples: list[ParsedSample] = []
    seen_names: set[str] = set()
    tech_map: dict[str, ParsedTechnology] = {}
    tech_order: list[str] = []

    for line_num, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        parts = line.split("\t")

        sample_name = parts[col_idx["sample_name"]].strip()
        sex = parts[col_idx["sex"]].strip()
        tech_name = parts[col_idx["tech_name"]].strip()
        vcf_path_raw = parts[col_idx["vcf_path"]].strip()
        icd10_raw = parts[col_idx["icd10_codes"]].strip()

        if sex not in ("male", "female"):
            raise ManifestError(
                f"Line {line_num}: invalid sex '{sex}', must be 'male' or 'female'"
            )

        icd10_codes = [c.strip() for c in icd10_raw.split(",") if c.strip()]
        if not icd10_codes:
            raise ManifestError(
                f"Line {line_num}: icd10_codes is empty for sample '{sample_name}'"
            )

        if os.path.isabs(vcf_path_raw):
            vcf_path = vcf_path_raw
        else:
            vcf_path = os.path.join(manifest_dir, vcf_path_raw)
        if not os.path.exists(vcf_path):
            raise ManifestError(
                f"Line {line_num}: vcf_path does not exist: {vcf_path}"
            )

        if sample_name in seen_names:
            raise ManifestError(f"Duplicate sample_name: '{sample_name}'")
        seen_names.add(sample_name)

        if tech_name not in tech_map:
            if tech_name.upper() == "WGS":
                bed_path = None
            else:
                if bed_dir is None:
                    raise ManifestError(
                        f"bed_dir required for non-WGS technology '{tech_name}'"
                    )
                bed_path = os.path.join(bed_dir, f"{tech_name}.bed")
                if not os.path.exists(bed_path):
                    raise ManifestError(
                        f"BED file not found for technology '{tech_name}': {bed_path}"
                    )
            tech_map[tech_name] = ParsedTechnology(tech_name=tech_name, bed_path=bed_path)
            tech_order.append(tech_name)

        samples.append(ParsedSample(
            sample_name=sample_name,
            sex=sex,
            tech_name=tech_name,
            vcf_path=vcf_path,
            icd10_codes=icd10_codes,
        ))

    technologies = [tech_map[name] for name in tech_order]
    return samples, technologies
