from ..capture import CaptureIndex
from ..models import Technology


def build_capture_indices(
    technologies: list[Technology],
    capture_dir: str,
) -> None:
    for tech in technologies:
        if tech.bed_path is None:
            idx = CaptureIndex.wgs()
        else:
            idx = CaptureIndex.from_bed(tech.bed_path)
        idx.save(f"{capture_dir}/tech_{tech.tech_id}.pickle")
