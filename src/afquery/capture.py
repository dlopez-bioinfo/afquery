import pickle
import warnings
import pyranges as pr
import pandas as pd
from .models import Technology


class CaptureIndex:
    _always_covered: bool
    _pr: "pr.PyRanges | None"

    def __init__(self, *, always_covered: bool = False, pyranges_obj=None):
        self._always_covered = always_covered
        self._pr = pyranges_obj

    @classmethod
    def wgs(cls) -> "CaptureIndex":
        return cls(always_covered=True)

    @classmethod
    def from_bed(cls, bed_path: str) -> "CaptureIndex":
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)
            ranges = pr.read_bed(bed_path)
        return cls(pyranges_obj=ranges)

    def covers(self, chrom: str, pos: int) -> bool:
        """Return True if 1-based pos is within any region for this chrom.

        BED is 0-based half-open: [Start, End) → 1-based: Start < pos <= End
        """
        if self._always_covered:
            return True
        df = self._pr.df
        chrom_df = df[df["Chromosome"] == chrom]
        if chrom_df.empty:
            return False
        return bool(((chrom_df["Start"] < pos) & (pos <= chrom_df["End"])).any())

    def save(self, path: str) -> None:
        with open(path, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, path: str) -> "CaptureIndex":
        with open(path, "rb") as f:
            return pickle.load(f)


def load_capture_indices(
    technologies: list[Technology], capture_dir: str
) -> dict[int, "CaptureIndex"]:
    result = {}
    for tech in technologies:
        path = f"{capture_dir}/tech_{tech.tech_id}.pickle"
        result[tech.tech_id] = CaptureIndex.load(path)
    return result
