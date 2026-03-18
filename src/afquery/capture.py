import bisect
import pickle
import warnings
import pyranges as pr
import pandas as pd
from .models import Technology


class CaptureIndex:
    _always_covered: bool
    _pr: "pr.PyRanges | None"
    # Per-chrom sorted interval index: {chrom: (sorted_starts, corresponding_ends)}
    _index: "dict[str, tuple[list[int], list[int]]]"

    def __init__(self, *, always_covered: bool = False, pyranges_obj=None):
        self._always_covered = always_covered
        self._pr = pyranges_obj
        self._index = {}
        if pyranges_obj is not None:
            df = pyranges_obj.df
            for chrom, group in df.groupby("Chromosome", observed=True):
                pairs = sorted(zip(group["Start"].tolist(), group["End"].tolist()))
                self._index[str(chrom)] = (
                    [s for s, _ in pairs],
                    [e for _, e in pairs],
                )

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
        entry = self._index.get(chrom)
        if entry is None:
            return False
        starts, ends = entry
        # Find rightmost interval whose Start < pos (0-based start, so Start < pos_1based)
        # bisect_left gives insertion point for pos in starts; all starts[:idx] < pos
        idx = bisect.bisect_left(starts, pos) - 1
        # Check candidates from idx downward: intervals may overlap, so we scan
        # backwards until Start is so small it can't possibly reach pos
        while idx >= 0:
            if ends[idx] >= pos:
                return True
            # Since starts are sorted ascending, once we go far enough back,
            # no remaining interval can cover pos either
            idx -= 1
        return False

    def __setstate__(self, state: dict) -> None:
        """Rebuild _index when loading pickles saved before the binary-search refactor."""
        self.__dict__.update(state)
        if not hasattr(self, "_index"):
            self._index = {}
            if getattr(self, "_pr", None) is not None:
                df = self._pr.df
                for chrom, group in df.groupby("Chromosome", observed=True):
                    pairs = sorted(zip(group["Start"].tolist(), group["End"].tolist()))
                    self._index[str(chrom)] = (
                        [s for s, _ in pairs],
                        [e for _, e in pairs],
                    )

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
