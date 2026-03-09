import json
from pathlib import Path

from .models import QueryParams, QueryResult
from .query import QueryEngine


class Database:
    def __init__(self, db_path: str):
        self._path = Path(db_path)
        self._engine = QueryEngine(db_path)
        self._manifest = json.loads((self._path / "manifest.json").read_text())

    def _reload(self) -> None:
        """Reinitialize engine and manifest after a mutating operation."""
        self._engine = QueryEngine(str(self._path))
        self._manifest = json.loads((self._path / "manifest.json").read_text())

    def query(
        self,
        chrom: str,
        pos: int,
        icd10: list[str],
        sex: str = "both",
    ) -> list[QueryResult]:
        params = QueryParams(chrom=chrom, pos=pos, icd10_codes=icd10, sex_filter=sex)
        return self._engine.query(params)

    def query_batch(
        self,
        chrom: str,
        positions: list[int],
        icd10: list[str],
        sex: str = "both",
    ) -> list[QueryResult]:
        return self._engine.query_batch(chrom, positions, icd10, sex)

    def query_region(
        self,
        chrom: str,
        start: int,
        end: int,
        icd10: list[str],
        sex: str = "both",
    ) -> list[QueryResult]:
        return self._engine.query_region(chrom, start, end, icd10, sex)

    def annotate_vcf(
        self,
        input_vcf: str,
        output_vcf: str,
        icd10: list[str],
        sex: str = "both",
    ) -> dict:
        from .annotate import annotate_vcf as _annotate
        return _annotate(self._engine, input_vcf, output_vcf, icd10, sex)

    def add_samples(
        self,
        manifest_path: str,
        threads: int = 8,
        tmp_dir: str | None = None,
        bed_dir: str | None = None,
        genome_build: str | None = None,
    ) -> dict:
        from .preprocess.update import add_samples as _add
        result = _add(
            str(self._path), manifest_path, threads=threads,
            tmp_dir=tmp_dir, bed_dir=bed_dir, genome_build=genome_build,
        )
        self._reload()
        return result

    def remove_samples(self, sample_names: list[str]) -> dict:
        from .preprocess.update import remove_samples as _remove
        result = _remove(str(self._path), sample_names)
        self._reload()
        return result

    def compact(self) -> dict:
        from .preprocess.compact import compact_database
        result = compact_database(self._path)
        self._reload()
        return result

    def check(self) -> list:
        from .preprocess.update import check_database
        return check_database(str(self._path))

    def info(self) -> dict:
        return dict(self._manifest)
