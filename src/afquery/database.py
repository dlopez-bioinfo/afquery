import json
import sqlite3
from pathlib import Path

from .models import QueryParams, QueryResult, SampleFilter
from .preprocess.migrate import migrate_sqlite
from .query import QueryEngine


class Database:
    def __init__(self, db_path: str):
        self._path = Path(db_path)
        migrate_sqlite(self._path / "metadata.sqlite")
        self._engine = QueryEngine(db_path)
        self._manifest = json.loads((self._path / "manifest.json").read_text())

    def _reload(self) -> None:
        """Reinitialize engine and manifest after a mutating operation."""
        self._engine = QueryEngine(str(self._path))
        self._manifest = json.loads((self._path / "manifest.json").read_text())

    def _make_filter(self, phenotype, sex, tech=None) -> SampleFilter:
        return SampleFilter.parse(
            phenotype_tokens=phenotype or [],
            tech_tokens=tech or [],
            sex=sex,
        )

    def query(
        self,
        chrom: str,
        pos: int,
        phenotype: list[str] | None = None,
        sex: str = "both",
        ref: str | None = None,
        alt: str | None = None,
        tech: list[str] | None = None,
    ) -> list[QueryResult]:
        sf = self._make_filter(phenotype, sex, tech)
        params = QueryParams(chrom=chrom, pos=pos, filter=sf, ref=ref, alt=alt)
        return self._engine.query(params)

    def query_batch(
        self,
        chrom: str,
        variants: list[tuple[int, str, str]],
        phenotype: list[str] | None = None,
        sex: str = "both",
        tech: list[str] | None = None,
    ) -> list[QueryResult]:
        sf = self._make_filter(phenotype, sex, tech)
        return self._engine.query_batch(chrom, variants, sf)

    def query_region(
        self,
        chrom: str,
        start: int,
        end: int,
        phenotype: list[str] | None = None,
        sex: str = "both",
        tech: list[str] | None = None,
    ) -> list[QueryResult]:
        sf = self._make_filter(phenotype, sex, tech)
        return self._engine.query_region(chrom, start, end, sf)

    def annotate_vcf(
        self,
        input_vcf: str,
        output_vcf: str,
        phenotype: list[str] | None = None,
        sex: str = "both",
        n_workers: int | None = None,
        tech: list[str] | None = None,
    ) -> dict:
        sf = self._make_filter(phenotype, sex, tech)
        from .annotate import annotate_vcf as _annotate
        return _annotate(self._engine, input_vcf, output_vcf, sf, n_workers=n_workers)

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
        con = sqlite3.connect(str(self._path / "metadata.sqlite"))
        try:
            total = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
            by_sex: dict[str, int] = {}
            for row in con.execute(
                "SELECT sex, COUNT(*) FROM samples GROUP BY sex"
            ).fetchall():
                by_sex[row[0]] = row[1]
            by_tech: dict[str, int] = {}
            for row in con.execute(
                "SELECT t.tech_name, COUNT(*) FROM samples s"
                " JOIN technologies t ON s.tech_id = t.tech_id GROUP BY t.tech_name"
            ).fetchall():
                by_tech[row[0]] = row[1]
            by_phenotype: dict[str, int] = {}
            for row in con.execute(
                "SELECT phenotype_code, COUNT(*) FROM sample_phenotype GROUP BY phenotype_code"
            ).fetchall():
                by_phenotype[row[0]] = row[1]
            tables = {r[0] for r in con.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()}
            changelog_recent: list[dict] = []
            if "changelog" in tables:
                for r in con.execute(
                    "SELECT event_id, event_type, event_time, sample_names, notes"
                    " FROM changelog ORDER BY event_id DESC LIMIT 5"
                ).fetchall():
                    changelog_recent.append({
                        "event_id": r[0],
                        "event_type": r[1],
                        "event_time": r[2],
                        "sample_names": json.loads(r[3]) if r[3] else None,
                        "notes": r[4],
                    })
        finally:
            con.close()

        m = dict(self._manifest)
        return {
            "db_path": str(self._path),
            "db_version": m.get("db_version", "unknown"),
            "genome_build": m.get("genome_build", "unknown"),
            "schema_version": m.get("schema_version", "unknown"),
            "created_at": m.get("created_at"),
            "updated_at": m.get("updated_at"),
            "sample_count": total,
            "by_sex": by_sex,
            "by_tech": by_tech,
            "by_phenotype": by_phenotype,
            "changelog_recent": changelog_recent,
        }

    def list_samples(self) -> list[dict]:
        """Return all samples with full metadata."""
        con = sqlite3.connect(str(self._path / "metadata.sqlite"))
        try:
            rows = con.execute(
                "SELECT s.sample_id, s.sample_name, s.sex, t.tech_name,"
                " s.vcf_path, s.ingested_at"
                " FROM samples s JOIN technologies t ON s.tech_id = t.tech_id"
                " ORDER BY s.sample_id"
            ).fetchall()
            pheno_map: dict[int, list[str]] = {}
            for sid, code in con.execute(
                "SELECT sample_id, phenotype_code FROM sample_phenotype ORDER BY sample_id"
            ).fetchall():
                pheno_map.setdefault(sid, []).append(code)
        finally:
            con.close()

        return [
            {
                "sample_id": r[0],
                "sample_name": r[1],
                "sex": r[2],
                "tech": r[3],
                "phenotypes": pheno_map.get(r[0], []),
                "vcf_path": r[4],
                "ingested_at": r[5],
            }
            for r in rows
        ]

    def changelog(self, limit: int | None = None) -> list[dict]:
        """Return changelog entries ordered by event_id descending."""
        con = sqlite3.connect(str(self._path / "metadata.sqlite"))
        try:
            q = "SELECT event_id, event_type, event_time, sample_names, notes FROM changelog ORDER BY event_id DESC"
            if limit is not None:
                q += f" LIMIT {int(limit)}"
            rows = con.execute(q).fetchall()
        finally:
            con.close()

        return [
            {
                "event_id": r[0],
                "event_type": r[1],
                "event_time": r[2],
                "sample_names": json.loads(r[3]) if r[3] else None,
                "notes": r[4],
            }
            for r in rows
        ]

    def set_db_version(self, version: str) -> None:
        """Set the db_version label in manifest.json atomically."""
        import os
        from datetime import datetime, timezone
        manifest_path = self._path / "manifest.json"
        try:
            manifest = json.loads(manifest_path.read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            manifest = {}
        old_version = manifest.get("db_version", "unknown")
        manifest["db_version"] = version
        tmp_path = str(manifest_path) + ".tmp"
        with open(tmp_path, "w") as f:
            json.dump(manifest, f, indent=2)
        os.replace(tmp_path, str(manifest_path))
        self._reload()
        # Log to changelog
        con = sqlite3.connect(str(self._path / "metadata.sqlite"))
        try:
            con.execute(
                "INSERT INTO changelog (event_type, event_time, sample_names, notes) VALUES (?, ?, ?, ?)",
                (
                    "version_set",
                    datetime.now(timezone.utc).isoformat(),
                    None,
                    f"version changed from {old_version} to {version}",
                ),
            )
            con.commit()
        finally:
            con.close()

    def get_all_phenotypes(self) -> list[str]:
        """Return all available phenotype codes in the database."""
        import sqlite3
        con = sqlite3.connect(str(self._path / "metadata.sqlite"))
        phenotypes = [
            row[0] for row in con.execute(
                "SELECT DISTINCT phenotype_code FROM sample_phenotype ORDER BY phenotype_code"
            ).fetchall()
        ]
        con.close()
        return phenotypes
