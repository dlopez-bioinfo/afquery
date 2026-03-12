import sqlite3
from pathlib import Path


def migrate_sqlite(db_path: str | Path) -> None:
    """Idempotently apply schema migrations to metadata.sqlite."""
    con = sqlite3.connect(str(db_path))
    try:
        col_info = con.execute("PRAGMA table_info(samples)").fetchall()
        existing_cols = {row[1] for row in col_info}

        if "vcf_path" not in existing_cols:
            con.execute("ALTER TABLE samples ADD COLUMN vcf_path TEXT")
        if "ingested_at" not in existing_cols:
            con.execute("ALTER TABLE samples ADD COLUMN ingested_at TEXT")

        con.execute("""
            CREATE TABLE IF NOT EXISTS changelog (
                event_id     INTEGER PRIMARY KEY AUTOINCREMENT,
                event_type   TEXT NOT NULL,
                event_time   TEXT NOT NULL,
                sample_names TEXT,
                notes        TEXT
            )
        """)
        con.commit()
    finally:
        con.close()
