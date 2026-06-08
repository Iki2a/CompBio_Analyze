from pathlib import Path
from urllib.request import urlopen
import tarfile
import tempfile

import pandas as pd
import rdata


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_DIR = PROJECT_ROOT / "data" / "external" / "gct"
SEQINR_URL = "https://cran.r-project.org/src/contrib/Archive/seqinr/seqinr_4.2-36.tar.gz"


def classify_temperature(topt):
    return "Termofil" if topt > 45 else "Non-termofil"


def clean_species_table(df):
    out = df.copy()
    out["species"] = (
        out["Genus"].astype(str).str.title()
        + " "
        + out["species"].astype(str).str.lower()
    )
    out = out.rename(
        columns={
            "Genus": "genus",
            "GC": "GCw",
        }
    )
    out["species_epithet"] = out["species"].str.split(" ", n=1).str[1]
    out["thermal_group"] = out["Topt"].apply(classify_temperature)
    return out[["species", "genus", "species_epithet", "GCw", "Topt", "thermal_group"]]


def clean_rna_table(df, source_name):
    out = df.copy()
    out = out.rename(columns={"Genus": "genus", "GC": "stem_gc_fraction"})
    out["source_table"] = source_name
    return out[["source_table", "genus", "Topt", "stem_gc_fraction"]]


def load_gct_from_seqinr():
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        archive_path = tmp_dir / "seqinr.tar.gz"
        archive_path.write_bytes(urlopen(SEQINR_URL, timeout=60).read())

        with tarfile.open(archive_path, "r:gz") as archive:
            member = next(m for m in archive.getmembers() if m.name.endswith("data/gcT.rda"))
            rda_path = tmp_dir / "gcT.rda"
            with archive.extractfile(member) as source:
                rda_path.write_bytes(source.read())

        parsed = rdata.parser.parse_file(rda_path)
        converted = rdata.conversion.convert(parsed)
        return converted["gcT"]


def write_metadata():
    metadata = """# gcT External Validation Dataset

Source: `gcT` dataset from the R package `seqinr` version 4.2-36.

Primary citation:
Galtier, N. & Lobry, J.R. (1997). Relationships between genomic G+C content,
RNA secondary structures, and optimal growth temperature in prokaryotes.
Journal of Molecular Evolution 44:632-636.

Dataset notes:
- `gct_species.csv` contains 772 bacterial species with optimal growth
  temperature (`Topt`) and genomic G+C content (`GCw`).
- RNA tables contain stem G+C fractions from the original `gcT` object and are
  not equivalent to the locally extracted `GC_tRNA`, `GC_16S`, or `GC_23S`
  features used in the TEMPURA/NCBI pipeline.
- This dataset is used only as external validation and is not merged into
  `data/processed/real_gc_data.csv`.

Download source:
https://cran.r-project.org/src/contrib/Archive/seqinr/seqinr_4.2-36.tar.gz
"""
    (EXTERNAL_DIR / "gct_metadata.md").write_text(metadata, encoding="utf-8")


def main():
    EXTERNAL_DIR.mkdir(parents=True, exist_ok=True)
    gct = load_gct_from_seqinr()

    species = clean_species_table(gct["species"])
    species.to_csv(EXTERNAL_DIR / "gct_species.csv", index=False)

    for table_name in ["gc16S", "gctRNA", "gc23S", "gc5S"]:
        clean_rna_table(gct[table_name], table_name).to_csv(
            EXTERNAL_DIR / f"{table_name}.csv", index=False
        )

    write_metadata()

    print(f"Saved gcT species table: {EXTERNAL_DIR / 'gct_species.csv'}")
    print(f"Rows: {len(species)}")
    print("Saved RNA stem GC tables: gc16S.csv, gctRNA.csv, gc23S.csv, gc5S.csv")


if __name__ == "__main__":
    main()
