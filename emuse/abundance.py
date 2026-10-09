import pandas as pd

TABLE_COLUMNS = ["abundance", "species", "genus", "family", "tax_id", "estimated_counts"]


def read_rel_abundance(path):
    df = pd.read_csv(path, sep="\t", dtype={"tax_id": str})
    df = df.rename(columns={"estimated counts": "estimated_counts"})
    missing = [column for column in TABLE_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(missing)}")
    return df.sort_values(by="abundance", ascending=False).reset_index(drop=True)


def abundance_table(df):
    return df[TABLE_COLUMNS].copy()


def to_records(df):
    return df.astype(object).where(df.notna(), None).to_dict(orient="records")
