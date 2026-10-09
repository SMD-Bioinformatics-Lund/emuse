import pandas as pd


def read_assignment_stats(path):
    assignment = pd.read_csv(path, sep="\t")
    summary = assignment.iloc[:, 1:].agg(["median", "mean"]).T.reset_index()
    summary.columns = ["tax_id", "median_probability", "mean_probability"]
    summary["tax_id"] = summary["tax_id"].astype(str)
    return summary
