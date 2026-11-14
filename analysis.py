"""
FGA vs immune gene expression in TCGA-PAAD (cBioPortal paad_tcga_pan_can_atlas_2018).

Usage:
  python analysis.py --clinical clinical.tsv --expr expression.txt [--purity abs.txt --purity-id-col sample --purity-col purity] [--leuk leuk.tsv]

Purity file: any CSV/TSV with a TCGA sample barcode column (default 'Sample ID')
and a purity column (e.g. 'CPE' from Aran et al. 2015, Supp Table). Merged on the
first 15 characters of the barcode (TCGA-XX-XXXX-01).
Outputs CSVs to --out (default ./results).
"""
import argparse, os
import numpy as np, pandas as pd
from scipy import stats

GENES = ["CD274", "PDCD1", "CTLA4", "CD8A", "GZMB", "PRF1", "LAG3", "HAVCR2"]


def load_expr(path):
    """cBioPortal export that may have lost newlines: 10 rows x 185 comma tokens
    (label + 184 samples). Falls back to normal parsing if newlines exist."""
    s = open(path).read()
    if "\n" in s.strip():
        df = pd.read_csv(path, header=None, index_col=0, sep=None, engine="python").T
        return df
    t = s.strip().split(",")
    n = 185
    rows = [t[i:i + n] for i in range(0, len(t), n)]
    d = {r[0]: r[1:] for r in rows}
    df = pd.DataFrame({k: d[k] for k in d if k not in ("STUDY_ID", "SAMPLE_ID")})
    df.index = d["SAMPLE_ID"]
    return df.apply(pd.to_numeric, errors="coerce")


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q); out[o] = np.minimum(q, 1)
    return out


def partial_spearman(x, y, covs):
    """Rank-based partial Spearman: rank everything, residualize on covariates, Pearson."""
    d = pd.concat([x, y] + covs, axis=1).dropna()
    r = d.rank()
    Z = np.column_stack([np.ones(len(r))] + [r.iloc[:, i] for i in range(2, r.shape[1])])
    res = [r.iloc[:, i] - Z @ np.linalg.lstsq(Z, r.iloc[:, i], rcond=None)[0] for i in (0, 1)]
    rho = np.corrcoef(res[0], res[1])[0, 1]
    df = len(r) - 2 - len(covs)
    tstat = rho * np.sqrt(df / (1 - rho**2))
    return rho, 2 * stats.t.sf(abs(tstat), df), len(r)


def corr_table(df, xcol, cols, covs=None, label=""):
    rows = []
    for g in cols:
        if covs:
            rho, p, n = partial_spearman(df[xcol], df[g], [df[c] for c in covs])
        else:
            d = df[[xcol, g]].dropna()
            rho, p = stats.spearmanr(d[xcol], d[g]); n = len(d)
        rows.append(dict(analysis=label, gene=g, rho=rho, p=p, n=n))
    t = pd.DataFrame(rows)
    t["q"] = bh(t["p"])
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinical", required=True)
    ap.add_argument("--expr", required=True)
    ap.add_argument("--purity")
    ap.add_argument("--purity-col", default="CPE")
    ap.add_argument("--purity-id-col", default="Sample ID")
    ap.add_argument("--leuk", help="TCGA_all_leuk_estimate.masked.20170107.tsv (no header: cancer, barcode, fraction)")
    ap.add_argument("--leuk-cancer", default="PAAD")
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    clin = pd.read_csv(a.clinical, sep="\t").rename(columns={
        "Sample ID": "sample", "Fraction Genome Altered": "FGA", "Mutation Count": "mut",
        "TMB (nonsynonymous)": "TMB", "Overall Survival (Months)": "os_m",
        "Overall Survival Status": "os_s", "Diagnosis Age": "age", "Subtype": "subtype",
        "Neoplasm Disease Stage American Joint Committee on Cancer Code": "stage", "Sex": "sex"})
    expr = load_expr(a.expr)
    df = clin.merge(expr, left_on="sample", right_index=True, how="left")
    df["CYT"] = np.sqrt(df["GZMB"] * df["PRF1"])  # geometric mean of 2 values
    cols = GENES + ["CYT"]
    df["annotated"] = df["subtype"].notna()
    df = df[df[GENES].notna().all(axis=1)].copy()
    print(f"expression cohort n={len(df)}")

    # --- purity merge ---
    if a.purity:
        pu = pd.read_csv(a.purity, sep=None, engine="python")
        pu["key"] = pu[a.purity_id_col].astype(str).str[:15]
        df["key"] = df["sample"].str[:15]
        df = df.merge(pu[["key", a.purity_col]].rename(columns={a.purity_col: "purity"}).drop_duplicates("key"),
                      on="key", how="left")
        print(f"purity available for {df['purity'].notna().sum()} / {len(df)}")
    else:
        print("No purity file given: skipping purity-adjusted analyses.")

    if a.leuk:
        lf = pd.read_csv(a.leuk, sep="\t", header=None, names=["cancer", "barcode", "LF"])
        lf = lf[lf.cancer == a.leuk_cancer].copy()
        lf["key"] = lf["barcode"].str[:15]
        if "key" not in df: df["key"] = df["sample"].str[:15]
        df = df.merge(lf.drop_duplicates("key")[["key", "LF"]], on="key", how="left")
        print(f"leukocyte fraction available for {df['LF'].notna().sum()} / {len(df)}")
        lrows = []
        for lab, sub in [("all", df), ("PDAC-annotated", df[df.annotated]),
                         ("FGA>0.02", df[df.FGA > 0.02]),
                         ("annotated & FGA>0.02", df[df.annotated & (df.FGA > 0.02)])]:
            d = sub[["FGA", "LF"]].dropna()
            r, p = stats.spearmanr(d.FGA, d.LF)
            lrows.append(dict(set=lab, rho=r, p=p, n=len(d)))
        lt = pd.DataFrame(lrows)
        lt.to_csv(f"{a.out}/fga_vs_leukocyte_fraction.csv", index=False)
        print("FGA vs leukocyte fraction:\n", lt.round(3).to_string(index=False))

    res = []
    # published analyses
    res.append(corr_table(df, "FGA", cols, label="all"))
    res.append(corr_table(df[df.annotated], "FGA", cols, label="PDAC-annotated"))
    # threshold robustness
    for t in (0, 0.005, 0.01, 0.02, 0.05, 0.10):
        res.append(corr_table(df[df.FGA > t], "FGA", cols, label=f"FGA>{t}"))
    res.append(corr_table(df, "FGA", cols, covs=["mut"], label="partial|mut"))

    if "purity" in df:
        dp = df.dropna(subset=["purity"])
        res.append(corr_table(dp, "FGA", cols, label="all (purity subset)"))
        res.append(corr_table(dp, "FGA", cols, covs=["purity"], label="partial|purity"))
        res.append(corr_table(dp, "FGA", cols, covs=["purity", "mut"], label="partial|purity+mut"))
        res.append(corr_table(dp, "purity", cols, label="purity vs expression"))
        r, p = stats.spearmanr(dp["FGA"], dp["purity"], nan_policy="omit")
        print(f"FGA vs purity: rho={r:.2f}, p={p:.2g}")

    out = pd.concat(res)
    out.to_csv(f"{a.out}/correlations.csv", index=False)
    print(out.pivot(index="gene", columns="analysis", values="rho").round(2).to_string())

    # --- survival (optional) ---
    try:
        from lifelines import CoxPHFitter
        d = df.copy()
        d["event"] = d["os_s"].astype(str).str.startswith("1").astype(int)
        d["log2CYT"] = np.log2(d["CYT"] + 1)
        d["stage34"] = d["stage"].astype(str).str.contains("III|IV").astype(int)
        d["male"] = (d["sex"] == "Male").astype(int)
        keep = ["os_m", "event", "log2CYT", "age", "male", "stage34"]
        if "purity" in d: keep.append("purity")
        c = CoxPHFitter().fit(d[keep].dropna(), "os_m", "event")
        c.summary.to_csv(f"{a.out}/cox_cyt.csv")
        print(c.summary[["exp(coef)", "p"]].round(3).to_string())
    except ImportError:
        print("lifelines not installed; skipping Cox model")


if __name__ == "__main__":
    main()
