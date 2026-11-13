"""Shared helpers for the single_cell_notebook notebooks."""
from __future__ import annotations

import gzip
import os
import re
import shutil
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def load_config(path: str | Path | None = None) -> dict:
    """Read config.yaml, expand env vars and resolve relative paths from the project root."""
    path = Path(path) if path else PROJECT_ROOT / "config.yaml"
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    for key, value in cfg["paths"].items():
        if key in ("counts_mtx", "genes_tsv", "cell_metadata_tsv", "author_h5ad_gz"):
            continue                                # file names inside data_dir
        expanded = os.path.expandvars(str(value))
        if "$" in expanded:
            raise ValueError(f"paths.{key}: environment variable not set ({value})")
        p = Path(expanded)
        cfg["paths"][key] = p if p.is_absolute() else PROJECT_ROOT / p
    for key in ("results_dir", "plots_dir", "tables_dir", "figures_dir"):
        cfg["paths"][key].mkdir(parents=True, exist_ok=True)
    return cfg


def setup_plotting() -> None:
    sc.settings.verbosity = 1
    sc.set_figure_params(dpi=80, dpi_save=300, facecolor="white", frameon=False)


def save_fig(fig: plt.Figure, cfg: dict, name: str) -> Path:
    """Save a figure as PNG in results/plots/ and return the path."""
    out = cfg["paths"]["plots_dir"] / f"{name}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    return out


def read_counts(cfg: dict, samples=None, chunksize: int = 5_000_000) -> ad.AnnData:
    """Build an AnnData of raw UMI counts from the GEO counts files (mtx + genes + cell metadata).

    The MatrixMarket file is cells x genes with 1-based integer coordinates. It is streamed in chunks
    and only the cells of `samples` (biosample IDs; None = all) are kept, so memory stays low.
    obs keeps only `qc.metadata_columns` from the author metadata; QC metrics are recomputed later.
    """
    paths, qc = cfg["paths"], cfg["qc"]
    data_dir = paths["data_dir"]
    genes = pd.read_csv(data_dir / paths["genes_tsv"], sep="\t", index_col=0)
    meta = pd.read_csv(data_dir / paths["cell_metadata_tsv"], sep="\t", index_col=0)

    with gzip.open(data_dir / paths["counts_mtx"], "rt") as fh:
        header = fh.readline()
        assert header.startswith("%%MatrixMarket matrix coordinate"), header
        assert "integer" in header, f"counts are not integer: {header.strip()}"
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        n_rows, n_cols, nnz = map(int, line.split())
        assert (n_rows, n_cols) == (len(meta), len(genes)), (
            f"mtx is {n_rows} x {n_cols}, metadata/genes have {len(meta)} cells / {len(genes)} genes")

        keep = np.ones(n_rows, bool) if samples is None else meta["biosample"].isin(list(samples)).to_numpy()
        new_row = np.full(n_rows, -1, dtype=np.int64)
        new_row[keep] = np.arange(keep.sum())

        rows, cols, vals, n_read = [], [], [], 0
        reader = pd.read_csv(fh, sep=" ", header=None, names=["r", "c", "v"],
                             dtype={"r": np.int32, "c": np.int32, "v": np.int32}, chunksize=chunksize)
        for chunk in reader:
            n_read += len(chunk)
            r = new_row[chunk["r"].to_numpy() - 1]
            m = r >= 0
            rows.append(r[m].astype(np.int32))
            cols.append(chunk["c"].to_numpy()[m] - 1)
            vals.append(chunk["v"].to_numpy()[m])
    assert n_read == nnz, f"read {n_read} entries, header says {nnz}"
    v = np.concatenate(vals)
    assert (v >= 0).all(), "negative counts"
    X = sp.csr_matrix((v.astype(np.float32), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(int(keep.sum()), n_cols))
    X.sum_duplicates()

    obs = meta.loc[keep, list(qc["metadata_columns"])].copy()
    for col in obs.select_dtypes("object"):
        obs[col] = obs[col].astype("category")
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes.index.astype(str)))
    if not adata.var_names.is_unique:
        adata.var_names_make_unique()
    return adata


def add_qc_metrics(adata: ad.AnnData, qc: dict) -> None:
    """Flag mt/ribo/hb genes and compute per-cell QC metrics in place."""
    names = adata.var_names
    adata.var["mt"] = names.str.startswith(qc["mt_prefix"])
    adata.var["ribo"] = names.str.startswith(tuple(qc["ribo_prefixes"]))
    adata.var["hb"] = names.str.contains(qc["hb_regex"], regex=True)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=False, inplace=True
    )


def filter_cells(adata: ad.AnnData, qc: dict) -> tuple[ad.AnnData, dict]:
    """Apply per-cell QC filters to one sample. Returns the filtered copy and the cells removed per step."""
    removed = {}
    keep = np.ones(adata.n_obs, dtype=bool)

    def step(name, mask):
        nonlocal keep
        removed[name] = int((keep & ~mask).sum())
        keep &= mask

    step("min_genes", adata.obs["n_genes_by_counts"].to_numpy() >= qc["min_genes"])
    upper = np.quantile(adata.obs["n_genes_by_counts"], qc["max_genes_quantile"])
    step("max_genes_quantile", adata.obs["n_genes_by_counts"].to_numpy() < upper)
    pct_mt = adata.obs["pct_counts_mt"].to_numpy()
    mt_cut = mad_upper(pct_mt, qc["mt_nmads"])
    if qc.get("max_pct_mt") is not None:
        mt_cut = min(mt_cut, qc["max_pct_mt"])
    step("mt_outlier", pct_mt <= mt_cut)
    step("max_pct_ribo", adata.obs["pct_counts_ribo"].to_numpy() < qc["max_pct_ribo"])
    step("doublet", ~adata.obs[qc["doublet_column"]].astype(bool).to_numpy())
    removed["max_genes_cutoff"] = round(float(upper), 1)
    removed["mt_cutoff"] = round(float(mt_cut), 2)
    return adata[keep].copy(), removed


def mad_upper(values: np.ndarray, nmads: float) -> float:
    """Upper outlier bound: median + nmads * MAD (unscaled median absolute deviation)."""
    med = np.median(values)
    return float(med + nmads * np.median(np.abs(values - med)))


def fraction_table(obs: pd.DataFrame, groupby: str, by: str) -> pd.DataFrame:
    """Fraction of cells of each `by` category within each `groupby` group (rows sum to 1)."""
    counts = pd.crosstab(obs[groupby], obs[by])
    return counts.div(counts.sum(axis=1), axis=0)


def plot_fractions(obs: pd.DataFrame, groupby: str, by: str, ax=None, colors=None):
    """Stacked bar chart of `by` fractions per `groupby` group."""
    frac = fraction_table(obs, groupby, by)
    ax = frac.plot(kind="bar", stacked=True, width=0.85, ax=ax, color=colors, legend=False)
    ax.set_ylabel("fraction of cells")
    ax.set_xlabel(groupby)
    ax.set_ylim(0, 1)
    ax.legend(title=by, bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
    ax.grid(False)
    return ax


def run_harmony(adata: ad.AnnData, key: str, basis: str = "X_pca",
                adjusted_basis: str = "X_pca_harmony", random_state: int = 0, **kwargs) -> None:
    """Run harmonypy on adata.obsm[basis] and store a C-ordered float32 (cells x PCs) result.

    Calls harmonypy directly: scanpy's harmony_integrate transposes Z_corr, which is wrong for
    harmonypy >= 2 (it already returns cells x PCs), and the Fortran-ordered result breaks UMAP.
    """
    import harmonypy

    ho = harmonypy.run_harmony(adata.obsm[basis], adata.obs, key,
                               random_state=random_state, verbose=False, **kwargs)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != adata.n_obs:
        Z = Z.T
    assert Z.shape == adata.obsm[basis].shape, Z.shape
    adata.obsm[adjusted_basis] = np.ascontiguousarray(Z, dtype=np.float32)


def broad_class(labels: pd.Series, rules: list) -> pd.Series:
    """Map fine labels (e.g. Allen subclasses from CellTypist) to broad classes with ordered [regex, class] rules."""
    compiled = [(re.compile(rx), cls) for rx, cls in rules]

    def one(label: str) -> str:
        for rx, cls in compiled:
            if rx.search(label):
                return cls
        return "Unmapped"

    mapping = {lab: one(str(lab)) for lab in pd.unique(labels)}
    return labels.map(mapping).astype("category")


def majority(obs: pd.DataFrame, groupby: str, col: str) -> pd.DataFrame:
    """Most frequent `col` value per `groupby` group, with its fraction."""
    grp = obs.groupby(groupby, observed=True)[col]
    return pd.DataFrame({
        f"{col}_majority": grp.agg(lambda s: s.value_counts().index[0]),
        f"{col}_majority_frac": grp.agg(lambda s: s.value_counts(normalize=True).iloc[0]),
    })


def bh_fdr(pvals) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN kept as NaN)."""
    from statsmodels.stats.multitest import multipletests

    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def ols_by_column(values: pd.DataFrame, meta: pd.DataFrame, formula: str = "~ batch + condition",
                  term: str = "condition") -> pd.DataFrame:
    """Fit `y ~ formula` (OLS) for every column of `values` (samples x features; NaN samples dropped per column).
    `meta` holds per-sample covariates (index = samples) with `term` categorical, first category = reference.
    Returns per feature: group means, coefficient of `term`, p-value, BH FDR, number of samples."""
    import statsmodels.formula.api as smf

    ref, alt = list(meta[term].cat.categories)[:2]
    coef_name = f"{term}[T.{alt}]"
    rows = []
    for col in values.columns:
        df = meta.join(values[col].rename("y"), how="inner").dropna(subset=["y"])
        ok = df[term].nunique() == 2 and len(df) > df.drop(columns="y").nunique().clip(upper=2).sum()
        fit = smf.ols(f"y {formula}", data=df).fit() if ok else None
        rows.append({"feature": col,
                     f"mean_{ref}": df.loc[df[term] == ref, "y"].mean(),
                     f"mean_{alt}": df.loc[df[term] == alt, "y"].mean(),
                     "coef": fit.params.get(coef_name, np.nan) if fit is not None else np.nan,
                     "pval": fit.pvalues.get(coef_name, np.nan) if fit is not None else np.nan,
                     "n_samples": len(df)})
    res = pd.DataFrame(rows).set_index("feature")
    res["fdr"] = bh_fdr(res["pval"])
    return res


def composition_test(counts: pd.DataFrame, meta: pd.DataFrame, formula: str = "~ batch + condition",
                     term: str = "condition") -> pd.DataFrame:
    """Propeller-style composition test.

    counts: samples x cell types (cell numbers); meta: per-sample covariates (index = samples), with `term`
    categorical and its first category as reference. Proportions are arcsin(sqrt)-transformed and an OLS model
    is fit per cell type (no variance moderation, unlike limma; with n = 3 vs 3 treat p-values as exploratory).
    """
    props = counts.div(counts.sum(axis=1), axis=0)
    ref, alt = list(meta[term].cat.categories)[:2]
    fit = ols_by_column(np.arcsin(np.sqrt(props)), meta, formula, term)
    mean_ref = props.loc[meta[term] == ref].mean()
    mean_alt = props.loc[meta[term] == alt].mean()
    res = pd.DataFrame({
        f"mean_prop_{ref}": mean_ref,
        f"mean_prop_{alt}": mean_alt,
        "log2_ratio": np.log2((mean_alt + 1e-6) / (mean_ref + 1e-6)),
        "coef": fit["coef"],
        "pval": fit["pval"],
        "fdr": fit["fdr"],
    })
    res.index.name = "cell_type"
    return res.sort_values("pval")


def slug(name: str) -> str:
    """File-name friendly version of a label ('Mural/Fibroblast' -> 'Mural_Fibroblast')."""
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")


def pseudobulk(adata: ad.AnnData, sample_key: str, group_key: str, layer: str = "counts",
               obs_cols: tuple = ()) -> ad.AnnData:
    """Sum raw counts per (sample, group). Returns an AnnData with one row per pair present in the data,
    obs = sample_key, group_key, n_cells and any per-sample `obs_cols` (must be constant within a sample)."""
    import scipy.sparse as sp

    keys = adata.obs[[sample_key, group_key]].astype(str)
    pair = keys[sample_key] + "|" + keys[group_key]
    codes, uniques = pd.factorize(pair, sort=True)
    ind = sp.csr_matrix((np.ones(adata.n_obs, dtype=np.float32), (codes, np.arange(adata.n_obs))),
                        shape=(len(uniques), adata.n_obs))
    X = adata.layers[layer] if layer else adata.X
    summed = ind @ X
    obs = pd.DataFrame([u.split("|", 1) for u in uniques], columns=[sample_key, group_key], index=uniques)
    obs["n_cells"] = np.bincount(codes)
    per_sample = adata.obs.groupby(sample_key, observed=True)[list(obs_cols)].first() if obs_cols else None
    for col in obs_cols:
        obs[col] = obs[sample_key].map(per_sample[col].astype(str))
    return ad.AnnData(X=sp.csr_matrix(summed, dtype=np.float32), obs=obs, var=pd.DataFrame(index=adata.var_names))


def plot_volcano(res: pd.DataFrame, ax, lfc: float, padj: float, n_label: int = 10, title: str = "",
                 flag_col: str | None = "flag") -> None:
    """Volcano plot of a PyDESeq2 results table (log2FoldChange, padj). Flagged genes (ambient/stress/technical) are drawn
    as open grey circles and not labelled; the top `n_label` unflagged significant genes are labelled."""
    from adjustText import adjust_text

    d = res.dropna(subset=["padj"]).copy()
    d["y"] = -np.log10(d["padj"].clip(lower=1e-300))
    flagged = d[flag_col].astype(str).ne("") if flag_col and flag_col in d else pd.Series(False, index=d.index)
    sig = (d["padj"] < padj) & (d["log2FoldChange"].abs() >= lfc)
    up, down = sig & (d["log2FoldChange"] > 0) & ~flagged, sig & (d["log2FoldChange"] < 0) & ~flagged
    ax.scatter(d.loc[~sig & ~flagged, "log2FoldChange"], d.loc[~sig & ~flagged, "y"], s=4, c="0.75", rasterized=True)
    ax.scatter(d.loc[up, "log2FoldChange"], d.loc[up, "y"], s=8, c="#c0392b", label=f"up ({up.sum()})")
    ax.scatter(d.loc[down, "log2FoldChange"], d.loc[down, "y"], s=8, c="#2c7bb6", label=f"down ({down.sum()})")
    ax.scatter(d.loc[flagged, "log2FoldChange"], d.loc[flagged, "y"], s=10, facecolors="none", edgecolors="0.4",
               linewidths=0.6, label=f"flagged ({(flagged & sig).sum()} sig.)")
    ax.axhline(-np.log10(padj), ls="--", lw=0.7, c="0.4")
    for x in (-lfc, lfc):
        ax.axvline(x, ls="--", lw=0.7, c="0.4")
    top = d[(up | down)].nsmallest(n_label, "padj")
    texts = [ax.text(r.log2FoldChange, r.y, g, fontsize=7) for g, r in top.iterrows()]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="0.5", lw=0.5))
    ax.set_xlabel("log2 fold change (chemotherapy / control)")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title(title)
    ax.grid(False)
    ax.legend(fontsize=7, frameon=False, loc="upper left")


def read_gmt(path: str | Path) -> dict:
    """Read a GMT file into {set_name: [genes]}."""
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def write_gmt(sets: dict, path: str | Path) -> None:
    with open(path, "w") as fh:
        for name, genes in sets.items():
            fh.write("\t".join([name, ""] + list(genes)) + "\n")


def load_gene_sets(cfg: dict) -> dict:
    """Enrichr libraries (downloaded once, cached as GMT in results/gene_sets/) plus custom GMT files.
    Returns {library: {set: [GENES]}} with upper-case symbols, so mouse genes match Enrichr's upper-case sets."""
    import gseapy as gp

    e = cfg["enrichment"]
    cache = cfg["paths"]["results_dir"] / "gene_sets"
    cache.mkdir(exist_ok=True)
    libs = {}
    for name in e["libraries"]:
        f = cache / f"{name}.gmt"
        if not f.exists():
            write_gmt(gp.get_library(name, organism=e["organism"]), f)
        libs[name] = read_gmt(f)
    for rel in e.get("custom_gmt", []):
        f = PROJECT_ROOT / rel
        libs[f.stem] = read_gmt(f)
    return {lib: {k: sorted({g.upper() for g in v}) for k, v in sets.items()} for lib, sets in libs.items()}


def senepy_signatures(cfg: dict, cell_types) -> tuple[dict, dict]:
    """SenePy signatures for this dataset: one 'universal' signature (all mouse cell-type signatures merged with
    SenePy's permutation threshold) and, per cell type in `senescence.matched`, a merged matched signature.
    Genes matching `senescence.exclude_gene_regex` or `annotation.ambient_genes` are removed.
    Returns ({name: [(gene, weight)]}, {cell_type: matched_name})."""
    import contextlib
    import io

    import senepy

    p = cfg["senescence"]
    drop_rx = re.compile(p["exclude_gene_regex"])
    drop = set(cfg["annotation"]["ambient_genes"])
    hubs = senepy.load_hubs(species=p["senepy_species"], sig_type="cell_type")
    meta = hubs.metadata
    with contextlib.redirect_stdout(io.StringIO()):
        hubs.merge_hubs(meta, new_name="universal", calculate_thresh=True)
    sigs = {"universal": hubs.hubs["universal"]}
    matched = {}
    for ct, rx in p["matched"].items():
        if ct not in set(cell_types):
            continue
        sel = meta[meta["cell"].str.contains(rx, regex=True)]
        name = f"matched_{slug(ct)}"
        hubs.merge_hubs(sel, new_name=name, overlap_threshold=p["matched_overlap"] if len(sel) > 1 else 1)
        sigs[name] = hubs.hubs[name]
        matched[ct] = name
    sigs = {k: [(g, w) for g, w in v if not drop_rx.search(g) and g not in drop] for k, v in sigs.items()}
    return sigs, matched


def decompress_author_h5ad(cfg: dict) -> Path:
    """Decompress the authors' h5ad.gz from $DATA into results/ once (only used for spliced/unspliced, 09/10)."""
    src = cfg["paths"]["data_dir"] / cfg["paths"]["author_h5ad_gz"]
    dst = cfg["paths"]["results_dir"] / src.name.removesuffix(".gz")
    if not dst.exists():
        tmp = dst.with_name(dst.name + ".part")
        with gzip.open(src, "rb") as fin, open(tmp, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=16 * 1024 * 1024)
        tmp.rename(dst)
    return dst


def read_h5ad_layers(path: str | Path, layers: list[str], cells: pd.Index, genes: pd.Index) -> dict:
    """Read selected CSR layers of an .h5ad for a subset of cells and genes, without loading the whole object
    (anndata's backed mode loads layers fully, and the authors' file also has large dense layers).
    Returns {layer: int32 csr_matrix (cells x genes)} in the order of `cells` and `genes`; all must exist in the file."""
    import h5py
    import scipy.sparse as sp
    from anndata.io import read_elem

    with h5py.File(path, "r") as f:
        obs_names = pd.Index(read_elem(f["obs"]).index.astype(str))
        var_names = pd.Index(read_elem(f["var"]).index.astype(str))
        rows, cols = obs_names.get_indexer(cells), var_names.get_indexer(genes)
        if (rows < 0).any() or (cols < 0).any():
            raise KeyError(f"{(rows < 0).sum()} cells / {(cols < 0).sum()} genes not found in {path}")
        out = {}
        for name in layers:
            g = f["layers"][name]
            assert g.attrs["encoding-type"] == "csr_matrix", g.attrs["encoding-type"]
            m = sp.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g.attrs["shape"]))
            out[name] = m[rows][:, cols].astype(np.int32)   # keep integers: scVelo treats float layers as normalized
            del m
    return out
