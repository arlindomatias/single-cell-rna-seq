# Single-cell RNA-seq analyses

Scanpy-based single-cell RNA-seq analyses of public GEO datasets on chemotherapy-induced
changes in the mouse nervous system. Each project lives in its own folder with its own
workflow, conda environments and README.

| Folder | Dataset | What it does | Status |
|---|---|---|---|
| [`single_cell_cisp/`](single_cell_cisp) | [GSE216146](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216146) | QC, Harmony integration, Leiden clustering, CellTypist annotation, cell-type composition, pseudobulk DE, enrichment (chemotherapy vs control) | In progress: exploratory script in `notebooks/`, Snakemake port pending |

## Environment

```bash
conda env create -f environment.yaml
conda activate single-cell
```

`environment.yaml` provides Snakemake plus the analysis stack (scanpy, harmonypy, leidenalg,
celltypist, pydeseq2, decoupler, gseapy). Rules use the per-project envs in `<project>/envs/`.

Raw data are not versioned. They are read from the directory in the `DATA` environment variable.

## License

MIT
