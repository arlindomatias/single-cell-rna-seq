# Single-cell RNA-seq analyses

Scanpy-based single-cell RNA-seq analyses of public GEO datasets on chemotherapy-induced
changes in the mouse nervous system. Each project lives in its own folder with its own
notebooks, configuration and README.

| Folder | Dataset | What it does | Status |
|---|---|---|---|
| [`single_cell_notebook/`](single_cell_notebook) | [GSE216146](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216146) | Cisplatin vs control mouse brain: QC, Harmony integration, CellTypist annotation, composition, pseudobulk DE, GSEA incl. SenMayo, senescence scores (SenePy, SenMayo), cell-cell communication (LIANA+), pseudotime and RNA velocity of oligodendrocyte maturation | Complete (10 notebooks) |

## Environment

```bash
conda env create -f environment.yaml
conda activate single-cell
python -m ipykernel install --user --name single-cell
```

`environment.yaml` provides the analysis stack (scanpy, harmonypy, leidenalg, celltypist, pydeseq2,
decoupler, gseapy, statsmodels; senepy and liana via pip). The analyses are run as Jupyter notebooks.
A project may add its own environment in `<project>/envs/` when a tool needs incompatible versions
(e.g. `single_cell_notebook/envs/velocity.yaml` for scVelo).

Raw data are not versioned. They are read from the directory in the `DATA` environment variable.

## License

MIT
