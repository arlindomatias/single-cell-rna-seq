"""Copy the key plots from results/plots/ to the tracked figures/ folder, downscaled for the README.

Usage (from single_cell_notebook/):  python scripts/export_figures.py
The list of figures and the maximum width are in config.yaml (`readme_figures`).
"""
from pathlib import Path

from PIL import Image

import sc_utils as scu


def main() -> None:
    cfg = scu.load_config()
    p = cfg["readme_figures"]
    src_dir, dst_dir = cfg["paths"]["plots_dir"], cfg["paths"]["figures_dir"]
    for name in p["files"]:
        src = src_dir / f"{name}.png"
        if not src.exists():
            raise FileNotFoundError(f"{src} (run the notebook that makes it first)")
        img = Image.open(src)
        if img.width > p["max_width"]:
            img = img.resize((p["max_width"], round(img.height * p["max_width"] / img.width)), Image.LANCZOS)
        dst = dst_dir / f"{name}.png"
        img.save(dst, optimize=True)
        print(f"{dst.relative_to(scu.PROJECT_ROOT)}  {img.width}x{img.height}  {dst.stat().st_size / 1e3:.0f} kB")


if __name__ == "__main__":
    main()
