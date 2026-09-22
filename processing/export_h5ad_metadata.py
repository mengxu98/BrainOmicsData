#!/usr/bin/env python3

"""Export complete H5AD observation and feature metadata without loading X."""

from __future__ import annotations

import argparse
import json
import os
import platform
from pathlib import Path

import anndata
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--obs-output", required=True)
    parser.add_argument("--var-output", required=True)
    parser.add_argument("--audit-output", required=True)
    parser.add_argument(
        "--matrix-group",
        required=True,
        choices=("X", "raw/X"),
    )
    return parser.parse_args()


def write_table(table: pd.DataFrame, index_name: str, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    exported = table.copy()
    exported.insert(0, index_name, exported.index.astype(str))
    exported.to_csv(
        output,
        sep="\t",
        index=False,
        compression="gzip",
    )


def main() -> None:
    args = parse_args()
    input_path = Path(args.input).resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    adata = anndata.read_h5ad(input_path, backed="r")
    try:
        obs = adata.obs.copy()
        if args.matrix_group == "raw/X":
            if adata.raw is None:
                raise ValueError("raw/X requested but the H5AD has no raw group")
            var = adata.raw.var.copy()
            var_names = adata.raw.var_names.astype(str)
        else:
            var = adata.var.copy()
            var_names = adata.var_names.astype(str)

        obs.index = adata.obs_names.astype(str)
        var.index = var_names

        obs_output = Path(args.obs_output).resolve()
        var_output = Path(args.var_output).resolve()
        audit_output = Path(args.audit_output).resolve()
        write_table(obs, "Original_Cell_ID", obs_output)
        write_table(var, "Original_Feature_ID", var_output)

        if len(obs) != adata.n_obs:
            raise RuntimeError("observation metadata export lost cells")
        if len(var) != len(var_names):
            raise RuntimeError("feature metadata export lost features")

        audit = {
            "input_file": str(input_path),
            "input_bytes": input_path.stat().st_size,
            "matrix_group": args.matrix_group,
            "input_cells": int(adata.n_obs),
            "input_features": int(len(var_names)),
            "exported_cells": int(len(obs)),
            "exported_features": int(len(var)),
            "all_cells_exported": len(obs) == adata.n_obs,
            "all_features_exported": len(var) == len(var_names),
            "anndata_version": anndata.__version__,
            "pandas_version": pd.__version__,
            "python_version": platform.python_version(),
            "pid": os.getpid(),
        }
        audit_output.parent.mkdir(parents=True, exist_ok=True)
        with audit_output.open("w", encoding="utf-8") as handle:
            json.dump(audit, handle, indent=2, sort_keys=True)
            handle.write("\n")
    finally:
        adata.file.close()


if __name__ == "__main__":
    main()
