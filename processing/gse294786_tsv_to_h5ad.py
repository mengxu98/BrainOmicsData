#!/usr/bin/env python3

"""Convert the complete dense cell-by-gene TSV to a sparse H5AD matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit-output", required=True)
    parser.add_argument("--chunk-rows", type=int, default=256)
    parser.add_argument("--progress-every", type=int, default=10)
    parser.add_argument("--expected-cells", type=int)
    parser.add_argument("--expected-features", type=int)
    parser.add_argument("--expected-input-sha256")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def add_string_array(group: h5py.Group, name: str, values: list[str]) -> None:
    dtype = h5py.string_dtype(encoding="utf-8")
    dataset = group.create_dataset(
        name,
        data=np.asarray(values, dtype=object),
        dtype=dtype,
        maxshape=(None,),
    )
    dataset.attrs["encoding-type"] = "string-array"
    dataset.attrs["encoding-version"] = "0.2.0"


def create_dataframe_group(
    handle: h5py.File,
    name: str,
    index_values: list[str],
) -> h5py.Group:
    group = handle.create_group(name)
    group.attrs["encoding-type"] = "dataframe"
    group.attrs["encoding-version"] = "0.2.0"
    group.attrs["_index"] = "_index"
    group.attrs["column-order"] = np.asarray([], dtype=h5py.string_dtype())
    add_string_array(group, "_index", index_values)
    return group


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    if args.chunk_rows < 1:
        raise ValueError("--chunk-rows must be positive")
    if args.progress_every < 1:
        raise ValueError("--progress-every must be positive")
    counts_path = Path(args.counts).resolve()
    output_path = Path(args.output).resolve()
    audit_path = Path(args.audit_output).resolve()
    if not counts_path.is_file():
        raise FileNotFoundError(counts_path)
    if output_path.exists() and not args.overwrite:
        raise FileExistsError(output_path)
    input_sha256 = sha256_file(counts_path)
    if (
        args.expected_input_sha256
        and input_sha256 != args.expected_input_sha256
    ):
        raise ValueError("counts TSV SHA-256 differs from its download record")

    header = pd.read_csv(counts_path, sep="\t", nrows=0)
    if len(header.columns) < 2:
        raise ValueError("counts TSV must contain cell_id and gene columns")
    cell_column = str(header.columns[0])
    feature_names = [str(value) for value in header.columns[1:]]
    if len(set(feature_names)) != len(feature_names):
        raise ValueError("counts TSV contains duplicated feature names")
    if (
        args.expected_features is not None
        and len(feature_names) != args.expected_features
    ):
        raise ValueError(
            "feature count differs from the reviewed source manifest: "
            f"expected {args.expected_features}, observed {len(feature_names)}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = Path(str(output_path) + ".tmp")
    if temporary_path.exists():
        temporary_path.unlink()

    total_cells = 0
    total_nonzero = 0
    total_count_sum = 0
    cell_ids: list[str] = []
    with h5py.File(temporary_path, "w") as handle:
        handle.attrs["encoding-type"] = "anndata"
        handle.attrs["encoding-version"] = "0.1.0"

        x_group = handle.create_group("X")
        x_group.attrs["encoding-type"] = "csr_matrix"
        x_group.attrs["encoding-version"] = "0.1.0"
        data_ds = x_group.create_dataset(
            "data",
            shape=(0,),
            maxshape=(None,),
            chunks=True,
            dtype=np.int32,
        )
        indices_ds = x_group.create_dataset(
            "indices",
            shape=(0,),
            maxshape=(None,),
            chunks=True,
            dtype=np.int32,
        )
        indptr_ds = x_group.create_dataset(
            "indptr",
            shape=(1,),
            maxshape=(None,),
            chunks=True,
            dtype=np.int64,
        )
        indptr_ds[0] = 0

        reader = pd.read_csv(
            counts_path,
            sep="\t",
            chunksize=args.chunk_rows,
            low_memory=False,
        )
        for chunk_number, chunk in enumerate(reader, start=1):
            ids = chunk.pop(cell_column).astype(str).tolist()
            non_numeric = [
                name
                for name, dtype in chunk.dtypes.items()
                if not pd.api.types.is_numeric_dtype(dtype)
            ]
            if non_numeric:
                raise ValueError(
                    "counts TSV contains non-numeric values in features: "
                    + ", ".join(non_numeric[:10])
                )
            values_source = chunk.to_numpy(copy=False)
            if np.issubdtype(values_source.dtype, np.floating):
                if not np.all(np.isfinite(values_source)):
                    raise ValueError("counts TSV contains non-finite values")
                if not np.all(values_source == np.floor(values_source)):
                    raise ValueError("counts TSV contains fractional counts")
            values_64 = values_source.astype(np.int64, copy=False)
            if values_64.shape[1] != len(feature_names):
                raise RuntimeError("feature count changed within counts TSV")
            if np.any(values_64 < 0):
                raise ValueError("counts TSV contains negative values")
            if np.any(values_64 > np.iinfo(np.int32).max):
                raise ValueError("counts TSV exceeds the int32 count range")
            values = values_64.astype(np.int32, copy=False)

            matrix = sparse.csr_matrix(values)
            rows = matrix.shape[0]
            old_nonzero = total_nonzero
            new_nonzero = old_nonzero + matrix.nnz
            data_ds.resize((new_nonzero,))
            indices_ds.resize((new_nonzero,))
            data_ds[old_nonzero:new_nonzero] = matrix.data
            indices_ds[old_nonzero:new_nonzero] = matrix.indices

            old_cells = total_cells
            total_cells += rows
            indptr_ds.resize((total_cells + 1,))
            indptr_ds[old_cells + 1 : total_cells + 1] = (
                matrix.indptr[1:] + old_nonzero
            )
            total_nonzero = new_nonzero
            total_count_sum += int(values_64.sum(dtype=np.int64))
            cell_ids.extend(ids)
            if chunk_number % args.progress_every == 0:
                print(
                    "processed "
                    f"{total_cells} cells, {len(feature_names)} features, "
                    f"{total_nonzero} nonzero values",
                    flush=True,
                )

        if len(set(cell_ids)) != len(cell_ids):
            raise ValueError("counts TSV contains duplicated cell IDs")
        if (
            args.expected_cells is not None
            and total_cells != args.expected_cells
        ):
            raise ValueError(
                "cell count differs from the reviewed source manifest: "
                f"expected {args.expected_cells}, observed {total_cells}"
            )
        x_group.attrs["shape"] = np.asarray(
            [total_cells, len(feature_names)],
            dtype=np.int64,
        )
        create_dataframe_group(handle, "obs", cell_ids)
        create_dataframe_group(handle, "var", feature_names)
        handle.create_group("layers")
        handle.create_group("obsm")
        handle.create_group("obsp")
        handle.create_group("varm")
        handle.create_group("varp")
        handle.create_group("uns")

    os.replace(temporary_path, output_path)
    with h5py.File(output_path, "r") as handle:
        output_shape = tuple(int(value) for value in handle["X"].attrs["shape"])
        if output_shape != (total_cells, len(feature_names)):
            raise RuntimeError("installed H5AD matrix shape verification failed")
        if int(handle["X/indptr"][-1]) != total_nonzero:
            raise RuntimeError("installed H5AD sparse pointer verification failed")
    audit = {
        "input_file": str(counts_path),
        "input_bytes": counts_path.stat().st_size,
        "output_file": str(output_path),
        "output_bytes": output_path.stat().st_size,
        "input_cells_processed": total_cells,
        "input_features_processed": len(feature_names),
        "nonzero_values_processed": total_nonzero,
        "total_count_sum": total_count_sum,
        "input_sha256": input_sha256,
        "output_sha256": sha256_file(output_path),
        "all_rows_processed": True,
        "all_features_processed": True,
        "chunk_rows": args.chunk_rows,
    }
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    with audit_path.open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
