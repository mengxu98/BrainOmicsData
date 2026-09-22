#!/usr/bin/env python3

"""Convert the complete Allen human M1 CSV matrix to donor H5AD shards."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import multiprocessing as mp
import os
from pathlib import Path
import shutil

import anndata
import numpy as np
import pandas as pd
from scipy import sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--chunk-cells", type=int, default=256)
    parser.add_argument(
        "--workers",
        type=int,
        default=min(12, max(1, os.cpu_count() or 1)),
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_name(value: str) -> str:
    return "".join(
        character if character.isalnum() or character in "._-" else "_"
        for character in value
    )


_WORKER_COLUMNS: list[str] = []
_WORKER_DONOR_BY_CELL: pd.Series | None = None
_WORKER_CHUNK_CELLS = 0
_WORKER_DIR: Path | None = None


def initialise_worker(
    columns: list[str],
    donor_by_cell: pd.Series,
    chunk_cells: int,
    work_dir: str,
) -> None:
    global _WORKER_COLUMNS
    global _WORKER_DONOR_BY_CELL
    global _WORKER_CHUNK_CELLS
    global _WORKER_DIR
    _WORKER_COLUMNS = columns
    _WORKER_DONOR_BY_CELL = donor_by_cell
    _WORKER_CHUNK_CELLS = chunk_cells
    _WORKER_DIR = Path(work_dir)


def parse_partition(partition: tuple[int, int, int]) -> dict[str, object]:
    index, start, end = partition
    if _WORKER_DONOR_BY_CELL is None or _WORKER_DIR is None:
        raise RuntimeError("Allen CSV worker was not initialized")

    matrix_pieces: dict[str, list[sparse.csr_matrix]] = {}
    cell_pieces: dict[str, list[np.ndarray]] = {}
    row_pieces: dict[str, list[np.ndarray]] = {}
    observed_cells: set[str] = set()
    source_row = 0
    matrix_file = Path(_WORKER_COLUMNS[0])
    columns = _WORKER_COLUMNS[1:]

    with matrix_file.open("rb") as handle:
        handle.seek(start)
        while handle.tell() < end:
            lines: list[bytes] = []
            for _ in range(_WORKER_CHUNK_CELLS):
                if handle.tell() >= end:
                    break
                line = handle.readline()
                if not line:
                    break
                lines.append(line)
            if not lines:
                break
            chunk = pd.read_csv(
                io.BytesIO(b"".join(lines)),
                header=None,
                names=columns,
                low_memory=False,
            )
            if list(chunk.columns) != columns:
                raise ValueError("Allen matrix columns changed while streaming")
            cells = chunk.iloc[:, 0].astype(str).to_numpy()
            if len(set(cells)) != len(cells) or any(
                cell in observed_cells for cell in cells
            ):
                raise ValueError("Allen matrix contains duplicated cell IDs")
            observed_cells.update(cells)
            donors = _WORKER_DONOR_BY_CELL.reindex(cells)
            if donors.isna().any():
                missing = cells[donors.isna().to_numpy()][:5]
                raise ValueError(
                    f"Allen matrix cells lack metadata: {missing.tolist()}"
                )

            values = chunk.iloc[:, 1:].to_numpy(copy=False)
            if not np.issubdtype(values.dtype, np.number):
                raise ValueError("Allen matrix contains non-numeric expression values")
            if np.issubdtype(values.dtype, np.integer):
                if np.any(values < 0):
                    raise ValueError("Allen matrix contains negative values")
            else:
                if not np.isfinite(values).all() or np.any(values < 0):
                    raise ValueError(
                        "Allen matrix contains negative or non-finite values"
                    )
                if not np.equal(values, np.floor(values)).all():
                    raise ValueError(
                        "Allen expression matrix is not an integer count matrix"
                    )
            if values.size and values.max() > np.iinfo(np.uint32).max:
                raise ValueError("Allen count exceeds uint32 range")
            chunk_matrix = sparse.csr_matrix(values, dtype=np.uint32)
            chunk_rows = np.arange(
                source_row,
                source_row + len(cells),
                dtype=np.int64,
            )
            source_row += len(cells)

            donor_values = donors.astype(str).to_numpy()
            for donor in sorted(set(donor_values)):
                selected = np.flatnonzero(donor_values == donor)
                matrix_pieces.setdefault(donor, []).append(
                    chunk_matrix[selected, :]
                )
                cell_pieces.setdefault(donor, []).append(cells[selected])
                row_pieces.setdefault(donor, []).append(chunk_rows[selected])

    outputs: dict[str, dict[str, object]] = {}
    for donor in sorted(matrix_pieces):
        matrix = sparse.vstack(
            matrix_pieces[donor],
            format="csr",
            dtype=np.uint32,
        )
        cells = np.concatenate(cell_pieces[donor]).astype(str)
        local_rows = np.concatenate(row_pieces[donor])
        prefix = f"part_{index:02d}_{safe_name(donor)}"
        matrix_output = _WORKER_DIR / f"{prefix}.matrix.npz"
        metadata_output = _WORKER_DIR / f"{prefix}.metadata.npz"
        matrix_temporary = _WORKER_DIR / (
            f"{prefix}.matrix.tmp.{os.getpid()}.npz"
        )
        metadata_temporary = _WORKER_DIR / (
            f"{prefix}.metadata.tmp.{os.getpid()}.npz"
        )
        sparse.save_npz(matrix_temporary, matrix, compressed=False)
        np.savez(metadata_temporary, cells=cells, local_rows=local_rows)
        os.replace(matrix_temporary, matrix_output)
        os.replace(metadata_temporary, metadata_output)
        outputs[donor] = {
            "matrix": str(matrix_output),
            "metadata": str(metadata_output),
            "cells": int(matrix.shape[0]),
            "nonzero_values": int(matrix.nnz),
        }
    return {"index": index, "rows": source_row, "donors": outputs}


def matrix_partitions(
    matrix_file: Path,
    header_end: int,
    workers: int,
) -> list[tuple[int, int, int]]:
    file_size = matrix_file.stat().st_size
    boundaries = [header_end]
    with matrix_file.open("rb") as handle:
        for index in range(1, workers):
            target = header_end + ((file_size - header_end) * index // workers)
            handle.seek(target)
            handle.readline()
            boundaries.append(handle.tell())
    boundaries.append(file_size)
    boundaries = sorted(set(boundaries))
    return [
        (index, start, end)
        for index, (start, end) in enumerate(
            zip(boundaries[:-1], boundaries[1:]),
            start=1,
        )
        if end > start
    ]


def main() -> None:
    args = parse_args()
    if args.chunk_cells < 1:
        raise ValueError("--chunk-cells must be positive")
    if args.workers < 1:
        raise ValueError("--workers must be positive")

    matrix_file = Path(args.matrix).resolve()
    metadata_file = Path(args.metadata).resolve()
    output_dir = Path(args.output_dir).resolve()
    if not matrix_file.is_file() or not metadata_file.is_file():
        raise FileNotFoundError("Allen matrix.csv or metadata.csv is missing")
    output_dir.mkdir(parents=True, exist_ok=True)

    source_metadata = pd.read_csv(
        metadata_file,
        usecols=["sample_name", "external_donor_name_label"],
        dtype="string",
    )
    if (
        source_metadata["sample_name"].isna().any()
        or source_metadata["external_donor_name_label"].isna().any()
        or source_metadata["sample_name"].duplicated().any()
    ):
        raise ValueError("Allen source metadata has missing or duplicated cell IDs")
    donor_by_cell = source_metadata.set_index("sample_name")[
        "external_donor_name_label"
    ]
    expected_cells = set(source_metadata["sample_name"].astype(str))

    with matrix_file.open("rb") as handle:
        header_line = handle.readline()
        header_end = handle.tell()
    columns = next(csv.reader([header_line.decode("utf-8-sig").rstrip("\r\n")]))
    if len(columns) < 2 or columns[0] != "sample_name":
        raise ValueError("Allen matrix first column is not sample_name")
    genes = [str(value) for value in columns[1:]]
    if len(set(genes)) != len(genes):
        raise ValueError("Allen matrix contains duplicated feature names")

    partitions = matrix_partitions(
        matrix_file,
        header_end,
        min(args.workers, max(1, len(source_metadata))),
    )
    work_dir = output_dir / f".allen_csv_work.{os.getpid()}"
    work_dir.mkdir(parents=True, exist_ok=False)
    worker_columns = [str(matrix_file), *columns]
    context_name = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
    context = mp.get_context(context_name)
    results: list[dict[str, object]] = []
    with context.Pool(
        processes=len(partitions),
        initializer=initialise_worker,
        initargs=(
            worker_columns,
            donor_by_cell,
            args.chunk_cells,
            str(work_dir),
        ),
    ) as pool:
        for result in pool.imap_unordered(parse_partition, partitions):
            results.append(result)
            print(
                f"[{len(results)}/{len(partitions)}] complete Allen CSV partition",
                flush=True,
            )
    results.sort(key=lambda value: int(value["index"]))
    offsets: dict[int, int] = {}
    source_row = 0
    for result in results:
        offsets[int(result["index"])] = source_row
        source_row += int(result["rows"])
    if source_row != len(source_metadata):
        raise ValueError("Allen matrix row count differs from source metadata")

    feature_metadata = pd.DataFrame(index=pd.Index(genes, name="feature"))
    manifest_rows: list[dict[str, object]] = []
    observed_cells: set[str] = set()
    observed_source_rows = np.zeros(len(source_metadata), dtype=bool)
    donors = sorted(
        {
            donor
            for result in results
            for donor in result["donors"]
        }
    )
    for donor in donors:
        matrices: list[sparse.csr_matrix] = []
        cell_pieces: list[np.ndarray] = []
        row_pieces: list[np.ndarray] = []
        for result in results:
            donor_output = result["donors"].get(donor)
            if donor_output is None:
                continue
            matrices.append(sparse.load_npz(donor_output["matrix"]))
            with np.load(donor_output["metadata"], allow_pickle=False) as payload:
                cell_pieces.append(payload["cells"].astype(str))
                row_pieces.append(
                    payload["local_rows"].astype(np.int64)
                    + offsets[int(result["index"])]
                )
        matrix = sparse.vstack(matrices, format="csr", dtype=np.uint32)
        cells = np.concatenate(cell_pieces).astype(str)
        source_rows = np.concatenate(row_pieces)
        if matrix.shape != (len(cells), len(genes)) or len(set(cells)) != len(cells):
            raise RuntimeError(f"Allen donor shard is inconsistent: {donor}")
        if any(cell in observed_cells for cell in cells):
            raise ValueError("Allen donor shards contain duplicated cell IDs")
        if observed_source_rows[source_rows].any():
            raise ValueError("Allen donor shards contain duplicated source rows")
        observed_cells.update(cells)
        observed_source_rows[source_rows] = True
        obs = pd.DataFrame(
            {
                "Cells": cells,
                "Original_Source_Row": source_rows,
                "Original_Donor_ID": donor,
            },
            index=pd.Index(cells, name="Cells_index"),
        )
        adata = anndata.AnnData(X=matrix, obs=obs, var=feature_metadata)
        filename = f"{safe_name(donor)}.h5ad"
        output = output_dir / filename
        temporary = output.with_name(f"{output.name}.tmp.{os.getpid()}")
        adata.write_h5ad(temporary, compression="gzip")
        os.replace(temporary, output)
        manifest_rows.append(
            {
                "donor": donor,
                "file": filename,
                "cells": matrix.shape[0],
                "features": matrix.shape[1],
                "nonzero_values": matrix.nnz,
                "first_source_row": int(source_rows.min()),
                "last_source_row": int(source_rows.max()),
                "sha256": sha256(output),
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    if (
        manifest["cells"].sum() != len(source_metadata)
        or observed_cells != expected_cells
        or not observed_source_rows.all()
    ):
        raise RuntimeError("Allen donor shards do not cover every source cell")
    manifest_output = output_dir / "manifest.tsv"
    temporary_manifest = manifest_output.with_name(
        f"{manifest_output.name}.tmp.{os.getpid()}"
    )
    manifest.to_csv(temporary_manifest, sep="\t", index=False)
    os.replace(temporary_manifest, manifest_output)
    shutil.rmtree(work_dir)


if __name__ == "__main__":
    main()
