#!/usr/bin/env python3
"""Train scVI on every cell exported by the R reference pipeline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import tempfile

import anndata as ad
from lightning.pytorch.callbacks import ModelCheckpoint
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
import scvi
import torch


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_lines(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return [line.rstrip("\r\n") for line in handle]


def sha256_text_lines(values: list[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "Dataset",
        "Order",
        "Source_Matrix_Units",
        "Source_Provenance_Eligible",
        "Numeric_Count_Eligible",
        "ScVI_Eligible",
        "Cells",
        "Features",
        "Nonzero_Values",
        "Data_File",
        "Data_SHA256",
        "Indices_File",
        "Indices_SHA256",
        "Indptr_File",
        "Indptr_SHA256",
        "Cells_File",
        "Cells_SHA256",
        "Batch_File",
        "Batch_SHA256",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError("scVI input manifest is empty or incomplete")
    rows.sort(key=lambda row: int(row["Order"]))
    if [int(row["Order"]) for row in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("scVI input manifest order is not contiguous")
    return rows


def checked_path(root: Path, relative: str) -> Path:
    candidate = (root / relative).resolve()
    if root.resolve() not in candidate.parents:
        raise ValueError(f"unsafe scVI input path: {relative}")
    if not candidate.is_file():
        raise FileNotFoundError(candidate)
    return candidate


def load_input(input_dir: Path) -> tuple[ad.AnnData, list[str], list[str]]:
    manifest_file = input_dir / "manifest.tsv"
    audit_file = input_dir / "input_audit.tsv"
    feature_file = input_dir / "features.txt"
    rows = read_manifest(manifest_file)
    features = read_lines(feature_file)
    if len(features) < 2 or len(features) != len(set(features)):
        raise ValueError("scVI feature list is missing or duplicated")

    matrices: list[sparse.csr_matrix] = []
    cell_ids: list[str] = []
    batches: list[str] = []
    datasets: list[str] = []
    for row in rows:
        cells = int(row["Cells"])
        n_features = int(row["Features"])
        nnz = int(row["Nonzero_Values"])
        if n_features != len(features) or cells < 1 or nnz < 0:
            raise ValueError(f"invalid dimensions for {row['Dataset']}")
        if (
            row["Source_Matrix_Units"]
            != "raw non-negative integer gene counts"
            or row["Source_Provenance_Eligible"] != "TRUE"
            or row["Numeric_Count_Eligible"] != "TRUE"
            or row["ScVI_Eligible"] != "TRUE"
        ):
            raise ValueError(
                f"{row['Dataset']} did not pass the scVI source-and-value gate"
            )

        resolved: dict[str, Path] = {}
        for prefix in ("Data", "Indices", "Indptr", "Cells", "Batch"):
            path = checked_path(input_dir, row[f"{prefix}_File"])
            if sha256_file(path) != row[f"{prefix}_SHA256"]:
                raise ValueError(
                    f"{row['Dataset']} {prefix.lower()} SHA-256 differs"
                )
            resolved[prefix] = path

        data = np.fromfile(resolved["Data"], dtype="<i4")
        indices = np.fromfile(resolved["Indices"], dtype="<i4")
        indptr = np.fromfile(resolved["Indptr"], dtype="<i4")
        shard_cells = read_lines(resolved["Cells"])
        shard_batches = read_lines(resolved["Batch"])
        if (
            data.size != nnz
            or indices.size != nnz
            or indptr.size != cells + 1
            or len(shard_cells) != cells
            or len(shard_batches) != cells
        ):
            raise ValueError(f"binary dimensions differ for {row['Dataset']}")
        if (
            np.any(data < 0)
            or np.any(indices < 0)
            or np.any(indices >= n_features)
            or indptr[0] != 0
            or indptr[-1] != nnz
            or np.any(np.diff(indptr.astype(np.int64)) < 0)
        ):
            raise ValueError(f"invalid sparse structure for {row['Dataset']}")
        if len(set(shard_cells)) != cells or any(not value for value in shard_batches):
            raise ValueError(f"missing or duplicate metadata for {row['Dataset']}")

        matrix = sparse.csc_matrix(
            (data.astype(np.float32), indices, indptr),
            shape=(n_features, cells),
        ).transpose().tocsr()
        matrices.append(matrix)
        cell_ids.extend(shard_cells)
        batches.extend(shard_batches)
        datasets.extend([row["Dataset"]] * cells)

    if len(cell_ids) != len(set(cell_ids)):
        raise ValueError("scVI input contains duplicated cell identifiers")
    matrix = sparse.vstack(matrices, format="csr")
    del matrices
    obs = pd.DataFrame(
        {
            "Integration_Batch_ID": pd.Categorical(batches),
            "Dataset": pd.Categorical(datasets),
        },
        index=pd.Index(cell_ids, name="Cells"),
    )
    var = pd.DataFrame(index=pd.Index(features, name="Gene"))
    adata = ad.AnnData(X=matrix, obs=obs, var=var)

    audit = pd.read_csv(audit_file, sep="\t")
    if len(audit) != 1:
        raise ValueError("scVI input audit must contain one row")
    expected = audit.iloc[0]
    if (
        int(expected["Datasets"]) != len(rows)
        or int(expected["Reference_Datasets"])
        != int(expected["Datasets"]) + int(expected["Excluded_Datasets"])
        or int(expected["Cells"]) != adata.n_obs
        or int(expected["Features"]) != adata.n_vars
        or int(expected["Nonzero_Values"]) != adata.X.nnz
        or expected["Cell_Order_SHA256"] != sha256_text_lines(cell_ids)
        or expected["Feature_Order_SHA256"] != sha256_text_lines(features)
    ):
        raise ValueError("scVI input audit differs from reconstructed data")
    return adata, cell_ids, features


def write_tsv(path: Path, columns: list[str], values: list[object]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerow(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=20260730)
    parser.add_argument("--latent-dimensions", type=int, default=50)
    parser.add_argument("--max-epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=2048)
    parser.add_argument("--n-hidden", type=int, default=128)
    parser.add_argument("--n-layers", type=int, default=2)
    parser.add_argument("--dropout-rate", type=float, default=0.1)
    parser.add_argument(
        "--dispersion",
        choices=("gene", "gene-batch", "gene-label", "gene-cell"),
        default="gene",
    )
    parser.add_argument(
        "--gene-likelihood",
        choices=("zinb", "nb", "poisson"),
        default="nb",
    )
    parser.add_argument("--train-size", type=float, default=0.9)
    parser.add_argument("--validation-size", type=float, default=0.1)
    parser.add_argument("--early-stopping-patience", type=int, default=15)
    parser.add_argument("--kl-warmup-steps", type=int, default=400)
    parser.add_argument("--learning-rate", type=float, default=0.001)
    parser.add_argument("--accelerator", choices=("auto", "cpu", "gpu"), default="gpu")
    parser.add_argument("--devices", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and not args.overwrite:
        raise FileExistsError(output_dir)
    if (
        args.latent_dimensions < 2
        or args.max_epochs < 1
        or args.batch_size < 1
        or args.n_hidden < 1
        or args.n_layers < 1
        or not 0.0 <= args.dropout_rate < 1.0
        or not 0.0 < args.train_size < 1.0
        or not 0.0 < args.validation_size < 1.0
        or args.train_size + args.validation_size > 1.0
        or args.early_stopping_patience < 1
        or args.kl_warmup_steps < 1
        or args.learning_rate <= 0.0
    ):
        raise ValueError("invalid scVI training parameter")
    if args.accelerator == "gpu" and not torch.cuda.is_available():
        raise RuntimeError("GPU scVI was requested but CUDA is unavailable")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    training_dir = output_dir.with_name(output_dir.name + ".training")
    if args.overwrite and training_dir.exists():
        shutil.rmtree(training_dir)
    training_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = training_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    contract_file = training_dir / "run_contract.json"
    run_contract = {
        "input_manifest_sha256": sha256_file(input_dir / "manifest.tsv"),
        "input_audit_sha256": sha256_file(input_dir / "input_audit.tsv"),
        "features_sha256": sha256_file(input_dir / "features.txt"),
        "seed": args.seed,
        "latent_dimensions": args.latent_dimensions,
        "max_epochs": args.max_epochs,
        "batch_size": args.batch_size,
        "n_hidden": args.n_hidden,
        "n_layers": args.n_layers,
        "dropout_rate": args.dropout_rate,
        "dispersion": args.dispersion,
        "gene_likelihood": args.gene_likelihood,
        "train_size": args.train_size,
        "validation_size": args.validation_size,
        "early_stopping_patience": args.early_stopping_patience,
        "kl_warmup_steps": args.kl_warmup_steps,
        "learning_rate": args.learning_rate,
        "accelerator": args.accelerator,
        "devices": args.devices,
        "scvi_tools": scvi.__version__,
        "torch": torch.__version__,
    }
    if contract_file.exists():
        with contract_file.open("r", encoding="utf-8") as handle:
            previous_contract = json.load(handle)
        if previous_contract != run_contract:
            raise ValueError(
                "existing scVI training checkpoint contract differs; "
                "use --overwrite only after auditing the changed inputs"
            )
    else:
        temporary_contract = contract_file.with_suffix(".json.tmp")
        with temporary_contract.open("w", encoding="utf-8") as handle:
            json.dump(run_contract, handle, indent=2, sort_keys=True)
        temporary_contract.replace(contract_file)

    scvi.settings.seed = args.seed
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
        torch.set_float32_matmul_precision("high")

    adata, cell_ids, features = load_input(input_dir)
    scvi.model.SCVI.setup_anndata(adata, batch_key="Integration_Batch_ID")
    model = scvi.model.SCVI(
        adata,
        n_hidden=args.n_hidden,
        n_latent=args.latent_dimensions,
        n_layers=args.n_layers,
        dropout_rate=args.dropout_rate,
        dispersion=args.dispersion,
        gene_likelihood=args.gene_likelihood,
    )
    checkpoint_callback = ModelCheckpoint(
        dirpath=checkpoint_dir,
        filename="epoch-{epoch:02d}-step-{step}",
        monitor="elbo_validation",
        mode="min",
        save_top_k=1,
        save_last=True,
        every_n_epochs=1,
    )
    last_checkpoint = checkpoint_dir / "last.ckpt"
    resume_checkpoint = str(last_checkpoint) if last_checkpoint.is_file() else None
    model.train(
        max_epochs=args.max_epochs,
        batch_size=args.batch_size,
        train_size=args.train_size,
        validation_size=args.validation_size,
        accelerator=args.accelerator,
        devices=args.devices,
        early_stopping=True,
        early_stopping_monitor="elbo_validation",
        early_stopping_mode="min",
        early_stopping_patience=args.early_stopping_patience,
        check_val_every_n_epoch=1,
        plan_kwargs={
            "lr": args.learning_rate,
            "n_epochs_kl_warmup": None,
            "n_steps_kl_warmup": args.kl_warmup_steps,
        },
        enable_progress_bar=True,
        enable_checkpointing=True,
        callbacks=[checkpoint_callback],
        default_root_dir=str(training_dir),
        ckpt_path=resume_checkpoint,
    )
    history_values = getattr(model, "history", None)
    checkpoint_state = torch.load(
        last_checkpoint,
        map_location="cpu",
        weights_only=False,
    )
    completed_steps = int(checkpoint_state.get("global_step", -1))
    if completed_steps < 0:
        raise ValueError("scVI checkpoint does not record the global step")
    final_kl_weight = min(1.0, completed_steps / args.kl_warmup_steps)
    if not np.isclose(final_kl_weight, 1.0):
        raise ValueError(
            "scVI training ended before full KL weight: "
            f"{final_kl_weight}"
        )
    latent = np.asarray(
        model.get_latent_representation(batch_size=args.batch_size),
        dtype="<f4",
        order="C",
    )
    if latent.shape != (adata.n_obs, args.latent_dimensions) or not np.isfinite(latent).all():
        raise ValueError("scVI returned an invalid latent matrix")

    staging = Path(tempfile.mkdtemp(prefix=output_dir.name + ".tmp.", dir=output_dir.parent))
    installed = False
    try:
        latent_file = staging / "scvi_latent.float32.bin"
        cells_file = staging / "cell_ids.txt"
        latent.tofile(latent_file)
        with cells_file.open("w", encoding="utf-8", newline="") as handle:
            for value in cell_ids:
                handle.write(value + "\n")
        model.save(staging / "model", overwrite=True, save_anndata=False)
        if history_values is not None and len(history_values) > 0:
            history = pd.concat(history_values, axis=1)
            history.to_csv(staging / "training_history.tsv", sep="\t", index=True)

        write_tsv(
            staging / "scvi_output_audit.tsv",
            [
                "Cells",
                "Features",
                "Latent_Dimensions",
                "Seed",
                "Max_Epochs",
                "Batch_Size",
                "Hidden_Units",
                "Hidden_Layers",
                "Dropout_Rate",
                "Dispersion",
                "Gene_Likelihood",
                "Train_Size",
                "Validation_Size",
                "Early_Stopping_Patience",
                "KL_Warmup_Steps",
                "Final_KL_Weight",
                "Learning_Rate",
                "Batch_Key",
                "Accelerator",
                "Devices",
                "Input_Manifest_SHA256",
                "Cell_Order_SHA256",
                "Feature_Order_SHA256",
                "Latent_SHA256",
                "Resumed_From_Checkpoint",
                "Training_Checkpoint",
            ],
            [
                adata.n_obs,
                adata.n_vars,
                args.latent_dimensions,
                args.seed,
                args.max_epochs,
                args.batch_size,
                args.n_hidden,
                args.n_layers,
                args.dropout_rate,
                args.dispersion,
                args.gene_likelihood,
                args.train_size,
                args.validation_size,
                args.early_stopping_patience,
                args.kl_warmup_steps,
                final_kl_weight,
                args.learning_rate,
                "Integration_Batch_ID",
                args.accelerator,
                args.devices,
                sha256_file(input_dir / "manifest.tsv"),
                sha256_text_lines(cell_ids),
                sha256_text_lines(features),
                sha256_file(latent_file),
                str(resume_checkpoint is not None),
                str(last_checkpoint),
            ],
        )
        versions = {
            "Python": platform.python_version(),
            "anndata": ad.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scvi-tools": scvi.__version__,
            "torch": torch.__version__,
            "CUDA_available": str(torch.cuda.is_available()),
            "CUDA_version": str(torch.version.cuda),
            "GPU": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "NA",
        }
        with (staging / "scvi_versions.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("Package", "Version"))
            writer.writerows(versions.items())
        with (staging / "run_parameters.json").open("w", encoding="utf-8") as handle:
            json.dump(vars(args), handle, default=str, indent=2, sort_keys=True)
        shutil.copy2(contract_file, staging / "training_run_contract.json")

        backup = None
        if output_dir.exists():
            backup = output_dir.with_name(output_dir.name + f".previous.{os.getpid()}")
            output_dir.rename(backup)
        staging.rename(output_dir)
        installed = True
        if backup is not None:
            shutil.rmtree(backup)
    finally:
        if not installed and staging.exists():
            shutil.rmtree(staging)


if __name__ == "__main__":
    main()
