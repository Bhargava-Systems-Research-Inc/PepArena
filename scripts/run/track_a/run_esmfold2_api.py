#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import time
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

import pandas as pd
from tqdm import tqdm

SUMMARY_COLUMNS = [
    "pdb_id",
    "mhc_class",
    "model",
    "num_loops",
    "num_sampling_steps",
    "total_length",
    "mean_plddt",
    "ptm",
    "iptm",
    "status",
    "error",
    "runtime_seconds",
    "pred_path",
]

DEFAULT_ENV_FILE = Path(__file__).resolve().parents[1] / ".env"

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Submit TCR-pMHC complexes to Biohub ESMFold2.")
    parser.add_argument("--input-csv", type=Path, default=Path("data/esmfold2_inputs.csv"))
    parser.add_argument("--out-dir", type=Path, default=Path("predictions/esmfold2"))
    parser.add_argument("--summary-csv", type=Path, default=Path("results/esmfold2_api_results.csv"))
    parser.add_argument("--model", default="esmfold2-fast-2026-05")
    parser.add_argument("--num-loops", type=int, default=3)
    parser.add_argument("--num-sampling-steps", type=int, default=32)
    parser.add_argument("--limit", type=int, default=None, help="Only run the first N valid rows.")
    parser.add_argument(
        "--max-new",
        type=int,
        default=None,
        help="Stop after attempting N predictions that do not already have output files.",
    )
    parser.add_argument("--force", action="store_true", help="Rerun even when output already exists.")
    parser.add_argument("--sleep", type=float, default=0.0, help="Seconds to sleep between API calls.")
    return parser.parse_args()

def import_biohub_client():
    try:
        from esm.sdk.api import FoldingConfig
        from esm.sdk.forge import SequenceStructureForgeInferenceClient
        from esm.utils.structure.input_builder import ProteinInput, StructurePredictionInput

        return SequenceStructureForgeInferenceClient, FoldingConfig, ProteinInput, StructurePredictionInput
    except ImportError as first_exc:
        try:
            from esm.sdk.api import FoldingConfig, ProteinInput, StructurePredictionInput
            from esm.sdk.forge import SequenceStructureForgeInferenceClient

            return SequenceStructureForgeInferenceClient, FoldingConfig, ProteinInput, StructurePredictionInput
        except ImportError as second_exc:
            raise ImportError(
                "Could not import Biohub ESM SDK classes. Install with: "
                'pip install "esm@git+https://github.com/Biohub/esm.git@main"'
            ) from second_exc or first_exc

def read_env_file(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    if not path.exists():
        return values

    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip('"').strip("'")
        if key:
            values[key] = value
    return values

def require_token() -> str:
    token = os.environ.get("BIOHUB_API_TOKEN")
    if not token:
        token = read_env_file(DEFAULT_ENV_FILE).get("BIOHUB_API_TOKEN")
    if not token:
        raise SystemExit(
            "Missing BIOHUB_API_TOKEN.\n"
            'Set it with:\n  export BIOHUB_API_TOKEN="your_token_here"\n'
            f"or add BIOHUB_API_TOKEN=... to {DEFAULT_ENV_FILE}"
        )
    return token

def is_valid_input_row(row: pd.Series) -> bool:
    return str(row.get("status", "")).strip().lower() == "ok"

def chains_for_class(mhc_class: str) -> list[str]:
    if mhc_class == "class_II":
        return ["A", "B", "C", "D", "E"]
    if mhc_class == "class_I":
        return ["A", "C", "D", "E"]
    raise ValueError(f"Unsupported mhc_class: {mhc_class}")

def instantiate_protein_input(ProteinInput, chain_id: str, sequence: str):
    attempts = [
        {"chain_id": chain_id, "sequence": sequence},
        {"id": chain_id, "sequence": sequence},
        {"name": chain_id, "sequence": sequence},
    ]
    errors = []
    for kwargs in attempts:
        try:
            return ProteinInput(**kwargs)
        except TypeError as exc:
            errors.append(str(exc))
    raise TypeError(f"Could not construct ProteinInput for chain {chain_id}: {' | '.join(errors)}")

def instantiate_structure_input(StructurePredictionInput, proteins: list[Any]):
    attempts = [
        {"proteins": proteins},
        {"protein_inputs": proteins},
        {"sequences": proteins},
    ]
    errors = []
    for kwargs in attempts:
        try:
            return StructurePredictionInput(**kwargs)
        except TypeError as exc:
            errors.append(str(exc))
    try:
        return StructurePredictionInput(proteins)
    except TypeError as exc:
        errors.append(str(exc))
    raise TypeError(f"Could not construct StructurePredictionInput: {' | '.join(errors)}")

def build_prediction_input(row: pd.Series, ProteinInput, StructurePredictionInput):
    proteins = []
    for chain_id in ["A", "B", "C", "D", "E"]:
        sequence = str(row.get(f"chain_{chain_id}_sequence", "")).strip()
        if not sequence or sequence.lower() == "nan":
            continue
        proteins.append(instantiate_protein_input(ProteinInput, chain_id, sequence))
    return instantiate_structure_input(StructurePredictionInput, proteins)

def output_path(args: argparse.Namespace, pdb_id: str) -> Path:
    config_dir = f"{args.model}_loops{args.num_loops}_steps{args.num_sampling_steps}"
    return args.out_dir / config_dir / pdb_id / f"{pdb_id}_model_0.cif"

def read_existing_summary(path: Path) -> Dict[str, Dict[str, object]]:
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    existing: Dict[str, Dict[str, object]] = {}
    for _, row in df.iterrows():
        key = str(row["pdb_id"])
        existing[key] = {col: row.get(col, "") for col in SUMMARY_COLUMNS}
    return existing

def write_summary(path: Path, rows_by_pdb: Dict[str, Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_COLUMNS)
        writer.writeheader()
        for pdb_id in sorted(rows_by_pdb):
            row = {col: rows_by_pdb[pdb_id].get(col, "") for col in SUMMARY_COLUMNS}
            writer.writerow(row)

def base_summary_row(row: pd.Series, args: argparse.Namespace, pred_path: Path) -> Dict[str, object]:
    return {
        "pdb_id": str(row["pdb_id"]),
        "mhc_class": str(row["mhc_class"]),
        "model": args.model,
        "num_loops": args.num_loops,
        "num_sampling_steps": args.num_sampling_steps,
        "total_length": row.get("total_length", ""),
        "mean_plddt": "",
        "ptm": "",
        "iptm": "",
        "status": "",
        "error": "",
        "runtime_seconds": "",
        "pred_path": str(pred_path),
    }

def value_from_path(obj: Any, path: str) -> Optional[Any]:
    current = obj
    for part in path.split("."):
        if current is None:
            return None
        if isinstance(current, dict):
            current = current.get(part)
        else:
            current = getattr(current, part, None)
    return current

def first_metric(result: Any, names: Iterable[str]) -> Any:
    for name in names:
        value = value_from_path(result, name)
        if value is not None:
            return scalarize(value)
    return ""

def scalarize(value: Any) -> Any:
    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            pass
    if isinstance(value, (list, tuple)) and len(value) == 1:
        return scalarize(value[0])
    return value

def mean_scalar(value: Any) -> Any:
    if value is None:
        return ""
    if hasattr(value, "mean"):
        try:
            return scalarize(value.mean())
        except Exception:
            pass
    if isinstance(value, (list, tuple)):
        flat = []
        for item in value:
            if isinstance(item, (list, tuple)):
                flat.extend(item)
            else:
                flat.append(item)
        numeric = [float(item) for item in flat if item is not None]
        if numeric:
            return sum(numeric) / len(numeric)
    return scalarize(value)

def first_mean_metric(result: Any, names: Iterable[str]) -> Any:
    for name in names:
        value = value_from_path(result, name)
        if value is not None:
            return mean_scalar(value)
    return ""

def extract_metrics(result: Any) -> Dict[str, Any]:
    return {
        "mean_plddt": first_mean_metric(
            result,
            [
                "mean_plddt",
                "plddt",
                "confidence.mean_plddt",
                "confidence.plddt",
                "metrics.mean_plddt",
                "metrics.plddt",
            ],
        ),
        "ptm": first_metric(result, ["ptm", "pTM", "confidence.ptm", "confidence.pTM", "metrics.ptm"]),
        "iptm": first_metric(result, ["iptm", "ipTM", "confidence.iptm", "confidence.ipTM", "metrics.iptm"]),
    }

def mmcif_text_from_result(result: Any) -> str:
    error_msg = getattr(result, "error_msg", None)
    if error_msg:
        raise RuntimeError(str(error_msg))

    candidates = [
        "complex",
        "structure",
        "output",
        "prediction",
    ]
    objects = [result]
    for name in candidates:
        value = getattr(result, name, None)
        if value is not None:
            objects.append(value)
    if isinstance(result, dict):
        objects.extend(value for value in result.values() if value is not None)

    for obj in objects:
        for method_name in ("to_mmcif", "to_cif", "to_mmcif_string", "to_cif_string"):
            method = getattr(obj, method_name, None)
            if callable(method):
                text = method()
                if isinstance(text, bytes):
                    return text.decode("utf-8")
                if isinstance(text, str):
                    return text
        for attr_name in ("mmcif", "cif", "mmcif_string", "cif_string"):
            text = getattr(obj, attr_name, None)
            if isinstance(text, bytes):
                return text.decode("utf-8")
            if isinstance(text, str):
                return text

    available = sorted(
        set(dir(result))
        | {f"complex.{name}" for name in dir(getattr(result, "complex", object())) if not name.startswith("_")}
    )
    raise AttributeError(
        "Could not find an mmCIF serialization method on ESMFold2 result. "
        f"Inspected attributes include: {', '.join(available[:80])}"
    )

def save_prediction(result: Any, pred_path: Path) -> None:
    pred_path.parent.mkdir(parents=True, exist_ok=True)
    pred_path.write_text(mmcif_text_from_result(result))

def main() -> None:
    args = parse_args()
    token = require_token()
    SequenceStructureForgeInferenceClient, FoldingConfig, ProteinInput, StructurePredictionInput = import_biohub_client()

    if not args.input_csv.exists():
        raise FileNotFoundError(f"Input CSV does not exist: {args.input_csv}")

    input_df = pd.read_csv(args.input_csv)
    valid_df = input_df[input_df.apply(is_valid_input_row, axis=1)].copy()
    if args.limit is not None:
        valid_df = valid_df.head(args.limit)

    client = SequenceStructureForgeInferenceClient(
        model=args.model,
        url="https://biohub.ai",
        token=token,
    )
    config = FoldingConfig(
        num_loops=args.num_loops,
        num_sampling_steps=args.num_sampling_steps,
    )

    summary_rows = read_existing_summary(args.summary_csv)
    new_attempts = 0
    for _, row in tqdm(list(valid_df.iterrows()), total=len(valid_df), desc="Running ESMFold2"):
        pdb_id = str(row["pdb_id"])
        pred_path = output_path(args, pdb_id)
        summary = base_summary_row(row, args, pred_path)

        if pred_path.exists() and not args.force:
            previous = summary_rows.get(pdb_id, {})
            summary.update({key: previous.get(key, summary[key]) for key in SUMMARY_COLUMNS})
            summary["status"] = "skipped_existing"
            summary["error"] = ""
            summary["pred_path"] = str(pred_path)
            summary_rows[pdb_id] = summary
            write_summary(args.summary_csv, summary_rows)
            continue

        if args.max_new is not None and new_attempts >= args.max_new:
            print(f"Reached --max-new={args.max_new}; stopping this resumable pass.")
            break

        new_attempts += 1
        started = time.time()
        stop_after_error = False
        try:
            inp = build_prediction_input(row, ProteinInput, StructurePredictionInput)
            result = client.fold_all_atom(inp, config=config)
            save_prediction(result, pred_path)
            metrics = extract_metrics(result)
            summary.update(metrics)
            summary["status"] = "ok"
            summary["error"] = ""
        except Exception as exc:
            summary["status"] = "failed"
            summary["error"] = f"{type(exc).__name__}: {exc}"
            stop_after_error = "daily credit limit" in summary["error"].lower()
        finally:
            summary["runtime_seconds"] = round(time.time() - started, 3)
            summary_rows[pdb_id] = summary
            write_summary(args.summary_csv, summary_rows)
            if args.sleep > 0:
                time.sleep(args.sleep)

        if stop_after_error:
            print("Biohub daily credit limit reached; stopping this resumable pass.")
            break

    print(f"Wrote summary to {args.summary_csv}")

if __name__ == "__main__":
    main()
