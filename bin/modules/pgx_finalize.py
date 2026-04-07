#!/usr/bin/env python3
"""
Build pgx_meta.json and pgx_result.json after PharmCAT (or on failure).
Contract: gx-exome PGx add-on for service-daemon consumption.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from datetime import datetime, timezone
from typing import Any, Dict, Optional


def _read_json(path: str) -> Optional[Any]:
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except OSError:
        return None
    except json.JSONDecodeError:
        return None


def _tail(path: str, max_bytes: int = 16_384) -> str:
    try:
        with open(path, "rb") as f:
            data = f.read()
        if len(data) > max_bytes:
            data = data[-max_bytes:]
        return data.decode("utf-8", errors="replace")
    except OSError:
        return ""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--reference-genome", default="GRCh38")
    ap.add_argument("--image", default="pgkb/pharmcat:2.15.5")
    ap.add_argument("--tool-version", default="unknown")
    ap.add_argument("--exit-code", type=int, default=1)
    ap.add_argument("--stderr", default="pharmcat.stderr.log")
    ap.add_argument("--work-dir", default=".")
    args = ap.parse_args()
    wd = args.work_dir
    rc = int(args.exit_code)

    stderr_tail = _tail(os.path.join(wd, args.stderr)) if args.stderr else ""

    # PharmCAT emits files named <base>.match.json etc.; we use -bf sample_pgx
    base_glob = os.path.join(wd, f"{args.sample_id}_pgx*")
    paths = sorted(glob.glob(base_glob))
    match_path = os.path.join(wd, f"{args.sample_id}_pgx.match.json")
    phen_path = os.path.join(wd, f"{args.sample_id}_pgx.phenotype.json")
    rep_json = os.path.join(wd, f"{args.sample_id}_pgx.report.json")
    rep_html = os.path.join(wd, f"{args.sample_id}_pgx.report.html")
    if not os.path.isfile(rep_json):
        for p in sorted(glob.glob(os.path.join(wd, f"{args.sample_id}_pgx*.json"))):
            if p.endswith(".report.json"):
                rep_json = p
                break

    match_data = _read_json(match_path)
    phen_data = _read_json(phen_path)
    rep_data = _read_json(rep_json) if os.path.isfile(rep_json) else None

    resource_versions: Dict[str, Any] = {
        "pharmcat_docker_image": args.image,
        "pharmcat_cli_version": args.tool_version,
        "cpic_data_note": "CPIC allele definitions and recommendations bundled with PharmCAT (see PharmCAT release notes / phenotype JSON).",
        "pharmgkb_note": "PharmGKB-backed data via PharmCAT release.",
    }

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    meta: Dict[str, Any] = {
        "schema": "gx_exome_pgx_meta_v1",
        "tool": "PharmCAT",
        "tool_version": args.tool_version,
        "container_image": args.image,
        "reference_genome": args.reference_genome,
        "resource_versions": resource_versions,
        "sample_id": args.sample_id,
        "timestamp": ts,
        "exit_status": "success" if rc == 0 else "error",
    }
    if rc != 0:
        meta["message"] = (stderr_tail or "PharmCAT failed; see pharmcat.stderr.log in pgx/").strip()

    out_result: Dict[str, Any] = {
        "schema": "gx_exome_pgx_result_v1",
        "sample_id": args.sample_id,
        "reference_genome": args.reference_genome,
        "tool": "PharmCAT",
        "artifacts": {
            "named_allele_matcher_json": os.path.basename(match_path) if match_data is not None else None,
            "phenotype_json": os.path.basename(phen_path) if phen_data is not None else None,
            "reporter_json": os.path.basename(rep_json) if rep_data is not None else None,
            "reporter_html": os.path.basename(rep_html) if os.path.isfile(rep_html) else None,
            "stderr_log": os.path.basename(args.stderr) if args.stderr and os.path.isfile(os.path.join(wd, args.stderr)) else None,
        },
    }

    if rc == 0 and (match_data is not None or phen_data is not None or rep_data is not None):
        out_result["named_allele_matcher"] = match_data
        out_result["phenotype"] = phen_data
        if rep_data is not None:
            out_result["reporter"] = rep_data
    elif rc == 0:
        out_result["warning"] = "PharmCAT exited 0 but expected JSON outputs were not found."
        meta["exit_status"] = "error"
        meta["message"] = out_result["warning"]
    else:
        out_result["error"] = {
            "exit_code": rc,
            "stderr_tail": stderr_tail[-4000:] if stderr_tail else None,
            "artifact_glob_matches": [os.path.basename(p) for p in paths],
        }

    with open(os.path.join(wd, "pgx_meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
        f.write("\n")

    with open(os.path.join(wd, "pgx_result.json"), "w", encoding="utf-8") as f:
        json.dump(out_result, f, indent=2, ensure_ascii=False)
        f.write("\n")


if __name__ == "__main__":
    main()
