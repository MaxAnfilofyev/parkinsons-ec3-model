"""Shared utilities for organizing staged figure-generation runs.

This module mirrors the directory/run structure used in `code/run_ec3_all.py` so
that new supplementary scripts (e.g., for S1–S7) can save outputs in a
consistent, timestamped way.

Usage pattern
-------------

    from datetime import datetime
    from stage_utils import make_stage_run_dir, write_readme, STAGE_NAMES

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_dir = make_stage_run_dir(STAGE_NAMES["S1_S2_EarlierVariants"], timestamp)

    # Save figures / CSVs into `run_dir`...

    write_readme(
        run_dir,
        header="Stage S1–S2: Earlier Monostable / Feedback-Only Variants",
        body_lines=[
            f"Timestamp: {timestamp}",
            "",
            "Files:",
            "- ...",
            "",
            "Description:",
            "Figures corresponding to Supplementary S1–S2 (earlier model variants).",
        ],
    )

Each stage directory will have the structure:

    <stage_name>/
        runs/
            <timestamp_1>/
                ... outputs for that run ...
            <timestamp_2>/
                ...
        latest_run -> runs/<most_recent_timestamp>/   # symlink if possible

This is identical in spirit to the existing Stage_2/3/4 directories.
"""

from __future__ import annotations

import os
import shutil
from typing import Iterable


# Recommended stage names for supplementary / extended analyses.
# These are just string constants used as the `stage_name` argument to
# `make_stage_run_dir`.
STAGE_NAMES = {
    # Earlier monostable / feedback-only model variants (Supplementary S1–S2).
    "S1_S2_EarlierVariants": "Stage_S1_S2_EarlierVariants",
    # Effect of nonlinear energy amplification on nullcline geometry (S3).
    "S3_NonlinearAmplification": "Stage_S3_NonlinearAmplification",
    # EC3 robustness sweeps over k_M, L1, C (S4–S6).
    "S4_S6_RobustnessSweeps": "Stage_S4_S6_RobustnessSweeps",
    # Full continuation of equilibria vs A showing right-hand fold (S7).
    "S7_FullContinuation": "Stage_S7_FullContinuation",
    # Optional: EC2-specific supplementary figures (e.g., from `supp.py`).
    "EC2_Supplement": "Stage_EC2_Supplement",
}


def make_stage_run_dir(stage_name: str, timestamp: str) -> str:
    """Create `<stage_name>/runs/<timestamp>/` and update `latest_run`.

    This is modeled after `make_stage_run_dir` in `run_ec3_all.py`.

    Parameters
    ----------
    stage_name : str
        Name of the stage directory to create (e.g., "Stage_S3_NonlinearAmplification").
    timestamp : str
        Timestamp string used to name the run subdirectory, e.g.,
        `datetime.now().strftime("%Y-%m-%d_%H-%M-%S")`.

    Returns
    -------
    run_dir : str
        Absolute path to the created run directory.
    """
    root = os.getcwd()
    stage_dir = os.path.join(root, stage_name)
    runs_dir = os.path.join(stage_dir, "runs")
    os.makedirs(runs_dir, exist_ok=True)

    run_dir = os.path.join(runs_dir, timestamp)
    os.makedirs(run_dir, exist_ok=True)

    latest_link = os.path.join(stage_dir, "latest_run")

    # Remove any existing `latest_run` so we can recreate it as a symlink.
    if os.path.islink(latest_link) or os.path.isfile(latest_link):
        try:
            os.remove(latest_link)
        except OSError:
            pass
    elif os.path.isdir(latest_link):
        try:
            shutil.rmtree(latest_link)
        except OSError:
            pass

    # Try to create a symlink; if not allowed, write a text file instead.
    try:
        os.symlink(run_dir, latest_link, target_is_directory=True)
    except OSError:
        with open(os.path.join(stage_dir, "latest_run.txt"), "w") as f:
            f.write(run_dir)

    return run_dir


def write_readme(run_dir: str, header: str, body_lines: Iterable[str]) -> None:
    """Create a simple `README.txt` summarizing a run.

    Parameters
    ----------
    run_dir : str
        The run directory in which to place the README.
    header : str
        Title for the README.
    body_lines : Iterable[str]
        Lines to write after the header and underline.
    """
    path = os.path.join(run_dir, "README.txt")
    with open(path, "w") as f:
        f.write(header + "\n")
        f.write("=" * len(header) + "\n\n")
        for line in body_lines:
            f.write(line + "\n")
