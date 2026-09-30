"""Shared config loading: config/config.yaml is the single source of truth.

STATGEN_CONFIG env var overrides the path so alternate configs can
drive the same DAG entry points.
"""

from __future__ import annotations

import os
import subprocess

import yaml

DEFAULT_CONFIG = "config/config.yaml"
ENV_VAR = "STATGEN_CONFIG"


def config_path() -> str:
    return os.environ.get(ENV_VAR, DEFAULT_CONFIG)


def load_config() -> dict:
    with open(config_path()) as f:
        return yaml.safe_load(f)


def git_sha() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:
        return "unknown"
