#!/usr/bin/env python3
"""Standalone CLI entry point for exporting MicrobeTrace bundles from PhyloGAS runs.

Can be invoked either via the unified CLI:
    phylogas export-microbetrace --config configs/va.yaml --recipe 1S__surs

or directly via this script without installing the package:
    python3 scripts/utils/export_microbetrace.py --config configs/va.yaml --recipe 1S__surs
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from phylogas.microbetrace import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
