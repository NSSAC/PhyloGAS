"""Configuration loading for the PhyloGAS pipeline.

A single YAML file drives every stage. This module loads it, expands
``{placeholders}`` against previously-defined values, and provides dotted-path
lookup so callers can write ``cfg.get("genetic_painter.initial_viral_load")``.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "PyYAML is required to read the PhyloGAS config.\n"
        "  pip install pyyaml   (or use the bundled environment.yml)"
    ) from exc


_PLACEHOLDER = re.compile(r"\{([a-zA-Z0-9_.]+)\}")


class ConfigError(Exception):
    """Raised when the configuration is missing or malformed."""


class Config:
    """A loaded PhyloGAS configuration."""

    def __init__(self, data: dict, path: Path | None = None):
        self.data = data
        self.path = path
        self._expand_all()

    # ---------------------------------------------------------------- loading
    @classmethod
    def load(cls, path: str | os.PathLike) -> "Config":
        p = Path(path)
        if not p.is_file():
            raise ConfigError(
                f"Config file not found: {p}\n"
                "  Copy the template to get started:\n"
                "    cp config.template.yaml config.yaml"
            )
        with p.open() as fh:
            data = yaml.safe_load(fh) or {}
        if not isinstance(data, dict):
            raise ConfigError(f"Config root must be a mapping, got {type(data).__name__}")
        return cls(data, p)

    # ------------------------------------------------------------- expansion
    def _expand_all(self) -> None:
        """Resolve ``{dotted.key}`` placeholders, repeating until stable.

        Allows a config to say::

            data_dir: "data"
            population:
              persontrait_file: "{data_dir}/{population.state}_persontrait.csv"
        """
        for _ in range(10):  # bounded: guards against circular references
            if not self._expand_pass(self.data):
                return
        raise ConfigError(
            "Could not resolve placeholders after 10 passes - check for a "
            "circular reference such as a: '{b}' / b: '{a}'"
        )

    def _expand_pass(self, node: Any, _changed: list | None = None) -> bool:
        changed = _changed if _changed is not None else []
        if isinstance(node, dict):
            for k, v in node.items():
                if isinstance(v, str):
                    new = self._expand_str(v)
                    if new != v:
                        node[k] = new
                        changed.append(True)
                else:
                    self._expand_pass(v, changed)
        elif isinstance(node, list):
            for i, v in enumerate(node):
                if isinstance(v, str):
                    new = self._expand_str(v)
                    if new != v:
                        node[i] = new
                        changed.append(True)
                else:
                    self._expand_pass(v, changed)
        return bool(changed)

    def _expand_str(self, s: str) -> str:
        def sub(m: re.Match) -> str:
            val = self.get(m.group(1), default=None)
            # Leave unresolved placeholders alone so a later pass can fill them.
            return str(val) if val is not None and not isinstance(val, (dict, list)) else m.group(0)

        return _PLACEHOLDER.sub(sub, s)

    # ---------------------------------------------------------------- access
    def get(self, dotted: str, default: Any = _PLACEHOLDER) -> Any:
        """Look up ``a.b.c``. Raises ConfigError if absent and no default given."""
        node: Any = self.data
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                if default is _PLACEHOLDER:
                    raise ConfigError(
                        f"Missing required config key: '{dotted}'"
                        + (f"  (in {self.path})" if self.path else "")
                    )
                return default
            node = node[part]
        return node

    def require_path(self, dotted: str) -> Path:
        """Look up a key and assert the file it names exists."""
        raw = self.get(dotted)
        p = Path(str(raw)).expanduser()
        if not p.exists():
            raise ConfigError(
                f"Config key '{dotted}' points at a missing path:\n  {p}\n"
                "  Run `phylogas fetch-data` if this is a Dataverse-hosted input."
            )
        return p

    def results_dir(self, stage: str) -> Path:
        """Return (and create) the output directory for a named pipeline stage."""
        root = Path(str(self.get("results_dir", default="results")))
        d = root / stage
        d.mkdir(parents=True, exist_ok=True)
        return d

    def __repr__(self) -> str:
        name = self.get("project_name", default="<unnamed>")
        return f"<Config project_name={name!r} path={self.path}>"
