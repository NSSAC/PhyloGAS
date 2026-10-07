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


_AUTO = ("", "auto", None)


def abm_tick_zero(data: dict, base: Path | None = None, expand=None):
    """Calendar date of the ABM's tick 0, or None if the config gives none.

    `abm.tick_zero` wins; otherwise `tickZero` is read from the ABM's own
    config (`abm.config`, EpiHiper's config.json). Reading it rather than
    copying it keeps the date in one place -- the file the ABM itself runs on.
    """
    import datetime as _dt
    import json as _json

    abm = data.get("abm") or {}
    raw = abm.get("tick_zero")
    if raw in _AUTO:
        cfg_path = abm.get("config")
        if cfg_path in _AUTO:
            return None
        p = Path(str((expand or (lambda x: x))(str(cfg_path)))).expanduser()
        if not p.is_absolute() and base is not None and not p.exists():
            p = base / p
        if not p.is_file():
            raise ConfigError(f"abm.config not found: {cfg_path}")
        js = _json.loads(p.read_text())
        raw = js.get("tickZero")
        if raw is None:
            raise ConfigError(f"{cfg_path} has no 'tickZero'; set abm.tick_zero instead")
    return _dt.date.fromisoformat(str(raw)[:10])


def _first_import_tick(path) -> int:
    """First tick in an ABM seed schedule (`tick,count` CSV)."""
    import csv
    with open(path, newline="") as fh:
        ticks = [int(float(r["tick"])) for r in csv.DictReader(fh)
                 if r.get("tick") not in (None, "") and float(r.get("count", 1) or 0) > 0]
    if not ticks:
        raise ConfigError(f"seed schedule {path} has no importations")
    return min(ticks)


def resolve_location(data: dict) -> None:
    """Default `genetic_painter.location` from the population block.

    `division` is the state name and `divisionAbbr` its postal code, which is
    `population.state` upper-cased -- a transform `{placeholder}` substitution
    cannot express, so it is done here. "auto" or empty means derive.
    """
    gp = data.get("genetic_painter")
    pop = data.get("population") or {}
    if not isinstance(gp, dict):
        return
    loc = gp.get("location")
    if loc is None:
        loc = gp["location"] = {}
    if not isinstance(loc, dict):
        return              # a JSON string: the caller's responsibility
    if loc.get("division") in _AUTO and pop.get("state_name"):
        loc["division"] = str(pop["state_name"])
    if loc.get("divisionAbbr") in _AUTO and pop.get("state"):
        loc["divisionAbbr"] = str(pop["state"]).upper()


def resolve_derived(data: dict, base: Path | None = None, expand=None) -> None:
    """Fill every derived setting: location, start_tick, start_date.

    `expand` resolves `{placeholder}` paths for callers (the Snakefile) that
    run this before their own expansion pass.
    """
    expand = expand or (lambda s: s)
    resolve_location(data)
    gp = data.get("genetic_painter")
    if isinstance(gp, dict) and gp.get("start_tick") in _AUTO:
        sched = (data.get("abm") or {}).get("seed_schedule")
        if sched in _AUTO:
            raise ConfigError(
                "genetic_painter.start_tick is 'auto' but abm.seed_schedule is not "
                "set; point it at the schedule the ABM consumed, or give a tick.")
        p = Path(str(expand(str(sched)))).expanduser()
        if not p.is_absolute() and base is not None and not p.exists():
            p = base / p
        if not p.is_file():
            raise ConfigError(f"abm.seed_schedule not found: {p}")
        gp["start_tick"] = _first_import_tick(p)
    resolve_calendar(data, base, expand)


def resolve_calendar(data: dict, base: Path | None = None, expand=None) -> None:
    """Fill or check `genetic_painter.start_date` against the ABM's tick 0.

    The start date is not an independent setting: it is tick 0 plus
    `start_tick`. Giving both by hand let them drift -- a Massachusetts config
    labelled tick 128 as 2021-05-23 when the simulation's calendar puts it at
    2021-04-07, silently shifting every painted and line-list date 46 days.

    So `start_date: auto` (or empty) is derived, and an explicit date that
    disagrees with the ABM is an error rather than a quiet relabelling.
    Without an ABM tick 0 the explicit date is used as before.
    """
    import datetime as _dt

    gp = data.get("genetic_painter")
    if not isinstance(gp, dict):
        return
    t0 = abm_tick_zero(data, base, expand)
    given = gp.get("start_date")
    tick = gp.get("start_tick")
    if t0 is None:
        if given in _AUTO:
            raise ConfigError(
                "genetic_painter.start_date is 'auto' but no ABM tick 0 is "
                "configured; set abm.config (or abm.tick_zero), or give a date.")
        return
    if tick is None:
        raise ConfigError("genetic_painter.start_tick is required to place start_date")
    derived = t0 + _dt.timedelta(days=int(tick))
    if given in _AUTO:
        gp["start_date"] = derived.isoformat()
        return
    if str(given)[:10] != derived.isoformat():
        raise ConfigError(
            f"genetic_painter.start_date {given} disagrees with the ABM calendar: "
            f"tick {tick} is {derived.isoformat()} (tick 0 = {t0.isoformat()}).\n"
            f"  Set start_date: auto to derive it, or change start_tick to "
            f"{(_dt.date.fromisoformat(str(given)[:10]) - t0).days} if "
            f"{given} is the date you meant.")


class Config:
    """A loaded PhyloGAS configuration."""

    def __init__(self, data: dict, path: Path | None = None):
        self.data = data
        self.path = path
        self._expand_all()
        resolve_derived(self.data, Path(path).parent if path else None)

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

    @staticmethod
    def resolve_variant(raw) -> "Path | None":
        """Find a configured path, tolerating compression-extension drift.

        A config may name ``va_person.csv.xz`` while the file on disk is the
        decompressed ``va_person.csv`` (or vice versa) -- common when inputs
        are copied or symlinked from a cluster rather than fetched. Every
        consumer reads both forms via pandas, so the checkers should too.

        Returns the first existing variant, or None.
        """
        if raw is None:
            return None
        p = Path(str(raw)).expanduser()
        if p.exists():
            return p
        # configured compressed -> try the decompressed form
        if p.suffix in (".xz", ".gz", ".bz2", ".zst"):
            plain = p.with_suffix("")
            if plain.exists():
                return plain
        else:
            # configured plain -> try the compressed forms
            for ext in (".xz", ".gz", ".bz2", ".zst"):
                cand = p.with_suffix(p.suffix + ext)
                if cand.exists():
                    return cand
        return None

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

    # ------------------------------------------------------- derived paths
    def ascertainment_outputs(self) -> dict:
        """Paths TwinSampler derives from ``ascertainment.output``.

        simulate_linelist.py strips the extension and appends a suffix, so
        these are deterministic rather than guessed::

            linelist.csv  ->  linelist_allevents.csv.xz
                              linelist_mugration.json
        """
        out = self.get("ascertainment.output", default=None)
        if not out:
            return {}
        base = re.sub(r"\.csv(\.gz|\.xz)?$", "", str(out))
        return {
            "allevents": f"{base}_allevents.csv.xz",
            "mugration": f"{base}_mugration.json",
        }

    def resolve_benchmark_input(self, key: str, kind: str):
        """Resolve a benchmark input that may be 'auto', 'none', or a path.

        Returns ``(path_or_None, reason)``. ``auto`` derives the path from
        ``ascertainment.output`` and keeps it only if the file is actually
        there, so a pipeline that has not produced it yet simply skips that
        benchmark rather than failing.
        """
        raw = self.get(key, default="auto")
        raw = "auto" if raw in (None, "") else str(raw)

        if raw.lower() in ("none", "skip", "off", "false"):
            return None, "disabled in config"

        if raw.lower() == "auto":
            derived = self.ascertainment_outputs().get(kind)
            if not derived:
                return None, "auto: ascertainment.output is not configured"
            found = self.resolve_variant(derived)
            if found is None:
                return None, (f"auto: not there yet; the pipeline's line-list "
                              f"step writes it ({derived})")
            return found, f"auto from ascertainment.output"

        found = self.resolve_variant(raw)
        if found is None:
            return None, f"configured path not found: {raw}"
        return found, "explicit"

    def results_dir(self, stage: str) -> Path:
        """Return (and create) the output directory for a named pipeline stage."""
        root = Path(str(self.get("results_dir", default="results")))
        d = root / stage
        d.mkdir(parents=True, exist_ok=True)
        return d

    def __repr__(self) -> str:
        name = self.get("project_name", default="<unnamed>")
        return f"<Config project_name={name!r} path={self.path}>"
