"""Locating tools that live inside the Nextstrain runtime.

augur and nextclade are installed in the Nextstrain runtime, not in the
PhyloGAS environment, so a plain `subprocess.run(["augur", ...])` fails with
"not on PATH" even when `nextstrain build` works perfectly.

`--exec` is the documented way to run one command inside the runtime, but it
belongs to newer releases; an older CLI has no such flag, and putting it on
`nextstrain build` gets it forwarded to Snakemake, which rejects it as an
invalid `--executor` choice. So two routes are tried, in this order:

1. The tool's own path, on PATH or in the conda runtime's bin directory:

       $NEXTSTRAIN_HOME/runtimes/conda/env/bin/<tool>

   (NEXTSTRAIN_HOME defaults to ~/.nextstrain). Direct, and no dependence on
   which flags the installed CLI happens to have.

2. Piped into `nextstrain shell <dir>`, which is not interactive-only: the
   shell it starts reads stdin when stdin is not a terminal. Slower and it
   prints a banner, but it is the public interface and it is the only route
   that works for the docker and singularity runtimes, where there is no
   bin directory on the host to point at.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

# Runtimes that put their tools in a directory we can reach directly. The
# container runtimes (docker, singularity) do not, and need the CLI.
_RUNTIME_BINS = ("runtimes/conda/env/bin",)


def _homes():
    env = os.environ.get("NEXTSTRAIN_HOME")
    if env:
        yield Path(env).expanduser()
    yield Path("~/.nextstrain").expanduser()


def runtime_bin(tool: str) -> "Path | None":
    """The tool's path inside the Nextstrain runtime, if it is installed."""
    seen = set()
    for home in _homes():
        if home in seen:
            continue
        seen.add(home)
        for rel in _RUNTIME_BINS:
            candidate = home / rel / tool
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return candidate
    return None


def run(tool: str, args, ncov_dir=None, **kw):
    """Run `tool` with `args` inside the Nextstrain runtime.

    Returns (CompletedProcess, source) or (None, reason). `kw` goes to
    subprocess.run; capture_output and text default on.
    """
    import subprocess

    kw.setdefault("capture_output", True)
    kw.setdefault("text", True)
    args = [str(a) for a in args]

    cmd, source = tool_command(tool)
    if cmd is not None:
        return subprocess.run(cmd + args, **kw), source

    if ncov_dir is not None and shutil.which("nextstrain"):
        # Quote every argument: this becomes a shell line, and paths from a
        # config can contain spaces.
        import shlex
        line = " ".join(shlex.quote(a) for a in [tool, *args])
        kw.pop("input", None)
        return (subprocess.run(["nextstrain", "shell", str(ncov_dir)],
                               input=line + "\n", **kw),
                "nextstrain shell")

    return None, source


def tool_command(tool: str, ncov_dir=None) -> "tuple[list, str] | tuple[None, str]":
    """How to run `tool`, and a one-word label for where it came from.

    Returns ([argv prefix], source) or (None, reason). Order: the ambient
    environment first (so an explicit install wins), then the runtime.
    """
    direct = shutil.which(tool)
    if direct:
        return [direct], "PATH"

    inside = runtime_bin(tool)
    if inside is not None:
        return [str(inside)], "runtime"

    homes = ", ".join(str(h) for h in dict.fromkeys(_homes()))
    return None, (f"{tool} is neither on PATH nor in a Nextstrain conda runtime "
                  f"(looked under {homes}; set NEXTSTRAIN_HOME if it is elsewhere), "
                  f"and the Nextstrain CLI is not available to run it")
