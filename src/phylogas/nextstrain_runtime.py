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
    # Check common shared-storage Nextstrain home on Rivanna (`/project/biocomplexity/$USER/nextstrain`)
    # before `~/.nextstrain`, since `~/.nextstrain` often holds an older container config
    # when the active conda runtime was relocated off `$HOME` quota.
    user = os.environ.get("USER") or ""
    if user:
        proj_ns = Path(f"/project/biocomplexity/{user}/nextstrain")
        if proj_ns.is_dir():
            yield proj_ns
    yield Path("~/.nextstrain").expanduser()


def runtime_bin(tool: str, ncov_dir=None) -> "Path | None":
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
    if ncov_dir is not None:
        snk_conda = Path(ncov_dir) / ".snakemake" / "conda"
        if snk_conda.is_dir():
            for env_dir in sorted(snk_conda.iterdir()):
                cand = env_dir / "bin" / tool
                if cand.is_file() and os.access(cand, os.X_OK):
                    return cand
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

    cmd, source = tool_command(tool, ncov_dir=ncov_dir)
    if cmd is not None:
        # When invoking `augur` from a Nextstrain conda runtime's `bin/`, prepend
        # that `bin/` directory to PATH so child tools (`iqtree3`, `iqtree2`,
        # `fasttree`, `mafft`) that `augur tree` shells out to are found too.
        run_kw = dict(kw)
        if source == "runtime":
            bin_dir = str(Path(cmd[0]).parent)
            env_map = dict(run_kw.get("env") or os.environ)
            cur_path = env_map.get("PATH", "")
            env_map["PATH"] = f"{bin_dir}:{cur_path}" if cur_path else bin_dir
            run_kw["env"] = env_map
        return subprocess.run(cmd + args, **run_kw), source

    if ncov_dir is not None and shutil.which("nextstrain"):
        # Quote every argument: this becomes a shell line, and paths from a
        # config can contain spaces.
        import shlex
        line = " ".join(shlex.quote(a) for a in [tool, *args])
        shell_kw = dict(kw)
        shell_kw.pop("input", None)
        shell_kw.pop("cwd", None)

        # Try default runtime first ONLY if docker binary exists or singularity isn't needed;
        # on HPC clusters (e.g. Rivanna) where `~/.nextstrain/config` may point to a
        # retired `/opt/apptainer/current/bin/docker` wrapper while `apptainer`/`singularity`
        # is loaded via `module load apptainer`, try `--singularity` as automatic fallback
        # (or primary when `docker` is absent and `singularity`/`apptainer` is on PATH).
        flag_attempts: list[tuple[list[str], str]] = []
        has_singularity = bool(shutil.which("singularity") or shutil.which("apptainer"))
        has_docker = bool(shutil.which("docker"))
        if has_singularity and not has_docker:
            flag_attempts.append((["--singularity"], "nextstrain shell --singularity"))
            flag_attempts.append(([], "nextstrain shell"))
        else:
            flag_attempts.append(([], "nextstrain shell"))
            if has_singularity:
                flag_attempts.append((["--singularity"], "nextstrain shell --singularity"))
        flag_attempts.append((["--conda"], "nextstrain shell --conda"))

        last_proc = None
        last_src = "nextstrain shell"
        for rt_flags, label in flag_attempts:
            proc = subprocess.run(
                ["nextstrain", "shell", *rt_flags, str(ncov_dir)],
                input=line + "\n",
                **shell_kw,
            )
            last_proc, last_src = proc, label
            if proc.returncode == 0:
                return proc, label
        return last_proc, last_src

    return None, source


def tool_command(tool: str, ncov_dir=None) -> "tuple[list, str] | tuple[None, str]":
    """How to run `tool`, and a one-word label for where it came from.

    Returns ([argv prefix], source) or (None, reason). Order: the ambient
    environment first (so an explicit install wins), then the runtime.
    """
    direct = shutil.which(tool)
    if direct:
        return [direct], "PATH"

    inside = runtime_bin(tool, ncov_dir=ncov_dir)
    if inside is not None:
        return [str(inside)], "runtime"

    homes = ", ".join(str(h) for h in dict.fromkeys(_homes()))
    return None, (f"{tool} is neither on PATH nor in a Nextstrain conda runtime "
                  f"(looked under {homes}; set NEXTSTRAIN_HOME if it is elsewhere), "
                  f"and the Nextstrain CLI is not available to run it")
