#!/usr/bin/env bash
# Prepare & verify R 4.6 outbreak reconstruction packages on a Rivanna login node.
#
# Reuses the pre-compiled `juniper0`, `TransPhylo`, `ape`, `igraph`, and `Rcpp`
# stack in `/project/bii_nssac/people/bl4zc/Rlibs/juniper_R4.6` when readable,
# and installs `outbreaker2` (plus any missing packages) into your personal
# `$HOME/Rlibs/phylogas_R4.6` directory (compute nodes have no outbound internet).
#
# Usage (run once on a Rivanna login node before submitting sbatch):
#   bash scripts/setup_outbreak_rlibs.sh

set -euo pipefail

USER_RLIB="${HOME}/Rlibs/phylogas_R4.6"
SHARED_RLIB="/project/bii_nssac/people/bl4zc/Rlibs/juniper_R4.6"
mkdir -p "${USER_RLIB}"

if command -v module >/dev/null 2>&1; then
    module purge
    module load gcc/14.2.0 openmpi/5.0.7 R/4.6.0
fi

if [[ -d "${SHARED_RLIB}" && -r "${SHARED_RLIB}" ]]; then
    export R_LIBS_USER="${USER_RLIB}:${SHARED_RLIB}"
    echo "Using layered R_LIBS_USER:"
    echo "  user  : ${USER_RLIB}"
    echo "  shared: ${SHARED_RLIB}"
else
    export R_LIBS_USER="${USER_RLIB}"
    echo "Shared library ${SHARED_RLIB} not readable; installing all packages into ${USER_RLIB}"
fi

Rscript -e '
raw_libs <- Sys.getenv("R_LIBS_USER", unset = "")
parts <- strsplit(raw_libs, ":", fixed = TRUE)[[1]]
parts <- parts[nzchar(parts) & dir.exists(parts)]
.libPaths(unique(c(parts, .libPaths())))
user_lib <- parts[1]

options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)

need_cran <- c("ape", "igraph", "Rcpp", "RcppArmadillo", "extraDistr", "remotes", "outbreaker2")
missing_cran <- need_cran[!vapply(need_cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  cat("Installing missing CRAN packages into", user_lib, ":", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, lib = user_lib)
}

if (!requireNamespace("TransPhylo", quietly = TRUE)) {
  cat("Installing archived TransPhylo 1.4.5 from CRAN archive...\n")
  install.packages(
    "https://cloud.r-project.org/src/contrib/Archive/TransPhylo/TransPhylo_1.4.5.tar.gz",
    repos = NULL, type = "source", lib = user_lib
  )
}

if (!requireNamespace("juniper0", quietly = TRUE)) {
  cat("Installing broadinstitute/juniper0@c6efa8a...\n")
  remotes::install_github(
    "broadinstitute/juniper0@c6efa8ac6b855abdced4c363b7a7364f6111ea82",
    lib = user_lib, upgrade = "never"
  )
}

cat("\n=== R Package Verification ===\n")
for (pkg in c("ape", "TransPhylo", "outbreaker2", "juniper0")) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  ver <- if (ok) as.character(packageVersion(pkg)) else "MISSING"
  cat(sprintf("  %-14s : %s\n", pkg, ver))
}
'
