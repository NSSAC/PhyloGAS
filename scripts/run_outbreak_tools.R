#!/usr/bin/env Rscript
# Run Bayesian outbreak / transmission reconstruction tools on one prepared
# PhyloGAS focal outbreak directory:
#   1. outbreaker2  (Campbell et al. 2018; skipped gracefully if not installed)
#   2. TransPhylo   (Didelot et al. 2017; uses TreeTime dated_tree_years.nwk)
#   3. JUNIPER      (Specht et al. 2025; broadinstitute/juniper0 consensus mode)
#
# IMPORTANT: `juniper0` exports an `initialize()` function that masks
# `methods::initialize`, so `juniper0` MUST be attached last after `outbreaker2`
# and `TransPhylo` have already completed.
#
# Usage:
#   Rscript scripts/run_outbreak_tools.R --workdir <outbreak_dir>
#       [--tools outbreaker2,transphylo,juniper]
#       [--juniper-global 1500] [--juniper-local 100] [--juniper-thin 15]
#       [--tp-iter 50000] [--tp-thin 50]
#       [--ob-iter 20000] [--ob-thin 40]
#       [--seed 42]

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL, as_fn = as.character) {
  i <- match(flag, args)
  if (is.na(i) || i >= length(args)) return(default)
  as_fn(args[i + 1])
}

workdir        <- getarg("--workdir")
tools_raw      <- getarg("--tools", "outbreaker2,transphylo,juniper")
juniper_global <- getarg("--juniper-global", 1500L, as.integer)
juniper_local  <- getarg("--juniper-local", 100L, as.integer)
juniper_thin   <- getarg("--juniper-thin", 15L, as.integer)
tp_iter        <- getarg("--tp-iter", 50000L, as.integer)
tp_thin        <- getarg("--tp-thin", 50L, as.integer)
ob_iter        <- getarg("--ob-iter", 20000L, as.integer)
ob_thin        <- getarg("--ob-thin", 40L, as.integer)
seed           <- getarg("--seed", 42L, as.integer)

if (is.null(workdir) || !nzchar(workdir)) {
  stop("--workdir <outbreak_dir> is required")
}
workdir <- normalizePath(workdir, mustWork = TRUE)
setwd(workdir)

# Honour colon-separated R_LIBS_USER so a user library layered over a shared
# group library (e.g. ~/Rlibs:/project/bii_nssac/people/bl4zc/Rlibs/juniper_R4.6)
# resolves packages from both.
raw_libs <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(raw_libs)) {
  parts <- strsplit(raw_libs, ":", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts) & dir.exists(parts)]
  if (length(parts) > 0) .libPaths(unique(c(parts, .libPaths())))
}

want_tools <- trimws(tolower(strsplit(tools_raw, ",", fixed = TRUE)[[1]]))
set.seed(seed)

# Simple base-R ESS estimator (avoids requiring `coda`)
calc_ess <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 10 || stats::sd(x) == 0) return(NA_real_)
  ac <- stats::acf(x, lag.max = min(50L, n - 2L), plot = FALSE)$acf[-1]
  pos <- which(ac <= 0)
  k <- if (length(pos) > 0) pos[1] - 1L else length(ac)
  tau <- if (k >= 1L) 1.0 + 2.0 * sum(ac[seq_len(k)]) else 1.0
  max(1.0, min(as.numeric(n), n / max(1.0, tau)))
}

write_json_simple <- function(lst, path) {
  fmt_val <- function(v) {
    if (is.null(v) || length(v) == 0 || is.na(v[1])) return("null")
    if (is.logical(v[1])) return(ifelse(v[1], "true", "false"))
    if (is.numeric(v[1])) {
      if (!is.finite(v[1])) return("null")
      return(format(v[1], digits = 8, scientific = FALSE))
    }
    sprintf('"%s"', gsub('"', '\\"', as.character(v[1]), fixed = TRUE))
  }
  keys <- names(lst)
  lines <- vapply(seq_along(lst), function(i) {
    sprintf('  "%s": %s', keys[i], fmt_val(lst[[i]]))
  }, character(1))
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), con = path)
}

write_donor_probs <- function(sample_names, direct_mat, direct_unobs,
                              indirect_mat, indirect_unobs, n_post, outfile) {
  # direct_mat[i, j] = count of post-burnin states where sample i directly infected j
  n_s <- length(sample_names)
  rows <- list()
  idx_r <- 0L
  for (j in seq_len(n_s)) {
    rec <- sample_names[j]
    for (i in seq_len(n_s)) {
      if (i == j) next
      dp <- direct_mat[i, j] / n_post
      ip <- indirect_mat[i, j] / n_post
      if (dp > 0 || ip > 0) {
        idx_r <- idx_r + 1L
        rows[[idx_r]] <- data.frame(
          recipient_id = rec, donor_id = sample_names[i],
          direct_prob = round(dp, 6), indirect_prob = round(ip, 6),
          stringsAsFactors = FALSE
        )
      }
    }
    dp_u <- direct_unobs[j] / n_post
    ip_u <- indirect_unobs[j] / n_post
    if (dp_u > 0 || ip_u > 0) {
      idx_r <- idx_r + 1L
      rows[[idx_r]] <- data.frame(
        recipient_id = rec, donor_id = "UNOBSERVED",
        direct_prob = round(dp_u, 6), indirect_prob = round(ip_u, 6),
        stringsAsFactors = FALSE
      )
    }
  }
  out_df <- if (length(rows) > 0) do.call(rbind, rows) else data.frame(
    recipient_id = character(0), donor_id = character(0),
    direct_prob = numeric(0), indirect_prob = numeric(0)
  )
  utils::write.csv(out_df, outfile, row.names = FALSE, quote = FALSE)
}

# Read outbreak_manifest.json for simulation-matched priors
manifest_txt <- paste(readLines("outbreak_manifest.json", warn = FALSE), collapse = "\n")
extract_num <- function(key, default) {
  m <- regmatches(manifest_txt, regexpr(sprintf('"%s"\\s*:\\s*[-+0-9.eE]+', key), manifest_txt))
  if (length(m) == 0) return(default)
  as.numeric(sub(sprintf('"%s"\\s*:\\s*', key), "", m))
}

init_mu  <- extract_num("init_mu_day", 1.873e-6)
a_g      <- extract_num("a_g", 9.968)
lambda_g <- extract_num("lambda_g", 2.133)
a_s      <- extract_num("a_s", 37.4)
lambda_s <- extract_num("lambda_s", 9.4)
init_pi  <- 0.45

meta_df <- utils::read.csv("input_data/metadata.csv", stringsAsFactors = FALSE)
n_obs_samples <- nrow(meta_df)
cat(sprintf("=== Outbreak tools runner in %s (N = %d samples) ===\n", workdir, n_obs_samples))
cat(sprintf("  Priors: init_mu=%.4g/day, GI Gamma(%.2f, rate=%.2f), Soj Gamma(%.2f, rate=%.2f)\n",
            init_mu, a_g, lambda_g, a_s, lambda_s))

# ===========================================================================
# 1. outbreaker2
# ===========================================================================
if ("outbreaker2" %in% want_tools) {
  if (!requireNamespace("outbreaker2", quietly = TRUE)) {
    cat("[outbreaker2] Package not installed in .libPaths(); skipping gracefully.\n")
  } else {
    cat("\n[outbreaker2] Starting MCMC...\n")
    tryCatch({
      t0 <- Sys.time()
      dna <- ape::read.dna("input_data/aligned.fasta", format = "fasta")
      # Align metadata row order to FASTA labels
      ord <- match(rownames(dna), meta_df$sample)
      if (any(is.na(ord))) stop("FASTA labels do not match metadata.csv sample column")
      m_ord <- meta_df[ord, , drop = FALSE]
      s_names <- m_ord$sample
      # outbreaker2 requires DNAbin rownames to be "1".."N" (matching integer case indices)
      rownames(dna) <- as.character(seq_len(nrow(dna)))
      dates_d <- as.Date(m_ord$date)
      d0 <- min(dates_d)
      int_dates <- as.integer(dates_d - d0)

      k_grid <- 1:40
      w_dens <- stats::pgamma(k_grid, shape = a_g, rate = lambda_g) -
                stats::pgamma(k_grid - 1, shape = a_g, rate = lambda_g)
      w_dens <- pmax(w_dens, 1e-12); w_dens <- w_dens / sum(w_dens)
      f_dens <- stats::pgamma(k_grid, shape = a_s, rate = lambda_s) -
                stats::pgamma(k_grid - 1, shape = a_s, rate = lambda_s)
      f_dens <- pmax(f_dens, 1e-12); f_dens <- f_dens / sum(f_dens)

      ob_dat <- outbreaker2::outbreaker_data(
        dates = int_dates, dna = dna, w_dens = w_dens, f_dens = f_dens
      )
      ob_cfg <- outbreaker2::create_config(
        n_iter = ob_iter, sample_every = ob_thin,
        init_mu = init_mu, init_pi = init_pi,
        move_kappa = TRUE, max_kappa = 5L
      )
      ob_res <- outbreaker2::outbreaker(data = ob_dat, config = ob_cfg)
      t_ob <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      n_tot <- nrow(ob_res)
      keep <- seq(max(1L, floor(n_tot / 2) + 1L), n_tot)
      n_post <- length(keep)
      n_s <- length(s_names)

      dir_mat <- matrix(0.0, n_s, n_s)
      ind_mat <- matrix(0.0, n_s, n_s)
      dir_unobs <- numeric(n_s)
      ind_unobs <- numeric(n_s)
      tinf_med <- numeric(n_s)

      for (j in seq_len(n_s)) {
        al_col <- ob_res[[paste0("alpha_", j)]][keep]
        kp_col <- ob_res[[paste0("kappa_", j)]][keep]
        ti_col <- ob_res[[paste0("t_inf_", j)]][keep]
        tinf_med[j] <- stats::median(ti_col, na.rm = TRUE)

        for (s_i in seq_along(al_col)) {
          anc <- al_col[s_i]
          kap <- kp_col[s_i]
          if (is.na(anc) || anc < 1L || anc > n_s) {
            dir_unobs[j] <- dir_unobs[j] + 1.0
            ind_unobs[j] <- ind_unobs[j] + 1.0
          } else {
            ind_mat[anc, j] <- ind_mat[anc, j] + 1.0
            if (!is.na(kap) && kap == 1L) {
              dir_mat[anc, j] <- dir_mat[anc, j] + 1.0
            } else {
              dir_unobs[j] <- dir_unobs[j] + 1.0
            }
          }
        }
      }

      write_donor_probs(s_names, dir_mat, dir_unobs, ind_mat, ind_unobs,
                        n_post, "outbreaker2_donor_probs.csv")

      # Consensus MPC network (breaking any cycles greedily by lowest credibility)
      mpc_parent <- rep(NA_integer_, n_s)
      mpc_prob <- numeric(n_s)
      for (j in seq_len(n_s)) {
        best_i <- which.max(ind_mat[, j])
        if (length(best_i) == 1L && ind_mat[best_i, j] >= ind_unobs[j]) {
          mpc_parent[j] <- best_i
          mpc_prob[j] <- ind_mat[best_i, j] / n_post
        }
      }
      # Break cycles if present
      for (j in seq_len(n_s)) {
        cur <- mpc_parent[j]
        seen <- j
        while (!is.na(cur)) {
          if (cur %in% seen) {
            weakest <- seen[which.min(mpc_prob[seen])]
            mpc_parent[weakest] <- NA_integer_
            break
          }
          seen <- c(seen, cur)
          cur <- mpc_parent[cur]
        }
      }

      nodes_ob <- data.frame(
        idx = seq_len(n_s),
        infector_idx = mpc_parent,
        infection_date = as.character(d0 + round(tinf_med)),
        observed = TRUE,
        sample = s_names,
        stringsAsFactors = FALSE
      )
      utils::write.csv(nodes_ob, "outbreaker2_nodes.csv", row.names = FALSE, na = "")

      mu_post <- ob_res$mu[keep]
      pi_post <- ob_res$pi[keep]
      write_json_simple(list(
        tool = "outbreaker2",
        runtime_secs = round(t_ob, 2),
        n_post_samples = n_post,
        mu_est_day = stats::median(mu_post, na.rm = TRUE),
        mu_est_yr = stats::median(mu_post, na.rm = TRUE) * 365.25,
        pi_est = stats::median(pi_post, na.rm = TRUE),
        R_est = NA_real_,
        ess_min = round(min(calc_ess(ob_res$post[keep]), calc_ess(mu_post), na.rm = TRUE), 1)
      ), "outbreaker2_params.json")
      cat(sprintf("[outbreaker2] Done in %.1fs\n", t_ob))
    }, error = function(e) {
      cat(sprintf("[outbreaker2] ERROR: %s\n", conditionMessage(e)))
    })
  }
}

# ===========================================================================
# 2. TransPhylo (on TreeTime dated_tree_years.nwk)
# ===========================================================================
if ("transphylo" %in% want_tools) {
  nwk_file <- "treetime/dated_tree_years.nwk"
  meta_json <- "treetime/dated_tree_meta.json"
  if (!requireNamespace("TransPhylo", quietly = TRUE)) {
    cat("[TransPhylo] Package not installed in .libPaths(); skipping.\n")
  } else if (!file.exists(nwk_file)) {
    cat(sprintf("[TransPhylo] Missing %s; skipping.\n", nwk_file))
  } else {
    cat("\n[TransPhylo] Starting MCMC on TreeTime dated tree...\n")
    tryCatch({
      t0 <- Sys.time()
      tr <- ape::read.tree(nwk_file)
      if (!ape::is.binary(tr)) tr <- ape::multi2di(tr)
      tr$edge.length[tr$edge.length <= 1e-6] <- 1e-5

      meta_tt_txt <- paste(readLines(meta_json, warn = FALSE), collapse = "\n")
      m_dl <- regmatches(meta_tt_txt, regexpr('"date_last_sample"\\s*:\\s*[0-9.]+', meta_tt_txt))
      date_last <- if (length(m_dl) > 0) as.numeric(sub('"date_last_sample"\\s*:\\s*', "", m_dl)) else 2021.75

      ptree <- TransPhylo::ptreeFromPhylo(tr, dateLastSample = date_last)
      tp_res <- TransPhylo::inferTTree(
        ptree,
        mcmcIterations = tp_iter,
        thinning = tp_thin,
        w.shape = a_g,
        w.scale = 1.0 / lambda_g,
        ws.shape = a_s,
        ws.scale = 1.0 / lambda_s,
        startPi = init_pi,
        updatePi = TRUE,
        dateT = date_last + (7.0 / 365.25),
        verbose = FALSE
      )
      t_tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      n_rec <- length(tp_res)
      keep <- seq(max(1L, floor(n_rec / 2) + 1L), n_rec)
      n_post <- length(keep)

      tt_first <- TransPhylo::extractTTree(tp_res[[keep[1]]]$ctree)
      s_names <- tt_first$nam
      n_s <- length(s_names)

      dir_mat <- matrix(0.0, n_s, n_s)
      ind_mat <- matrix(0.0, n_s, n_s)
      dir_unobs <- numeric(n_s)
      ind_unobs <- numeric(n_s)
      pi_vec <- numeric(n_post)
      r_vec  <- numeric(n_post)
      p_log  <- numeric(n_post)

      for (idx_k in seq_along(keep)) {
        rec_k <- tp_res[[keep[idx_k]]]
        pi_vec[idx_k] <- rec_k$pi
        r_vec[idx_k]  <- rec_k$off.r
        p_log[idx_k]  <- rec_k$pTTree + rec_k$pPTree

        tt <- TransPhylo::extractTTree(rec_k$ctree)$ttree
        parents <- tt[, 1]
        n_hosts <- nrow(tt)

        for (j in seq_len(n_s)) {
          p1 <- parents[j]
          if (!is.na(p1) && p1 >= 1L && p1 <= n_s) {
            dir_mat[p1, j] <- dir_mat[p1, j] + 1.0
          } else {
            dir_unobs[j] <- dir_unobs[j] + 1.0
          }
          cur <- p1
          steps <- 0L
          anc_found <- 0L
          while (!is.na(cur) && cur >= 1L && cur <= n_hosts && steps <= n_hosts) {
            if (cur <= n_s) {
              anc_found <- cur
              break
            }
            cur <- parents[cur]
            steps <- steps + 1L
          }
          if (anc_found >= 1L) {
            ind_mat[anc_found, j] <- ind_mat[anc_found, j] + 1.0
          } else {
            ind_unobs[j] <- ind_unobs[j] + 1.0
          }
        }
      }

      write_donor_probs(s_names, dir_mat, dir_unobs, ind_mat, ind_unobs,
                        n_post, "transphylo_donor_probs.csv")

      med_ctree <- TransPhylo::medTTree(tp_res, burnin = 0.5)
      tt_med <- TransPhylo::extractTTree(med_ctree)
      n_m <- nrow(tt_med$ttree)
      nodes_tp <- data.frame(
        idx = seq_len(n_m),
        infector_idx = ifelse(tt_med$ttree[, 1] == 0, NA_integer_, as.integer(tt_med$ttree[, 1])),
        infection_num_date = round(tt_med$ttree[, 2], 8),
        observed = seq_len(n_m) <= n_s,
        sample = c(tt_med$nam, rep("", max(0L, n_m - n_s)))[seq_len(n_m)],
        stringsAsFactors = FALSE
      )
      utils::write.csv(nodes_tp, "transphylo_nodes.csv", row.names = FALSE, na = "")

      write_json_simple(list(
        tool = "TransPhylo",
        runtime_secs = round(t_tp, 2),
        n_post_samples = n_post,
        mu_est_day = NA_real_,
        mu_est_yr = NA_real_,
        pi_est = stats::median(pi_vec, na.rm = TRUE),
        R_est = stats::median(r_vec, na.rm = TRUE),
        n_hosts_med = n_m,
        n_unsampled_med = max(0L, n_m - n_s),
        ess_min = round(min(calc_ess(p_log), calc_ess(pi_vec), calc_ess(r_vec), na.rm = TRUE), 1)
      ), "transphylo_params.json")
      cat(sprintf("[TransPhylo] Done in %.1fs (pi_hat=%.3f, R_hat=%.2f)\n",
                  t_tp, stats::median(pi_vec), stats::median(r_vec)))
    }, error = function(e) {
      cat(sprintf("[TransPhylo] ERROR: %s\n", conditionMessage(e)))
    })
  }
}

# ===========================================================================
# 3. JUNIPER (`juniper0` — loaded LAST because it masks `methods::initialize`)
# ===========================================================================
if ("juniper" %in% want_tools) {
  if (!requireNamespace("juniper0", quietly = TRUE)) {
    cat("[JUNIPER] Package `juniper0` not installed in .libPaths(); skipping.\n")
  } else {
    cat("\n[JUNIPER] Attaching juniper0 and starting MCMC...\n")
    tryCatch({
      suppressPackageStartupMessages(library(juniper0))
      cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
      if (is.na(cpus) || cpus < 1L) cpus <- 4L
      subtrees <- max(1L, min(cpus, floor(n_obs_samples / 25L)))

      cat(sprintf("[JUNIPER] samples=%d, subtrees=%d, n_global=%d, n_local=%d, sample_every=%d\n",
                  n_obs_samples, subtrees, juniper_global, juniper_local, juniper_thin))

      t0 <- Sys.time()
      init <- initialize(
        n_subtrees   = subtrees,
        n_global     = juniper_global,
        n_local      = juniper_local,
        sample_every = juniper_thin,
        init_mu      = init_mu,
        a_g          = a_g,
        lambda_g     = lambda_g,
        a_s          = a_s,
        lambda_s     = lambda_s,
        init_pi      = init_pi,
        indir        = "input_data"
      )
      res <- run_mcmc(init, noisy = FALSE)
      t_jun <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      save(res, file = "juniper_results.RData")

      loglik_all <- res[[1]]
      states_all <- res[[2]]
      names_obs  <- res[[3]]
      s_max      <- as.Date(res[[5]])
      n_obs_j    <- length(names_obs)

      # Identify real observed sample indices (excluding JUNIPER's synthetic "ref_genome" at idx 1)
      real_mask <- names_obs != "ref_genome" & nzchar(names_obs)
      real_indices <- which(real_mask)
      s_names <- names_obs[real_indices]
      n_s <- length(s_names)
      obs_to_pos <- rep(0L, n_obs_j)
      obs_to_pos[real_indices] <- seq_len(n_s)

      n_states <- length(states_all)
      keep <- seq(max(1L, floor(n_states / 2) + 1L), n_states)
      n_post <- length(keep)

      dir_mat <- matrix(0.0, n_s, n_s)
      ind_mat <- matrix(0.0, n_s, n_s)
      dir_unobs <- numeric(n_s)
      ind_unobs <- numeric(n_s)
      mu_vec <- rep(NA_real_, n_post)
      pi_vec <- rep(NA_real_, n_post)
      r_vec  <- rep(NA_real_, n_post)

      for (idx_k in seq_along(keep)) {
        st <- states_all[[keep[idx_k]]]
        h <- st$h
        n_h <- if (!is.null(st$n)) st$n else length(h)
        if (!is.null(st$mu)) mu_vec[idx_k] <- as.numeric(st$mu[1])
        if (!is.null(st$pi)) pi_vec[idx_k] <- as.numeric(st$pi[1])
        if (!is.null(st$R))  r_vec[idx_k]  <- as.numeric(st$R[1])

        for (pos_j in seq_len(n_s)) {
          j_host <- real_indices[pos_j]
          p1 <- if (j_host <= length(h)) h[j_host] else NA_integer_
          if (!is.na(p1) && p1 >= 1L && p1 <= n_obs_j && obs_to_pos[p1] > 0L) {
            dir_mat[obs_to_pos[p1], pos_j] <- dir_mat[obs_to_pos[p1], pos_j] + 1.0
          } else {
            dir_unobs[pos_j] <- dir_unobs[pos_j] + 1.0
          }

          cur <- p1
          steps <- 0L
          anc_pos <- 0L
          while (!is.na(cur) && cur >= 1L && cur <= length(h) && steps <= n_h) {
            if (cur <= n_obs_j && obs_to_pos[cur] > 0L) {
              anc_pos <- obs_to_pos[cur]
              break
            }
            cur <- h[cur]
            steps <- steps + 1L
          }
          if (anc_pos > 0L) {
            ind_mat[anc_pos, pos_j] <- ind_mat[anc_pos, pos_j] + 1.0
          } else {
            ind_unobs[pos_j] <- ind_unobs[pos_j] + 1.0
          }
        }
      }

      # Fallback parameter extraction from summarize() if state slots use different names
      if (all(is.na(mu_vec)) || all(is.na(pi_vec))) {
        try({
          s_obj <- summarize(res)
          if (is.list(s_obj)) {
            if (all(is.na(mu_vec)) && !is.null(s_obj$mu)) mu_vec <- as.numeric(s_obj$mu)
            if (all(is.na(pi_vec)) && !is.null(s_obj$pi)) pi_vec <- as.numeric(s_obj$pi)
            if (all(is.na(r_vec))  && !is.null(s_obj$R))  r_vec  <- as.numeric(s_obj$R)
          }
        }, silent = TRUE)
      }

      write_donor_probs(s_names, dir_mat, dir_unobs, ind_mat, ind_unobs,
                        n_post, "juniper_donor_probs.csv")

      # Export MAP post-burnin state to `juniper_nodes.csv`
      ll_keep <- if (length(loglik_all) >= max(keep)) loglik_all[keep] else rep(0, n_post)
      best_state <- states_all[[keep[which.max(ll_keep)]]]
      h_b <- best_state$h
      seqs_b <- best_state$seq
      n_b <- if (!is.null(best_state$n)) best_state$n else length(h_b)
      inf_off <- vapply(seq_len(n_b), function(j) {
        if (j <= length(seqs_b) && length(seqs_b[[j]]) >= 1L) seqs_b[[j]][1] else NA_real_
      }, numeric(1))

      nodes_jun <- data.frame(
        idx = seq_len(n_b),
        infector_idx = ifelse(seq_len(n_b) <= length(h_b), h_b[seq_len(n_b)], NA_integer_),
        infection_date = as.character(s_max + inf_off),
        observed = seq_len(n_b) <= n_obs_j,
        sample = c(names_obs, rep("", max(0L, n_b - n_obs_j)))[seq_len(n_b)],
        stringsAsFactors = FALSE
      )
      utils::write.csv(nodes_jun, "juniper_nodes.csv", row.names = FALSE, na = "")

      mu_med <- if (any(!is.na(mu_vec))) stats::median(mu_vec, na.rm = TRUE) else NA_real_
      pi_med <- if (any(!is.na(pi_vec))) stats::median(pi_vec, na.rm = TRUE) else NA_real_
      r_med  <- if (any(!is.na(r_vec)))  stats::median(r_vec, na.rm = TRUE)  else NA_real_

      write_json_simple(list(
        tool = "JUNIPER",
        runtime_secs = round(t_jun, 2),
        n_post_samples = n_post,
        subtrees = subtrees,
        mu_est_day = mu_med,
        mu_est_yr = if (!is.na(mu_med)) mu_med * 365.25 else NA_real_,
        pi_est = pi_med,
        R_est = r_med,
        n_hosts_map = n_b,
        n_unsampled_map = max(0L, n_b - n_obs_j),
        ess_min = round(calc_ess(ll_keep), 1)
      ), "juniper_params.json")

      # Clean up JUNIPER scratch tmp/ directory
      if (dir.exists("tmp")) unlink("tmp", recursive = TRUE)
      cat(sprintf("[JUNIPER] Done in %.1fs (%.2fs/global iter)\n",
                  t_jun, t_jun / max(1L, juniper_global)))
    }, error = function(e) {
      cat(sprintf("[JUNIPER] ERROR: %s\n", conditionMessage(e)))
    })
  }
}

cat("\n=== All requested outbreak tools finished ===\n")
