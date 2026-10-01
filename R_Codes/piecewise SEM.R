options(max.print = 1e7)
options(width = 250)
options(warn = 1)
options(scipen = 999)
options(stringsAsFactors = FALSE)

rm(list = ls())
load("R_Datasets/metrics.rdata")

if (requireNamespace("utils", quietly = TRUE)) {
  try(utils::assignInNamespace("print.console", NULL, ns = "utils"), silent = TRUE)
}

if (!requireNamespace("piecewiseSEM", quietly = TRUE)) install.packages("piecewiseSEM")
if (!requireNamespace("boot", quietly = TRUE)) install.packages("boot")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(piecewiseSEM)
library(boot)
library(openxlsx)

fit_model_1 <- function(dat) {
  piecewiseSEM::psem(
    lm(benthic_invertebrate_richness ~ fish_richness, dat),
    lm(insect_richness ~ benthic_invertebrate_richness, dat),
    lm(zooplankton_richness ~ insect_richness, dat),
    lm(fungi_richness ~ benthic_invertebrate_richness, dat),
    lm(fish_synchrony ~ fish_richness, dat),
    lm(benthic_invertebrate_synchrony ~ benthic_invertebrate_richness, dat),
    lm(insect_synchrony ~ insect_richness, dat),
    lm(zooplankton_synchrony ~ zooplankton_richness + benthic_invertebrate_richness, dat),
    lm(phytoplankton_synchrony ~ phytoplankton_richness + fish_synchrony + benthic_invertebrate_synchrony, dat),
    lm(fungi_synchrony ~ fungi_richness + benthic_invertebrate_synchrony, dat),
    lm(bacteria_synchrony ~ bacteria_richness + benthic_invertebrate_synchrony, dat),
    lm(fish_stability ~ fish_synchrony, dat),
    lm(benthic_invertebrate_stability ~ benthic_invertebrate_synchrony, dat),
    lm(insect_stability ~ insect_synchrony, dat),
    lm(zooplankton_stability ~ zooplankton_synchrony, dat),
    lm(phytoplankton_stability ~ phytoplankton_synchrony, dat),
    lm(fungi_stability ~ fungi_synchrony, dat),
    lm(bacteria_stability ~ bacteria_synchrony, dat),
    lm(multigroups_stability ~ fish_synchrony + benthic_invertebrate_synchrony + insect_synchrony + zooplankton_synchrony + phytoplankton_synchrony + fungi_synchrony + bacteria_synchrony, dat)
  )
}

fit_model_2 <- function(dat) {
  piecewiseSEM::psem(
    lm(benthic_invertebrate_richness ~ fish_richness, dat),
    lm(insect_richness ~ benthic_invertebrate_richness, dat),
    lm(zooplankton_richness ~ insect_richness, dat),
    lm(fungi_richness ~ benthic_invertebrate_richness, dat),
    lm(fish_synchrony ~ fish_richness, dat),
    lm(benthic_invertebrate_synchrony ~ benthic_invertebrate_richness, dat),
    lm(insect_synchrony ~ insect_richness, dat),
    lm(zooplankton_synchrony ~ zooplankton_richness + benthic_invertebrate_richness, dat),
    lm(phytoplankton_synchrony ~ phytoplankton_richness + fish_synchrony + benthic_invertebrate_synchrony, dat),
    lm(fungi_synchrony ~ fungi_richness + benthic_invertebrate_synchrony, dat),
    lm(bacteria_synchrony ~ bacteria_richness + benthic_invertebrate_synchrony, dat),
    lm(fish_stability ~ fish_synchrony, dat),
    lm(benthic_invertebrate_stability ~ benthic_invertebrate_synchrony, dat),
    lm(insect_stability ~ insect_synchrony, dat),
    lm(zooplankton_stability ~ zooplankton_synchrony, dat),
    lm(phytoplankton_stability ~ phytoplankton_synchrony, dat),
    lm(fungi_stability ~ fungi_synchrony, dat),
    lm(bacteria_stability ~ bacteria_synchrony, dat),
    lm(multigroups_stability ~ fish_richness + benthic_invertebrate_richness + insect_richness + zooplankton_richness + phytoplankton_richness + fungi_richness + bacteria_richness, dat)
  )
}

fit_model_3 <- function(dat) {
  piecewiseSEM::psem(
    lm(benthic_invertebrate_diversity ~ fish_diversity, dat),
    lm(insect_diversity ~ benthic_invertebrate_diversity, dat),
    lm(zooplankton_diversity ~ insect_diversity, dat),
    lm(fungi_diversity ~ benthic_invertebrate_diversity, dat),
    lm(fish_synchrony ~ fish_diversity, dat),
    lm(benthic_invertebrate_synchrony ~ benthic_invertebrate_diversity, dat),
    lm(insect_synchrony ~ insect_diversity, dat),
    lm(zooplankton_synchrony ~ zooplankton_diversity + benthic_invertebrate_diversity, dat),
    lm(phytoplankton_synchrony ~ phytoplankton_diversity + fish_synchrony + benthic_invertebrate_synchrony, dat),
    lm(fungi_synchrony ~ fungi_diversity + benthic_invertebrate_synchrony, dat),
    lm(bacteria_synchrony ~ bacteria_diversity + benthic_invertebrate_synchrony, dat),
    lm(fish_stability ~ fish_synchrony, dat),
    lm(benthic_invertebrate_stability ~ benthic_invertebrate_synchrony, dat),
    lm(insect_stability ~ insect_synchrony, dat),
    lm(zooplankton_stability ~ zooplankton_synchrony, dat),
    lm(phytoplankton_stability ~ phytoplankton_synchrony, dat),
    lm(fungi_stability ~ fungi_synchrony, dat),
    lm(bacteria_stability ~ bacteria_synchrony, dat),
    lm(multigroups_stability ~ fish_synchrony + benthic_invertebrate_synchrony + insect_synchrony + zooplankton_synchrony + phytoplankton_synchrony + fungi_synchrony + bacteria_synchrony, dat)
  )
}

fit_model_4 <- function(dat) {
  piecewiseSEM::psem(
    lm(benthic_invertebrate_diversity ~ fish_diversity, dat),
    lm(insect_diversity ~ benthic_invertebrate_diversity, dat),
    lm(zooplankton_diversity ~ insect_diversity, dat),
    lm(fungi_diversity ~ benthic_invertebrate_diversity, dat),
    lm(fish_synchrony ~ fish_diversity, dat),
    lm(benthic_invertebrate_synchrony ~ benthic_invertebrate_diversity, dat),
    lm(insect_synchrony ~ insect_diversity, dat),
    lm(zooplankton_synchrony ~ zooplankton_diversity + benthic_invertebrate_diversity, dat),
    lm(phytoplankton_synchrony ~ phytoplankton_diversity + fish_synchrony + benthic_invertebrate_synchrony, dat),
    lm(fungi_synchrony ~ fungi_diversity + benthic_invertebrate_synchrony, dat),
    lm(bacteria_synchrony ~ bacteria_diversity + benthic_invertebrate_synchrony, dat),
    lm(fish_stability ~ fish_synchrony, dat),
    lm(benthic_invertebrate_stability ~ benthic_invertebrate_synchrony, dat),
    lm(insect_stability ~ insect_synchrony, dat),
    lm(zooplankton_stability ~ zooplankton_synchrony, dat),
    lm(phytoplankton_stability ~ phytoplankton_synchrony, dat),
    lm(fungi_stability ~ fungi_synchrony, dat),
    lm(bacteria_stability ~ bacteria_synchrony, dat),
    lm(multigroups_stability ~ fish_diversity + benthic_invertebrate_diversity + insect_diversity + zooplankton_diversity + phytoplankton_diversity + fungi_diversity + bacteria_diversity, dat)
  )
}

fit_model_5 <- function(dat) {
  piecewiseSEM::psem(
    lm(mean_path_length ~ fish_richness + zooplankton_richness + benthic_invertebrate_richness + phytoplankton_richness, dat),
    lm(modularity ~ mean_path_length + fish_richness + zooplankton_richness + benthic_invertebrate_richness + phytoplankton_richness, dat),
    lm(vulnerability ~ mean_path_length + fish_richness + zooplankton_richness + benthic_invertebrate_richness + phytoplankton_richness, dat),
    lm(robustness ~ mean_path_length + modularity + fish_richness + zooplankton_richness + benthic_invertebrate_richness + phytoplankton_richness + vulnerability, dat)
  )
}

fit_model_6 <- function(dat) {
  piecewiseSEM::psem(
    lm(mean_path_length ~ fish_diversity + zooplankton_diversity + benthic_invertebrate_diversity + phytoplankton_diversity, dat),
    lm(modularity ~ mean_path_length + fish_diversity + zooplankton_diversity + benthic_invertebrate_diversity + phytoplankton_diversity, dat),
    lm(vulnerability ~ mean_path_length + fish_diversity + zooplankton_diversity + benthic_invertebrate_diversity + phytoplankton_diversity, dat),
    lm(robustness ~ mean_path_length + modularity + fish_diversity + zooplankton_diversity + benthic_invertebrate_diversity + phytoplankton_diversity + vulnerability, dat)
  )
}

psem_bootstrap_analysis <- function(label, fit_fun, dat, R = 1000, seed = 123) {
  cat("\n\n")
  cat("  ", label, "\n")
  
  standardize_df <- function(df) {
    df2 <- df
    num_cols <- vapply(df, is.numeric, logical(1))
    df2[num_cols] <- lapply(df[num_cols], function(x) {
      s <- stats::sd(x, na.rm = TRUE)
      if (is.na(s) || s == 0) {
        return(x - mean(x, na.rm = TRUE))
      }
      (x - mean(x, na.rm = TRUE)) / s
    })
    df2
  }
  
  model <- fit_fun(dat)
  
  cat("\n", "pSEM Result Summary", "\n")
  print(summary(model))
  cat("\n")
  
  cf0 <- piecewiseSEM::coefs(model)
  cf0 <- cf0[cf0$Predictor != "(Intercept)", , drop = FALSE]
  path_names <- paste0(cf0$Response, " ~ ", cf0$Predictor)
  n_paths <- length(path_names)
  dat_std <- standardize_df(dat)
  model_std <- tryCatch(
    suppressWarnings(fit_fun(dat_std)),
    error = function(e) NULL
  )
  
  std_est <- rep(NA_real_, n_paths)
  if (!is.null(model_std)) {
    cf0_std <- tryCatch(
      suppressWarnings(piecewiseSEM::coefs(model_std)),
      error = function(e) NULL
    )
    if (!is.null(cf0_std)) {
      cf0_std <- cf0_std[cf0_std$Predictor != "(Intercept)", , drop = FALSE]
      tmp <- cf0_std$Estimate
      names(tmp) <- paste0(cf0_std$Response, " ~ ", cf0_std$Predictor)
      std_est <- tmp[path_names]
    }
  }
  
  boot_stat <- function(d, idx) {
    dd <- d[idx, , drop = FALSE]
    
    m <- tryCatch(
      suppressWarnings(fit_fun(dd)),
      error = function(e) NULL
    )
    
    if (is.null(m)) {
      orig_vals <- rep(NA_real_, n_paths)
    } else {
      cf <- tryCatch(
        suppressWarnings(piecewiseSEM::coefs(m)),
        error = function(e) NULL
      )
      if (is.null(cf)) {
        orig_vals <- rep(NA_real_, n_paths)
      } else {
        cf <- cf[cf$Predictor != "(Intercept)", , drop = FALSE]
        out <- cf$Estimate
        names(out) <- paste0(cf$Response, " ~ ", cf$Predictor)
        orig_vals <- out[path_names]
      }
    }
    
    dd_std <- standardize_df(dd)
    m_std <- tryCatch(
      suppressWarnings(fit_fun(dd_std)),
      error = function(e) NULL
    )
    
    if (is.null(m_std)) {
      std_vals <- rep(NA_real_, n_paths)
    } else {
      cf_std <- tryCatch(
        suppressWarnings(piecewiseSEM::coefs(m_std)),
        error = function(e) NULL
      )
      if (is.null(cf_std)) {
        std_vals <- rep(NA_real_, n_paths)
      } else {
        cf_std <- cf_std[cf_std$Predictor != "(Intercept)", , drop = FALSE]
        out_std <- cf_std$Estimate
        names(out_std) <- paste0(cf_std$Response, " ~ ", cf_std$Predictor)
        std_vals <- out_std[path_names]
      }
    }
    
    c(orig_vals, std_vals)
  }
  
  set.seed(seed)
  boot_res <- boot(data = dat, statistic = boot_stat, R = R)
  
  boot_t_orig <- boot_res$t[, seq_len(n_paths), drop = FALSE]
  boot_t_std <- boot_res$t[, (n_paths + 1):(2 * n_paths), drop = FALSE]
  
  orig_sign <- sign(cf0$Estimate)
  
  prop_same_dir <- sapply(seq_len(n_paths), function(i) {
    v <- boot_t_orig[, i]
    v <- v[!is.na(v)]
    if (length(v) == 0L) {
      return(NA_real_)
    }
    mean(sign(v) == orig_sign[i])
  })
  
  boot_summary <- data.frame(
    Path = path_names,
    Estimate = cf0$Estimate,
    Boot_SE = apply(boot_t_orig, 2, sd, na.rm = TRUE),
    Bias = apply(boot_t_orig, 2, mean, na.rm = TRUE) - cf0$Estimate,
    Boot_Median = apply(boot_t_orig, 2, median, na.rm = TRUE),
    CI_lower = apply(boot_t_orig, 2, quantile, probs = 0.025, na.rm = TRUE),
    CI_upper = apply(boot_t_orig, 2, quantile, probs = 0.975, na.rm = TRUE),
    N_ok = apply(!is.na(boot_t_orig), 2, sum),
    Prop_same_dir = prop_same_dir,
    stringsAsFactors = FALSE
  )
  
  boot_summary$p_boot <- 2 * pmin(
    colMeans(boot_t_orig <= 0, na.rm = TRUE),
    colMeans(boot_t_orig >= 0, na.rm = TRUE)
  )
  
  boot_summary$Signif <- ifelse(
    boot_summary$CI_lower > 0 | boot_summary$CI_upper < 0,
    "*",
    ""
  )
  
  bca_lo <- bca_hi <- rep(NA_real_, n_paths)
  for (i in seq_len(n_paths)) {
    ci <- tryCatch(
      boot.ci(boot_res, type = "bca", index = i),
      error = function(e) NULL
    )
    if (!is.null(ci) && !is.null(ci$bca)) {
      bca_lo[i] <- ci$bca[4]
      bca_hi[i] <- ci$bca[5]
    }
  }
  
  boot_summary$BCa_lower <- bca_lo
  boot_summary$BCa_upper <- bca_hi
  boot_summary$Signif_BCa <- ifelse(
    !is.na(bca_lo) & (bca_lo > 0 | bca_hi < 0),
    "*",
    ""
  )
  
  boot_summary$Std_Estimate <- std_est
  boot_summary$Std_Boot_SE <- apply(boot_t_std, 2, sd, na.rm = TRUE)
  boot_summary$Std_Bias <- apply(boot_t_std, 2, mean, na.rm = TRUE) - std_est
  boot_summary$Std_CI_lower <- apply(
    boot_t_std, 2, quantile,
    probs = 0.025, na.rm = TRUE
  )
  boot_summary$Std_CI_upper <- apply(
    boot_t_std, 2, quantile,
    probs = 0.975, na.rm = TRUE
  )
  boot_summary$Std_p_boot <- 2 * pmin(
    colMeans(boot_t_std <= 0, na.rm = TRUE),
    colMeans(boot_t_std >= 0, na.rm = TRUE)
  )
  boot_summary$Std_Signif <- ifelse(
    boot_summary$Std_CI_lower > 0 | boot_summary$Std_CI_upper < 0,
    "*",
    ""
  )
  
  std_bca_lo <- std_bca_hi <- rep(NA_real_, n_paths)
  for (i in seq_len(n_paths)) {
    ci <- tryCatch(
      boot.ci(boot_res, type = "bca", index = n_paths + i),
      error = function(e) NULL
    )
    if (!is.null(ci) && !is.null(ci$bca)) {
      std_bca_lo[i] <- ci$bca[4]
      std_bca_hi[i] <- ci$bca[5]
    }
  }
  
  boot_summary$Std_BCa_lower <- std_bca_lo
  boot_summary$Std_BCa_upper <- std_bca_hi
  boot_summary$Std_Signif_BCa <- ifelse(
    !is.na(std_bca_lo) & (std_bca_lo > 0 | std_bca_hi < 0),
    "*",
    ""
  )
  
  cat("Bootstrap Result Summary", "\n")
  cat("\n")
  cat("  Prop_same_dir = the proportion of bootstrap estimates retaining\n")
  cat("                  the original effect direction\n")
  cat("\n")
  print(boot_summary, row.names = FALSE, digits = 4)
  cat("\n")
  
  invisible(list(model = model, boot = boot_res, summary = boot_summary))
}

export_psem_results <- function(res_list, labels, out_file) {
  wb <- createWorkbook()
  hdr_names <- c("Coefficients", "Individual R-squared",
                 "Global goodness-of-fit", "Tests of directed separation")
  for (i in seq_along(res_list)) {
    res <- res_list[[i]]
    label <- labels[i]
    sm <- summary(res$model)
    lines <- capture.output(print(sm))
    
    hdr_idx <- sapply(hdr_names, function(h) {
      idx <- grep(paste0("^", h, ":?"), lines)
      if (length(idx) == 0) NA_integer_ else idx[1]
    })
    
    get_lines <- function(h) {
      i_start <- hdr_idx[h]
      if (is.na(i_start)) return(character(0))
      nexts <- hdr_idx[!is.na(hdr_idx) & hdr_idx > i_start]
      i_end <- if (length(nexts) > 0) min(nexts) - 1 else length(lines)
      txt <- lines[(i_start + 1):i_end]
      while (length(txt) > 0 && !nzchar(trimws(txt[1]))) txt <- txt[-1]
      while (length(txt) > 0 && !nzchar(trimws(txt[length(txt)]))) txt <- txt[-length(txt)]
      txt
    }
    
    sh <- paste0(label, " psem")
    addWorksheet(wb, sh)
    row <- 1
    
    writeData(wb, sh, "Coefficients:", startRow = row, colNames = FALSE)
    row <- row + 1
    coef_df <- as.data.frame(sm$coefficients)
    rownames(coef_df) <- NULL
    writeData(wb, sh, coef_df, startRow = row, colNames = TRUE)
    row <- row + nrow(coef_df) + 2
    
    for (h in hdr_names[-1]) {
      writeData(wb, sh, paste0(h, ":"), startRow = row, colNames = FALSE)
      row <- row + 1
      txt <- get_lines(h)
      if (length(txt) == 0) {
        writeData(wb, sh, "(empty)", startRow = row, colNames = FALSE)
        row <- row + 2
      } else {
        for (ln in txt) {
          writeData(wb, sh, ln, startRow = row, colNames = FALSE)
          row <- row + 1
        }
        row <- row + 1
      }
    }
    
    sh2 <- paste0(label, " bootstrap")
    addWorksheet(wb, sh2)
    writeData(wb, sh2, res$summary, rowNames = FALSE)
  }
  saveWorkbook(wb, out_file, overwrite = TRUE)
}

R_boot <- 1000
res_1 <- psem_bootstrap_analysis("Model 1", fit_model_1, data_average, R = R_boot)
res_2 <- psem_bootstrap_analysis("Model 2", fit_model_2, data_average, R = R_boot)
res_3 <- psem_bootstrap_analysis("Model 3", fit_model_3, data_average, R = R_boot)
res_4 <- psem_bootstrap_analysis("Model 4", fit_model_4, data_average, R = R_boot)
res_5 <- psem_bootstrap_analysis("Model 5", fit_model_5, data_original, R = R_boot)
res_6 <- psem_bootstrap_analysis("Model 6", fit_model_6, data_original, R = R_boot)

export_psem_results(
  list(res_1, res_2, res_3, res_4, res_5, res_6),
  c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6"),
  "Result_Tables/pSEM Results.xlsx"
)
