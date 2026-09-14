
################################################################################
# GWAASB-WHEAT — REPRODUCIBLE GY-ONLY ANALYSIS PIPELINE
# Grain-yield-only reference implementation
#
# Included analyses:
#   1) GBLUP
#   2) GBLUP-GxE
#   3) CV1 and CV2
#   4) 100 independent repetitions, 5 folds each
#   5) Same folds for both models
#   6) CV split BEFORE phenotype adjustment
#   7) Accuracy calculated within environment
#   8) Mean, SD, SE, 95% CI, range
#   9) Paired Delta-r = GBLUP-GxE - GBLUP
#  10) Full-data genomic GE -> GWAASB
#  11) GWAASB vs phenotypic WAASB, GE-row SD, and Wricke
#
# Parents G101/G102 are excluded.
################################################################################


################################################################################
# 0. PACKAGES
################################################################################

pkgs <- c(
  "tidyverse",
  "writexl",
  "rrBLUP",
  "BGLR",
  "emmeans",
  "metan"
)

missing <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  install.packages(missing, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(writexl)
  library(rrBLUP)
  library(BGLR)
  library(emmeans)
  library(metan)
})

options(stringsAsFactors = FALSE)

################################################################################
# 1. SETTINGS
################################################################################

TRAIT <- "GY"

POP_FILES <- list(
  H3 = list(
    pheno = "phenotype_data_h3.csv",
    geno  = "SNP_clean_fixed_h3.csv",
    parents = c("G101", "G102")
  ),
  H4 = list(
    pheno = "phenotype_data_h4.csv",
    geno  = "SNP_clean_fixed_h4.csv",
    parents = c("G101", "G102")
  )
)

EXPECTED_DH <- 100

N_FOLDS <- 5
N_REPS  <- 100

BASE_SEED <- 20260824

# Full-data model for final GWAASB
FULL_NITER  <- 12000
FULL_BURNIN <- 4000
FULL_THIN   <- 5

# Repeated CV model
# This is much lighter than the prior script but still suitable for a final run.
CV_NITER  <- 5000
CV_BURNIN <- 1500
CV_THIN   <- 5

OUTDIR <- "results"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


################################################################################
# 2. GENERAL HELPERS
################################################################################

find_file_ci <- function(fname) {
  if (file.exists(fname)) return(fname)

  x <- list.files(".", full.names = FALSE)
  hit <- x[tolower(x) == tolower(fname)]

  if (length(hit) == 1) return(hit)

  stop("File not found: ", fname)
}


find_col <- function(dat, target) {
  i <- which(tolower(names(dat)) == tolower(target))

  if (length(i) != 1) {
    stop(
      "Column ", target, " not found uniquely.\n",
      paste(names(dat), collapse = ", ")
    )
  }

  names(dat)[i]
}


safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 3) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)

  cor(x[ok], y[ok], method = method)
}


rmse_fun <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((obs[ok] - pred[ok])^2))
}


mae_fun <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (!any(ok)) return(NA_real_)
  mean(abs(obs[ok] - pred[ok]))
}


ci95 <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) < 2) return(c(NA_real_, NA_real_))

  m <- mean(x)
  se <- sd(x) / sqrt(length(x))

  m + c(-1, 1) * qt(0.975, df = length(x) - 1) * se
}


################################################################################
# 3. READ PHENOTYPES
################################################################################

read_pheno <- function(path) {

  dat <- read.csv(path, check.names = FALSE)

  names(dat)[names(dat) == find_col(dat, "ENV")] <- "ENV"
  names(dat)[names(dat) == find_col(dat, "GEN")] <- "GEN"
  names(dat)[names(dat) == find_col(dat, "REP")] <- "REP"
  names(dat)[names(dat) == find_col(dat, TRAIT)] <- "GY"

  # Hard GY-only boundary: discard every column outside the approved schema
  # immediately after import so downstream objects and exports cannot retain
  # additional phenotypic variables.
  dat <- dat[, c("ENV", "GEN", "REP", "GY"), drop = FALSE]

  dat %>%
    mutate(
      ENV = factor(ENV),
      GEN = factor(GEN),
      REP = factor(REP),
      GY = as.numeric(GY)
    )
}


################################################################################
# 4. READ SNPs
################################################################################

read_geno <- function(path) {

  x <- read.csv(path, check.names = FALSE)

  gcol <- names(x)[tolower(names(x)) == "gen"]

  if (length(gcol) != 1) {
    stop("SNP file must contain one GEN column.")
  }

  ids <- as.character(x[[gcol]])
  x[[gcol]] <- NULL

  M <- as.matrix(x)
  storage.mode(M) <- "numeric"
  rownames(M) <- ids

  vals <- sort(unique(as.numeric(M[is.finite(M)])))

  if (all(vals %in% c(0, 1, 2))) {
    M <- M - 1
  } else if (!all(vals %in% c(-1, 0, 1))) {
    stop(
      "Unexpected SNP coding: ",
      paste(vals, collapse = ", ")
    )
  }

  M
}


################################################################################
# 5. EXCLUDE PARENTS AND MATCH DATA
################################################################################

exclude_parents_match <- function(pheno, M, parents, pop) {

  pheno <- pheno %>%
    filter(!as.character(GEN) %in% parents) %>%
    droplevels()

  M <- M[
    !rownames(M) %in% parents,
    ,
    drop = FALSE
  ]

  common <- intersect(
    unique(as.character(pheno$GEN)),
    rownames(M)
  )

  pheno <- pheno %>%
    filter(as.character(GEN) %in% common) %>%
    droplevels()

  M <- M[
    common,
    ,
    drop = FALSE
  ]

  n <- n_distinct(pheno$GEN)

  if (n != EXPECTED_DH) {
    stop(
      pop,
      ": expected 100 DH lines after excluding G101/G102; found ",
      n
    )
  }

  message(pop, ": 100 DH lines retained.")

  list(pheno = pheno, M = M)
}


################################################################################
# 6. MARKER QC SUMMARY
################################################################################

marker_qc <- function(M, pop) {

  D <- M + 1

  miss <- colMeans(is.na(D))
  p <- colMeans(D, na.rm = TRUE) / 2
  maf <- pmin(p, 1 - p)
  mono <- apply(
    D,
    2,
    function(z) length(unique(z[is.finite(z)])) <= 1
  )

  tibble(
    Population = pop,
    N_genotypes = nrow(M),
    N_markers = ncol(M),
    Monomorphic_markers = sum(mono),
    Maximum_missing_rate = max(miss, na.rm = TRUE),
    Median_missing_rate = median(miss, na.rm = TRUE),
    Minimum_MAF = min(maf, na.rm = TRUE),
    Median_MAF = median(maf, na.rm = TRUE)
  )
}


################################################################################
# 7. FULL-DATA ENVIRONMENT BLUEs
################################################################################

stage1_full <- function(pheno) {

  split(pheno, pheno$ENV) %>%
    imap_dfr(function(d, ee) {

      d <- droplevels(d)

      fit <- lm(
        GY ~ REP + GEN,
        data = d
      )

      em <- as.data.frame(
        emmeans(fit, ~ GEN)
      )

      tibble(
        ENV = ee,
        GEN = as.character(em$GEN),
        BLUE = em$emmean,
        BLUE_SE = em$SE
      )
    }) %>%
    arrange(ENV, GEN)
}


################################################################################
# 8. LEAKAGE-CONTROLLED NESTED STAGE 1
################################################################################

stage1_nested <- function(
    pheno,
    train_cells,
    test_cells
) {

  train_cells <- train_cells %>%
    mutate(
      GEN = as.character(GEN),
      ENV = as.character(ENV)
    )

  test_cells <- test_cells %>%
    mutate(
      GEN = as.character(GEN),
      ENV = as.character(ENV)
    )

  ph <- pheno %>%
    mutate(
      GEN_chr = as.character(GEN),
      ENV_chr = as.character(ENV),
      REP_chr = as.character(REP)
    )

  tr_plot <- ph %>%
    inner_join(
      train_cells,
      by = c(
        "GEN_chr" = "GEN",
        "ENV_chr" = "ENV"
      )
    )

  te_plot <- ph %>%
    inner_join(
      test_cells,
      by = c(
        "GEN_chr" = "GEN",
        "ENV_chr" = "ENV"
      )
    )

  envs <- sort(unique(as.character(pheno$ENV)))

  tr_out <- list()
  te_out <- list()

  for (ee in envs) {

    tr <- tr_plot %>%
      filter(ENV_chr == ee) %>%
      droplevels()

    te <- te_plot %>%
      filter(ENV_chr == ee) %>%
      droplevels()

    if (nrow(tr) == 0) next

    # Training only
    fit <- lm(
      GY ~ REP + GEN,
      data = tr
    )

    emg <- as.data.frame(
      emmeans(fit, ~ GEN)
    )

    tr_out[[ee]] <- tibble(
      ENV = ee,
      GEN = as.character(emg$GEN),
      BLUE = emg$emmean,
      BLUE_SE = emg$SE
    )

    if (nrow(te) > 0) {

      # Training-only REP adjustment
      emr <- as.data.frame(
        emmeans(fit, ~ REP)
      )

      radj <- tibble(
        REP_chr = as.character(emr$REP),
        REP_mean = emr$emmean
      )

      grand_rep <- mean(
        radj$REP_mean,
        na.rm = TRUE
      )

      radj <- radj %>%
        mutate(
          REP_dev = REP_mean - grand_rep
        )

      te2 <- te %>%
        left_join(
          radj,
          by = "REP_chr"
        )

      if (anyNA(te2$REP_dev)) {
        stop(
          "Validation REP level absent from training in ENV ",
          ee
        )
      }

      te_out[[ee]] <- te2 %>%
        mutate(
          GY_adjusted =
            GY - REP_dev
        ) %>%
        group_by(
          GEN_chr,
          ENV_chr
        ) %>%
        summarise(
          Observed =
            mean(
              GY_adjusted,
              na.rm = TRUE
            ),
          .groups = "drop"
        ) %>%
        transmute(
          GEN = GEN_chr,
          ENV = ENV_chr,
          Observed
        )
    }
  }

  list(
    train_BLUE =
      bind_rows(tr_out),

    test_observed =
      bind_rows(te_out)
  )
}


################################################################################
# 9. BUILD G + MULTI-ENVIRONMENT KERNELS
################################################################################

prepare_design <- function(M, pheno) {

  gens <- sort(
    intersect(
      rownames(M),
      unique(as.character(pheno$GEN))
    )
  )

  envs <- sort(
    unique(as.character(pheno$ENV))
  )

  M <- M[
    gens,
    ,
    drop = FALSE
  ]

  G <- rrBLUP::A.mat(
    M,
    impute.method = "mean"
  )

  G <- G[
    gens,
    gens,
    drop = FALSE
  ]

  grid <- expand.grid(
    GEN = gens,
    ENV = envs,
    stringsAsFactors = FALSE
  ) %>%
    arrange(ENV, GEN)

  grid$GEN <- factor(
    grid$GEN,
    levels = gens
  )

  grid$ENV <- factor(
    grid$ENV,
    levels = envs
  )

  Zg <- model.matrix(
    ~ 0 + GEN,
    data = grid
  )

  colnames(Zg) <- sub(
    "^GEN",
    "",
    colnames(Zg)
  )

  Zg <- Zg[
    ,
    gens,
    drop = FALSE
  ]

  Ze <- model.matrix(
    ~ 0 + ENV,
    data = grid
  )

  KG <- Zg %*% G %*% t(Zg)

  KE <- Ze %*% t(Ze)

  KGE <- KG * KE

  ETA_G <- list(
    ENV = list(
      X = Ze,
      model = "FIXED"
    ),
    G = list(
      K = KG,
      model = "RKHS"
    )
  )

  ETA_GE <- list(
    ENV = list(
      X = Ze,
      model = "FIXED"
    ),
    G = list(
      K = KG,
      model = "RKHS"
    ),
    GE = list(
      K = KGE,
      model = "RKHS"
    )
  )

  list(
    M = M,
    G = G,
    gens = gens,
    envs = envs,
    grid = grid,
    ETA_G = ETA_G,
    ETA_GE = ETA_GE
  )
}


################################################################################
# 10. MAP TRAINING BLUEs TO COMPLETE GRID
################################################################################

make_response <- function(grid, train_blue) {

  tab <- grid %>%
    transmute(
      GEN = as.character(GEN),
      ENV = as.character(ENV)
    ) %>%
    left_join(
      train_blue,
      by = c("GEN", "ENV")
    )

  y <- tab$BLUE

  se <- tab$BLUE_SE

  fallback <- median(
    se[
      is.finite(se) &
        se > 0
    ],
    na.rm = TRUE
  )

  if (!is.finite(fallback)) {
    fallback <- 1
  }

  se[
    !is.finite(se) |
      se <= 0
  ] <- fallback

  weights <- fallback / se

  weights[is.na(y)] <- 1

  list(
    y = y,
    weights = weights
  )
}


################################################################################
# 11. BGLR FITTER
################################################################################

fit_bglr <- function(
    y,
    ETA,
    weights,
    nIter,
    burnIn,
    thin,
    seed,
    tag
) {

  set.seed(seed)

  BGLR(
    y = y,
    ETA = ETA,
    response_type = "gaussian",
    weights = weights,
    nIter = nIter,
    burnIn = burnIn,
    thin = thin,
    saveAt = file.path(
      tempdir(),
      paste0(
        tag,
        "_",
        seed,
        "_"
      )
    ),
    rmExistingFiles = TRUE,
    verbose = FALSE
  )
}


################################################################################
# 12. CV1 PARTITIONS
################################################################################

make_cv1 <- function(
    gens,
    seed
) {

  set.seed(seed)

  g <- sample(gens)

  tibble(
    GEN = g,
    Fold = rep(
      seq_len(N_FOLDS),
      length.out = length(g)
    )
  )
}


################################################################################
# 13. CV2 PARTITIONS
################################################################################

make_cv2 <- function(
    gens,
    envs,
    seed
) {

  set.seed(seed)

  map_dfr(
    gens,
    function(g) {

      e <- sample(envs)

      tibble(
        GEN = g,
        ENV = e,
        Fold = rep(
          seq_len(N_FOLDS),
          length.out = length(e)
        )
      )
    }
  )
}


################################################################################
# 14. RUN ONE FOLD
################################################################################

run_one_fold <- function(
    pheno,
    design,
    cv_type,
    fold,
    rep_id,
    cv1_map,
    cv2_map,
    pop
) {

  all_cells <- expand.grid(
    GEN = design$gens,
    ENV = design$envs,
    stringsAsFactors = FALSE
  )

  if (cv_type == "CV1") {

    test_gens <- cv1_map %>%
      filter(Fold == fold) %>%
      pull(GEN)

    test_cells <- all_cells %>%
      filter(GEN %in% test_gens)

  } else {

    test_cells <- cv2_map %>%
      filter(Fold == fold) %>%
      select(GEN, ENV)
  }

  train_cells <- anti_join(
    all_cells,
    test_cells,
    by = c("GEN", "ENV")
  )

  s1 <- stage1_nested(
    pheno = pheno,
    train_cells = train_cells,
    test_cells = test_cells
  )

  resp <- make_response(
    design$grid,
    s1$train_BLUE
  )

  seed0 <- BASE_SEED +
    rep_id * 100000 +
    fold * 1000 +
    ifelse(
      cv_type == "CV2",
      50000000,
      0
    )

  # Model 1
  fit_G <- fit_bglr(
    y = resp$y,
    ETA = design$ETA_G,
    weights = resp$weights,
    nIter = CV_NITER,
    burnIn = CV_BURNIN,
    thin = CV_THIN,
    seed = seed0 + 1,
    tag = paste(
      pop,
      cv_type,
      rep_id,
      fold,
      "G",
      sep = "_"
    )
  )

  # Model 2
  fit_GE <- fit_bglr(
    y = resp$y,
    ETA = design$ETA_GE,
    weights = resp$weights,
    nIter = CV_NITER,
    burnIn = CV_BURNIN,
    thin = CV_THIN,
    seed = seed0 + 2,
    tag = paste(
      pop,
      cv_type,
      rep_id,
      fold,
      "GE",
      sep = "_"
    )
  )

  pred_grid <- design$grid %>%
    transmute(
      GEN = as.character(GEN),
      ENV = as.character(ENV),
      Pred_GBLUP = fit_G$yHat,
      Pred_GBLUP_GxE = fit_GE$yHat
    )

  pred <- s1$test_observed %>%
    left_join(
      pred_grid,
      by = c("GEN", "ENV")
    ) %>%
    mutate(
      POP = pop,
      CV = cv_type,
      Repetition = rep_id,
      Fold = fold,
      .before = 1
    )

  # Environment-specific prediction accuracy
  metrics <- pred %>%
    group_by(
      POP,
      CV,
      Repetition,
      Fold,
      ENV
    ) %>%
    summarise(
      N_validation =
        n_distinct(GEN),

      r_GBLUP =
        safe_cor(
          Observed,
          Pred_GBLUP
        ),

      r_GBLUP_GxE =
        safe_cor(
          Observed,
          Pred_GBLUP_GxE
        ),

      RMSE_GBLUP =
        rmse_fun(
          Observed,
          Pred_GBLUP
        ),

      RMSE_GBLUP_GxE =
        rmse_fun(
          Observed,
          Pred_GBLUP_GxE
        ),

      MAE_GBLUP =
        mae_fun(
          Observed,
          Pred_GBLUP
        ),

      MAE_GBLUP_GxE =
        mae_fun(
          Observed,
          Pred_GBLUP_GxE
        ),

      .groups = "drop"
    )

  list(
    predictions = pred,
    metrics = metrics
  )
}


################################################################################
# 15. REPEATED CV1 + CV2
################################################################################

run_cv <- function(
    pheno,
    design,
    pop
) {

  pred_list <- list()
  metric_list <- list()

  z <- 0

  for (r in seq_len(N_REPS)) {

    message(
      pop,
      " repetition ",
      r,
      "/",
      N_REPS
    )

    cv1 <- make_cv1(
      design$gens,
      BASE_SEED + r * 10 + 1
    )

    cv2 <- make_cv2(
      design$gens,
      design$envs,
      BASE_SEED + r * 10 + 2
    )

    for (cv_type in c("CV1", "CV2")) {

      for (f in seq_len(N_FOLDS)) {

        ans <- run_one_fold(
          pheno = pheno,
          design = design,
          cv_type = cv_type,
          fold = f,
          rep_id = r,
          cv1_map = cv1,
          cv2_map = cv2,
          pop = pop
        )

        z <- z + 1

        pred_list[[z]] <-
          ans$predictions

        metric_list[[z]] <-
          ans$metrics
      }
    }
  }

  pred <- bind_rows(
    pred_list
  )

  env_metrics <- bind_rows(
    metric_list
  )

  # Average first within fold across environments
  fold_metrics <- env_metrics %>%
    group_by(
      POP,
      CV,
      Repetition,      Fold
    ) %>%
    summarise(
      Accuracy_GBLUP =
        mean(
          r_GBLUP,
          na.rm = TRUE
        ),

      Accuracy_GBLUP_GxE =
        mean(
          r_GBLUP_GxE,
          na.rm = TRUE
        ),

      RMSE_GBLUP =
        mean(
          RMSE_GBLUP,
          na.rm = TRUE
        ),

      RMSE_GBLUP_GxE =
        mean(
          RMSE_GBLUP_GxE,
          na.rm = TRUE
        ),

      .groups = "drop"
    )

  # Then average over the 5 folds in each independent repetition
  rep_metrics <- fold_metrics %>%
    group_by(
      POP,
      CV,
      Repetition
    ) %>%
    summarise(
      Accuracy_GBLUP =
        mean(
          Accuracy_GBLUP,
          na.rm = TRUE
        ),

      Accuracy_GBLUP_GxE =
        mean(
          Accuracy_GBLUP_GxE,
          na.rm = TRUE
        ),

      Delta_r =
        Accuracy_GBLUP_GxE -
        Accuracy_GBLUP,

      RMSE_GBLUP =
        mean(
          RMSE_GBLUP,
          na.rm = TRUE
        ),

      RMSE_GBLUP_GxE =
        mean(
          RMSE_GBLUP_GxE,
          na.rm = TRUE
        ),

      .groups = "drop"
    )

  list(
    predictions = pred,
    env_metrics = env_metrics,
    fold_metrics = fold_metrics,
    rep_metrics = rep_metrics
  )
}


################################################################################
# 16. CV SUMMARY
################################################################################

summarize_metric <- function(
    dat,
    var,
    model
) {

  x <- dat[[var]]
  x <- x[is.finite(x)]

  ci <- ci95(x)

  tibble(
    Model = model,
    Mean = mean(x),
    SD = sd(x),
    SE = sd(x) / sqrt(length(x)),
    CI95_lower = ci[1],
    CI95_upper = ci[2],
    Min = min(x),
    Max = max(x),
    N_repetitions = length(x)
  )
}


summarize_cv <- function(
    cv,
    pop
) {

  acc <- cv$rep_metrics %>%
    group_split(CV) %>%
    map_dfr(function(d) {

      cvname <- unique(d$CV)

      bind_rows(
        summarize_metric(
          d,
          "Accuracy_GBLUP",
          "GBLUP"
        ),

        summarize_metric(
          d,
          "Accuracy_GBLUP_GxE",
          "GBLUP-GxE"
        )
      ) %>%
        mutate(
          CV = cvname,
          .before = 1
        )
    }) %>%
    mutate(
      POP = pop,
      .before = 1
    )

  delta <- cv$rep_metrics %>%
    group_by(
      POP,
      CV
    ) %>%
    summarise(
      Mean_Delta_r =
        mean(
          Delta_r,
          na.rm = TRUE
        ),

      SD_Delta_r =
        sd(
          Delta_r,
          na.rm = TRUE
        ),

      SE_Delta_r =
        SD_Delta_r /
        sqrt(
          sum(
            is.finite(Delta_r)
          )
        ),

      Proportion_GxE_better =
        mean(
          Delta_r > 0,
          na.rm = TRUE
        ),

      Median_Delta_r =
        median(
          Delta_r,
          na.rm = TRUE
        ),

      Min_Delta_r =
        min(
          Delta_r,
          na.rm = TRUE
        ),

      Max_Delta_r =
        max(
          Delta_r,
          na.rm = TRUE
        ),

      .groups = "drop"
    )

  env <- cv$env_metrics %>%
    group_by(
      POP,
      CV,
      ENV
    ) %>%
    summarise(
      GBLUP_mean_r =
        mean(
          r_GBLUP,
          na.rm = TRUE
        ),

      GBLUP_SD_r =
        sd(
          r_GBLUP,
          na.rm = TRUE
        ),

      GBLUP_GxE_mean_r =
        mean(
          r_GBLUP_GxE,
          na.rm = TRUE
        ),

      GBLUP_GxE_SD_r =
        sd(
          r_GBLUP_GxE,
          na.rm = TRUE
        ),

      .groups = "drop"
    )

  list(
    accuracy = acc,
    delta = delta,
    environment = env
  )
}


################################################################################
# 17. FULL-DATA GBLUP + GBLUP-GxE
################################################################################

fit_full_models <- function(
    full_blues,
    design,
    pop
) {

  resp <- make_response(
    design$grid,
    full_blues
  )

  fit_G <- fit_bglr(
    y = resp$y,
    ETA = design$ETA_G,
    weights = resp$weights,
    nIter = FULL_NITER,
    burnIn = FULL_BURNIN,
    thin = FULL_THIN,
    seed = BASE_SEED + 11,
    tag = paste0(pop, "_FULL_G")
  )

  fit_GE <- fit_bglr(
    y = resp$y,
    ETA = design$ETA_GE,
    weights = resp$weights,
    nIter = FULL_NITER,
    burnIn = FULL_BURNIN,
    thin = FULL_THIN,
    seed = BASE_SEED + 22,
    tag = paste0(pop, "_FULL_GE")
  )

  design$grid %>%
    transmute(
      GEN = as.character(GEN),
      ENV = as.character(ENV),

      Observed_BLUE =
        resp$y,

      Pred_GBLUP =
        fit_G$yHat,

      Pred_GBLUP_GxE =
        fit_GE$yHat,

      Genomic_main_effect =
        fit_GE$ETA$G$u,

      Predicted_GE_effect =
        fit_GE$ETA$GE$u,

      Predicted_genotypic_value =
        fit_GE$ETA$G$u +
        fit_GE$ETA$GE$u
    )
}


################################################################################
# 18. GWAASB
################################################################################

make_GE_matrix <- function(pred) {

  pred %>%
    select(
      GEN,
      ENV,
      Predicted_GE_effect
    ) %>%
    pivot_wider(
      names_from = ENV,
      values_from =
        Predicted_GE_effect
    ) %>%
    column_to_rownames(
      "GEN"
    ) %>%
    as.matrix()
}


double_center <- function(M) {

  M -
    rowMeans(M) -
    matrix(
      colMeans(M),
      nrow = nrow(M),
      ncol = ncol(M),
      byrow = TRUE
    ) +
    mean(M)
}


calc_gwaasb <- function(GE) {

  GEc <- double_center(GE)

  sv <- svd(GEc)

  keep <- which(
    sv$d >
      sqrt(.Machine$double.eps) *
      max(sv$d)
  )

  d <- sv$d[keep]
  U <- sv$u[
    ,
    keep,
    drop = FALSE
  ]

  scores <- sweep(
    U,
    2,
    d,
    "*"
  )

  weights <-
    d^2 /
    sum(d^2)

  gwa <- rowSums(
    sweep(
      abs(scores),
      2,
      weights,
      "*"
    )
  ) /
    sum(weights)

  index <- tibble(
    GEN =
      rownames(GE),

    GWAASB =
      as.numeric(gwa),

    GE_row_SD =
      apply(
        GEc,
        1,
        sd
      ),

    Wricke_from_predicted_GE =
      rowSums(
        GEc^2
      )
  )

  axes <- tibble(
    Axis =
      paste0(
        "GIPC",
        seq_along(d)
      ),

    Singular_value = d,

    Percent =
      100 * weights,

    Cumulative_percent =
      100 * cumsum(weights)
  )

  list(
    index = index,
    axes = axes
  )
}


################################################################################
# 19. CONVENTIONAL PHENOTYPIC WAASB
################################################################################

phenotypic_waasb <- function(pheno) {

  m <- metan::waasb(
    pheno,
    env = ENV,
    gen = GEN,
    rep = REP,
    resp = GY,
    verbose = FALSE
  )

  w <- metan::get_model_data(
    m,
    what = "WAASB",
    verbose = FALSE
  ) %>%
    as.data.frame()

  gen_col <- names(w)[
    tolower(names(w)) %in%
      c(
        "gen",
        "genotype",
        "geno"
      )
  ]

  value_col <- names(w)[
    tolower(names(w)) ==
      tolower(TRAIT)
  ]

  if (length(value_col) == 0) {
    value_col <- names(w)[
      tolower(names(w)) ==
        "waasb"
    ]
  }

  if (
    length(gen_col) != 1 ||
      length(value_col) != 1
  ) {
    stop(
      "Cannot extract WAASB. Columns: ",
      paste(
        names(w),
        collapse = ", "
      )
    )
  }

  tibble(
    GEN =
      as.character(
        w[[gen_col]]
      ),

    Phenotypic_WAASB =
      as.numeric(
        w[[value_col]]
      )
  )
}


################################################################################
# 20. STABILITY COMPARISON
################################################################################

compare_stability <- function(
    pred,
    gwa,
    phenwa,
    pop
) {

  mean_perf <- pred %>%
    group_by(
      GEN
    ) %>%
    summarise(
      Mean_predicted_genotypic_value =
        mean(
          Predicted_genotypic_value,
          na.rm = TRUE
        ),
      .groups = "drop"
    )

  tab <- gwa$index %>%
    left_join(
      phenwa,
      by = "GEN"
    ) %>%
    left_join(
      mean_perf,
      by = "GEN"
    )

  cors <- tibble(
    Population = pop,

    Comparison = c(
      "GWAASB vs Phenotypic WAASB",
      "GWAASB vs GE row SD",
      "GWAASB vs Wricke"
    ),

    Pearson_r = c(
      safe_cor(
        tab$GWAASB,
        tab$Phenotypic_WAASB
      ),

      safe_cor(
        tab$GWAASB,
        tab$GE_row_SD
      ),

      safe_cor(
        tab$GWAASB,
        tab$Wricke_from_predicted_GE
      )
    ),

    Spearman_rho = c(
      safe_cor(
        tab$GWAASB,
        tab$Phenotypic_WAASB,
        "spearman"
      ),

      safe_cor(
        tab$GWAASB,
        tab$GE_row_SD,
        "spearman"
      ),

      safe_cor(
        tab$GWAASB,
        tab$Wricke_from_predicted_GE,
        "spearman"
      )
    )
  )

  list(
    table = tab,
    correlations = cors
  )
}


################################################################################
# 21. RUN ONE POPULATION
################################################################################

run_population <- function(
    pop,
    files
) {

  message(
    "\n==================== ",
    pop,
    " ====================\n"
  )

  pheno <- read_pheno(
    find_file_ci(
      files$pheno
    )
  )

  M <- read_geno(
    find_file_ci(
      files$geno
    )
  )

  z <- exclude_parents_match(
    pheno,
    M,
    files$parents,
    pop
  )

  pheno <- z$pheno
  M <- z$M

  qc <- marker_qc(
    M,
    pop
  )

  design <- prepare_design(
    M,
    pheno
  )

  write.csv(
    design$G,
    file.path(
      OUTDIR,
      paste0(
        pop,
        "_G_matrix.csv"
      )
    )
  )

  message(
    pop,
    ": full-data BLUEs..."
  )

  full_blues <- stage1_full(
    pheno
  )

  message(
    pop,
    ": full genomic models..."
  )

  full_pred <- fit_full_models(
    full_blues,
    design,
    pop
  )

  message(
    pop,
    ": GWAASB..."
  )

  GE <- make_GE_matrix(
    full_pred
  )

  gwa <- calc_gwaasb(
    GE
  )

  message(
    pop,
    ": phenotypic WAASB..."
  )

  pwa <- phenotypic_waasb(
    pheno
  )

  stab <- compare_stability(
    full_pred,
    gwa,
    pwa,
    pop
  )

  message(
    pop,
    ": CV1/CV2 x 100 repetitions..."
  )

  cv <- run_cv(
    pheno,
    design,
    pop
  )

  cvsum <- summarize_cv(
    cv,
    pop
  )

  # Useful plot: repeated accuracy distributions
  pp <- cv$rep_metrics %>%
    select(
      POP,
      CV,
      Repetition,
      Accuracy_GBLUP,
      Accuracy_GBLUP_GxE
    ) %>%
    pivot_longer(
      cols =
        starts_with(
          "Accuracy_"
        ),
      names_to =
        "Model",
      values_to =
        "Accuracy"
    ) %>%
    mutate(
      Model =
        recode(
          Model,
          "Accuracy_GBLUP" =
            "GBLUP",
          "Accuracy_GBLUP_GxE" =
            "GBLUP-GxE"
        )
    )

  p <- ggplot(
    pp,
    aes(
      Model,
      Accuracy
    )
  ) +
    geom_boxplot() +
    facet_wrap(
      ~ CV
    ) +
    theme_bw() +
    labs(
      title =
        paste0(
          pop,
          " repeated genomic prediction"
        ),
      x = NULL,
      y =
        "Within-environment prediction accuracy"
    )

  ggsave(
    filename =
      file.path(
        OUTDIR,
        paste0(
          pop,
          "_CV_accuracy.png"
        )
      ),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )

  sheets <- list(
    Marker_QC =
      qc,

    Stage1_BLUEs =
      full_blues,

    Full_predictions =
      full_pred,

    GWAASB =
      gwa$index,

    GWAASB_axes =
      gwa$axes,

    Phenotypic_WAASB =
      pwa,

    Stability_comparison =
      stab$table,

    Stability_correlations =
      stab$correlations,

    CV_accuracy_summary =
      cvsum$accuracy,

    CV_delta =
      cvsum$delta,

    CV_environment_accuracy =
      cvsum$environment,

    CV_repetition_metrics =
      cv$rep_metrics
  )

  write_xlsx(
    sheets,
    file.path(
      OUTDIR,
      paste0(
        pop,
        "_GY_GWAASB_RESULTS.xlsx"
      )
    )
  )

  # Large table saved as CSV
  write.csv(
    cv$predictions,
    file.path(
      OUTDIR,
      paste0(
        pop,
        "_CV_validation_predictions.csv"
      )
    ),
    row.names = FALSE
  )

  list(
    qc = qc,
    design = design,
    full_blues = full_blues,
    full_pred = full_pred,
    gwa = gwa,
    pwa = pwa,
    stability = stab,
    cv = cv,
    cvsum = cvsum
  )
}


################################################################################
# 22. RUN BOTH POPULATIONS
################################################################################

results <- list()

for (pop in names(POP_FILES)) {

  results[[pop]] <- run_population(
    pop,
    POP_FILES[[pop]]
  )
}


################################################################################
# 23. COMBINED TABLES
################################################################################

combined_acc <- map_dfr(
  results,
  ~ .x$cvsum$accuracy
)

combined_delta <- map_dfr(
  results,
  ~ .x$cvsum$delta
)

combined_env <- map_dfr(
  results,
  ~ .x$cvsum$environment
)

combined_stab <- map_dfr(
  results,
  ~ .x$stability$correlations
)

write_xlsx(
  list(
    Prediction_accuracy =
      combined_acc,

    Delta_GxE_minus_GBLUP =
      combined_delta,

    Environment_accuracy =
      combined_env,

    Stability_correlations =
      combined_stab
  ),
  file.path(
    OUTDIR,
    "COMBINED_GY_GWAASB_RESULTS.xlsx"
  )
)


################################################################################
# 24. SESSION INFO
################################################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)

message(
  "\nDONE. Results are in: ",
  normalizePath(OUTDIR),
  "\n"
)