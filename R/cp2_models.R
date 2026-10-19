# Shared feature construction, fold-local preprocessing, and model functions.
# This file defines functions only. It does not read any evaluation data.

cp2_seed <- 7072026L
cp2_families <- c("Elastic net", "GAM", "Random forest")
cp2_lambdas <- exp(seq(log(0.2), log(0.0005), length.out = 20L))
cp2_gamma <- c(1, 1.4)
cp2_rf_grid <- data.table::CJ(mtry = c(14L, 48L), min_node_size = c(10L, 30L))
cp2_smooth_vars <- c(
  "Age", "HR__mean", "MAP__mean", "Temp__mean", "GCS__min",
  "Creatinine__max", "BUN__max", "Lactate__max", "pH__min",
  "PaO2__min", "Glucose__max", "Urine__mean"
)
cp2_gam_linear <- c(
  "Gender_male", "Gender_unknown", "ICUType_2", "ICUType_3",
  "ICUType_4", "ICUType_unknown", "Lactate__count", "MechVent__count"
)

cp2_features <- function(wide) {
  w <- data.table::copy(wide)
  series_vars <- sort(sub("__count$", "", grep("__count$", names(w), value = TRUE)))
  summaries <- unlist(lapply(series_vars, function(v) {
    paste(v, c("last", "min", "max", "mean", "count"), sep = "__")
  }), use.names = FALSE)
  stopifnot(length(series_vars) == 37L, all(summaries %in% names(w)))
  x <- as.data.frame(w[, c("Age", "Height", "Weight", summaries), with = FALSE])
  for (n in grep("__count$", names(x), value = TRUE)) x[[n]] <- log1p(x[[n]])
  x$Gender_male <- as.numeric(!is.na(w$Gender) & w$Gender == 1)
  x$Gender_unknown <- as.numeric(is.na(w$Gender))
  for (i in 2:4) x[[paste0("ICUType_", i)]] <- as.numeric(!is.na(w$ICUType) & w$ICUType == i)
  x$ICUType_unknown <- as.numeric(is.na(w$ICUType))
  stopifnot(!any(c("RecordID", "InHospitalDeath", "SAPS_I", "SOFA", "Survival", "LengthOfStay") %in% names(x)))
  x
}

cp2_folds <- function(y, k, seed) {
  set.seed(seed)
  f <- integer(length(y))
  for (value in c(0L, 1L)) {
    ids <- which(y == value)
    f[ids] <- sample(rep(seq_len(k), length.out = length(ids)))
  }
  f
}

cp2_prep <- function(x) {
  medians <- vapply(x, function(z) if (all(is.na(z))) 0 else median(z, na.rm = TRUE), numeric(1))
  filled <- as.matrix(x)
  for (j in seq_len(ncol(filled))) filled[is.na(filled[, j]), j] <- medians[[j]]
  centers <- colMeans(filled)
  scales <- apply(filled, 2L, sd)
  keep <- names(scales)[is.finite(scales) & scales > 1e-8]
  list(input_columns = names(x), medians = medians, centers = centers, scales = scales, keep = keep)
}

cp2_bake <- function(prep, x) {
  stopifnot(identical(names(x), prep$input_columns))
  a <- as.matrix(x)
  for (j in seq_len(ncol(a))) a[is.na(a[, j]), j] <- prep$medians[[j]]
  a <- a[, prep$keep, drop = FALSE]
  a <- sweep(sweep(a, 2L, prep$centers[prep$keep], "-"), 2L, prep$scales[prep$keep], "/")
  stopifnot(all(is.finite(a)))
  a
}

cp2_auc <- function(y, score) {
  stopifnot(length(y) == length(score), !anyNA(y), all(is.finite(score)))
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cp2_gam_fit <- function(x, y, gamma) {
  smooths <- intersect(cp2_smooth_vars, colnames(x))
  linear <- intersect(cp2_gam_linear, colnames(x))
  formula <- as.formula(paste(".outcome ~", paste(
    c(sprintf("s(%s, k = 4, bs = 'cr')", smooths), linear), collapse = " + "
  )))
  dat <- as.data.frame(x[, c(smooths, linear), drop = FALSE])
  dat$.outcome <- y
  # A fixed, small ridge penalty prevents rare indicator levels separating the outcome.
  penalties <- setNames(lapply(linear, function(term) list(matrix(1, 1, 1), sp = 1)), linear)
  fit <- mgcv::gam(formula, data = dat, family = binomial(), method = "REML",
            select = TRUE, gamma = gamma,
            paraPen = penalties, control = mgcv::gam.control(maxit = 200L))
  if (!isTRUE(fit$converged)) stop("The GAM did not converge; evaluation has stopped.")
  if (!is.null(fit$outer.info$conv) && fit$outer.info$conv != "full convergence") {
    stop(paste("GAM smoothing parameter fit:", fit$outer.info$conv))
  }
  fit
}

cp2_rf_fit <- function(x, y, mtry, min_node_size, seed, importance = "none") {
  ranger::ranger(
    x = as.data.frame(x), y = factor(y, levels = c(0, 1), labels = c("Survived", "Died")),
    probability = TRUE, num.trees = 300L, mtry = min(mtry, ncol(x)),
    min.node.size = min_node_size, seed = seed, num.threads = 2L,
    importance = importance, write.forest = TRUE
  )
}

cp2_tune <- function(x, y, seed, stage) {
  inner <- cp2_folds(y, 3L, seed)
  records <- list()
  add <- function(family, config, fold, loss) {
    records[[length(records) + 1L]] <<- data.table::data.table(
      Stage = stage, Family = family, Config = config, InnerFold = fold, Brier = loss
    )
  }
  for (fold in 1:3) {
    train <- which(inner != fold)
    valid <- which(inner == fold)
    prep <- cp2_prep(x[train, , drop = FALSE])
    tr <- cp2_bake(prep, x[train, , drop = FALSE])
    va <- cp2_bake(prep, x[valid, , drop = FALSE])
    en <- glmnet::glmnet(tr, y[train], family = "binomial", alpha = 0.5,
                         lambda = cp2_lambdas, standardize = FALSE, control = list(maxit = 100000L))
    ep <- as.matrix(predict(en, newx = va, s = cp2_lambdas, type = "response"))
    for (j in seq_along(cp2_lambdas)) add("Elastic net", j, fold, mean((ep[, j] - y[valid])^2))
    for (j in seq_along(cp2_gamma)) {
      fit <- cp2_gam_fit(tr, y[train], cp2_gamma[[j]])
      p <- as.numeric(predict(fit, newdata = as.data.frame(va), type = "response"))
      add("GAM", j, fold, mean((p - y[valid])^2))
    }
    for (j in seq_len(nrow(cp2_rf_grid))) {
      g <- cp2_rf_grid[j]
      fit <- cp2_rf_fit(tr, y[train], g$mtry, g$min_node_size, seed + 100L * fold + j)
      p <- predict(fit, data = as.data.frame(va))$predictions[, "Died"]
      add("Random forest", j, fold, mean((p - y[valid])^2))
    }
  }
  results <- data.table::rbindlist(records)
  averages <- results[, .(MeanBrier = mean(Brier)), by = .(Family, Config)]
  data.table::setorder(averages, Family, MeanBrier, Config)
  selected <- averages[, .SD[1L], by = Family]
  list(selected = selected, results = results)
}

cp2_fit_models <- function(x, y, tuning, seed, importance = FALSE) {
  prep <- cp2_prep(x)
  a <- cp2_bake(prep, x)
  config <- function(f) tuning[Family == f, Config][[1L]]
  lambda <- cp2_lambdas[[config("Elastic net")]]
  en <- glmnet::glmnet(a, y, family = "binomial", alpha = 0.5, lambda = cp2_lambdas,
                       standardize = FALSE, control = list(maxit = 100000L))
  gam <- cp2_gam_fit(a, y, cp2_gamma[[config("GAM")]])
  g <- cp2_rf_grid[config("Random forest")]
  rf <- cp2_rf_fit(a, y, g$mtry, g$min_node_size, seed,
                  if (importance) "permutation" else "none")
  list(preprocessing = prep, models = list("Elastic net" = en, "GAM" = gam, "Random forest" = rf),
       lambda = lambda, gamma = cp2_gamma[[config("GAM")]], forest_config = g)
}

cp2_predict_all <- function(fitted, x) {
  a <- cp2_bake(fitted$preprocessing, x)
  list(
    "Elastic net" = as.numeric(predict(fitted$models[["Elastic net"]], newx = a,
                                       s = fitted$lambda, type = "response")),
    "GAM" = as.numeric(predict(fitted$models[["GAM"]], newdata = as.data.frame(a), type = "response")),
    "Random forest" = as.numeric(predict(fitted$models[["Random forest"]], data = as.data.frame(a))$predictions[, "Died"])
  )
}

# Use the frozen pipeline on a modeling table constructed with the same cleaning code.
cp2_predict_final <- function(bundle, wide) {
  package <- switch(bundle$family, "Elastic net" = "glmnet", "GAM" = "mgcv", "Random forest" = "ranger")
  if (!requireNamespace(package, quietly = TRUE)) stop(paste("Missing required package:", package))
  x <- cp2_features(wide)
  a <- cp2_bake(bundle$preprocessing, x)
  switch(bundle$family,
    "Elastic net" = as.numeric(predict(bundle$model, newx = a, s = bundle$lambda, type = "response")),
    "GAM" = as.numeric(predict(bundle$model, newdata = as.data.frame(a), type = "response")),
    "Random forest" = as.numeric(predict(bundle$model, data = as.data.frame(a))$predictions[, "Died"])
  )
}

cp2_coefficients <- function(fitted, fold) {
  b <- as.matrix(coef(fitted$models[["Elastic net"]], s = fitted$lambda))
  data.table::data.table(Feature = rownames(b), Coefficient = b[, 1L], Fold = fold)[Feature != "(Intercept)"]
}

cp2_roc <- function(y, score) {
  d <- data.table::data.table(y = y, score = score)
  d <- d[, .(TP = sum(y == 1), FP = sum(y == 0)), by = score][order(-score)]
  data.table::rbindlist(list(data.table::data.table(FPR = 0, TPR = 0),
    d[, .(FPR = cumsum(FP) / sum(FP), TPR = cumsum(TP) / sum(TP))]))
}
