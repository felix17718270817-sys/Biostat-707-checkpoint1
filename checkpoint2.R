# BIOSTAT 707 Checkpoint 2. Run: pixi run --locked checkpoint2
# Rebuilds set-a features, performs nested CV, saves models, and creates the PDF.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mgcv)
})
here::i_am("checkpoint2.R")
root <- here::here()
out_dir <- file.path(root, "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
started <- proc.time()[["elapsed"]]
set.seed(7072026L)
data.table::setDTthreads(2L)
source(file.path(root, "R", "cp2_models.R"))
source(file.path(root, "R", "cp2_report.R"))

message("Rebuilding set-a tables with the Checkpoint 1 cleaning rules...")
cp1 <- new.env(parent = globalenv())
sys.source(file.path(root, "checkpoint1.R"), envir = cp1)
w <- fread(file.path(out_dir, "set-a_wide.csv"))
stopifnot(nrow(w) == 4000L, uniqueN(w$RecordID) == 4000L, !anyNA(w$InHospitalDeath))
stopifnot(all(w$InHospitalDeath %in% c(0L, 1L)))
x <- cp2_features(w)
y <- w$InHospitalDeath
folds <- cp2_folds(y, 5L, cp2_seed)
fwrite(data.table(RecordID = w$RecordID, Fold = folds, Outcome = y), file.path(out_dir, "cp2_folds.csv"))

predictions <- list()
coefs <- list()
tuning_results <- list()
tuning_choices <- list()
audits <- list()
for (fold in 1:5) {
  message(sprintf("Outer fold %d/5: inner tuning and held-out prediction...", fold))
  train <- which(folds != fold)
  test <- which(folds == fold)
  stopifnot(length(intersect(w$RecordID[train], w$RecordID[test])) == 0L)
  tuned <- cp2_tune(x[train, , drop = FALSE], y[train], cp2_seed + fold, paste0("Outer", fold))
  tuning_results[[fold]] <- tuned$results
  tuning_choices[[fold]] <- copy(tuned$selected)[, OuterFold := fold]
  fit <- cp2_fit_models(x[train, , drop = FALSE], y[train], tuned$selected, cp2_seed + 1000L + fold)
  p <- cp2_predict_all(fit, x[test, , drop = FALSE])
  for (benchmark in c("SAPS_I", "SOFA")) {
    score <- w[[benchmark]][test]
    score[is.na(score)] <- median(w[[benchmark]][train], na.rm = TRUE)
    p[[if (benchmark == "SAPS_I") "SAPS-I" else "SOFA"]] <- score
  }
  predictions[[fold]] <- rbindlist(lapply(names(p), function(model) {
    stopifnot(all(is.finite(p[[model]])))
    data.table(RecordID = w$RecordID[test], Fold = fold, Outcome = y[test], Model = model, Score = p[[model]])
  }))
  coefs[[fold]] <- cp2_coefficients(fit, fold)
  audits[[fold]] <- data.table(
    Fold = fold, TrainingStays = length(train), ValidationStays = length(test),
    TrainingDeaths = sum(y[train]), ValidationDeaths = sum(y[test]), Overlap = 0L,
    RetainedColumns = length(fit$preprocessing$keep)
  )
}
oof <- rbindlist(predictions)
stopifnot(nrow(oof) == 20000L, all(oof[, .N, by = .(RecordID, Model)]$N == 1L))
stopifnot(all(oof[Model %chin% cp2_families, Score] >= 0 & oof[Model %chin% cp2_families, Score] <= 1))
fold_metrics <- oof[, .(AUROC = cp2_auc(Outcome, Score)), by = .(Model, Fold)]
metrics <- fold_metrics[, .(MeanAUROC = mean(AUROC), SDAUROC = sd(AUROC)), by = Model]
metrics <- metrics[match(c(cp2_families, "SAPS-I", "SOFA"), Model)]
comparison <- merge(dcast(fold_metrics, Model ~ Fold, value.var = "AUROC"), metrics, by = "Model")
comparison <- comparison[match(c(cp2_families, "SAPS-I", "SOFA"), Model)]
fwrite(comparison, file.path(out_dir, "cp2_comparison.csv"))
fwrite(oof, file.path(out_dir, "cp2_oof_predictions.csv"))
fwrite(rbindlist(audits), file.path(out_dir, "cp2_cv_audit.csv"))

# A predeclared one-standard-error heuristic favors the simpler model family.
model_metrics <- metrics[Model %chin% cp2_families]
best <- model_metrics[which.max(MeanAUROC)]
cutoff <- best$MeanAUROC - best$SDAUROC / sqrt(5)
eligible <- model_metrics[MeanAUROC >= cutoff, Model]
chosen <- cp2_families[cp2_families %chin% eligible][[1L]]
message(sprintf("Selected family: %s. Tuning and refitting on all set-a...", chosen))
full_tuning <- cp2_tune(x, y, cp2_seed + 900L, "AllSetA")
all_fit <- cp2_fit_models(x, y, full_tuning$selected, cp2_seed + 10000L, importance = TRUE)
fwrite(rbindlist(c(tuning_results, list(full_tuning$results))), file.path(out_dir, "cp2_inner_tuning.csv"))
fwrite(rbindlist(c(tuning_choices, list(copy(full_tuning$selected)[, OuterFold := 0L]))),
       file.path(out_dir, "cp2_selected_parameters.csv"))

coefficient_table <- cp2_coefficients(all_fit, 0L)
fold_coefficients <- rbindlist(coefs)
stability <- fold_coefficients[, .(FoldsSelected = sum(abs(Coefficient) > 1e-8)), by = Feature]
stability <- merge(coefficient_table[, .(Feature, Coefficient)], stability, by = "Feature", all.x = TRUE)
stability[is.na(FoldsSelected), FoldsSelected := 0L]
stability[, SelectionRate := FoldsSelected / 5]
stability <- stability[order(-abs(Coefficient))]
fwrite(stability, file.path(out_dir, "cp2_elastic_net_coefficients.csv"))
fwrite(fold_coefficients, file.path(out_dir, "cp2_fold_coefficients.csv"))
fwrite(as.data.table(summary(all_fit$models[["GAM"]])$s.table, keep.rownames = "Term"),
       file.path(out_dir, "cp2_gam_smooths.csv"))

bundle <- list(
  family = chosen, model = all_fit$models[[chosen]], preprocessing = all_fit$preprocessing,
  lambda = all_fit$lambda, gamma = all_fit$gamma, forest_config = all_fit$forest_config,
  predictors = names(x), seed = cp2_seed, trained_on = "set-a only",
  training_record_ids = w$RecordID, feature_code = "R/cp2_models.R",
  cleaning_code = "checkpoint1.R", time_window_hours = c(0, 48), cv_metrics = metrics
)
saveRDS(bundle, file.path(out_dir, "cp2_final_model.rds"))
saveRDS(all_fit, file.path(out_dir, "cp2_all_models.rds"))
reloaded <- readRDS(file.path(out_dir, "cp2_final_model.rds"))
reload_predictions <- cp2_predict_final(reloaded, w[1:20])
original_predictions <- cp2_predict_all(all_fit, x[1:20, , drop = FALSE])[[chosen]]
stopifnot(isTRUE(all.equal(reload_predictions, original_predictions, tolerance = 1e-12)))

checks <- data.table(
  Check = c("Unique set-a stays", "OOF rows (5 methods x 4000 stays)",
            "One OOF prediction per method and stay", "Training/validation ID overlap",
            "Final model reload reproduces predictions"),
  Passed = c(nrow(w) == 4000L && uniqueN(w$RecordID) == 4000L, nrow(oof) == 20000L,
             all(oof[, .N, by = .(Model, RecordID)]$N == 1L),
             all(rbindlist(audits)$Overlap == 0L), TRUE)
)
fwrite(checks, file.path(out_dir, "cp2_validation.csv"))
stopifnot(all(checks$Passed))
cp2_make_report(w, x, oof, comparison, metrics, stability, all_fit, chosen, cutoff,
                rbindlist(audits), out_dir)
elapsed <- proc.time()[["elapsed"]] - started
writeLines(c(sprintf("Elapsed seconds: %.1f", elapsed), capture.output(sessionInfo())),
           file.path(out_dir, "cp2_session_info.txt"))
message(sprintf("Checkpoint 2 completed in %.1f minutes. Report: output/checkpoint2.pdf", elapsed / 60))
