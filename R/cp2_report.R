# Generate a standalone PDF report with R graphics; no TeX installation is needed.

cp2_make_report <- function(w, x, oof, comparison, metrics, stability, fitted,
                            chosen, cutoff, audit, out_dir) {
  theme_set(theme_minimal(base_size = 10, base_family = "sans"))
  colors <- c("Elastic net" = "#2563EB", "GAM" = "#B45309", "Random forest" = "#047857",
              "SAPS-I" = "#9333EA", "SOFA" = "#B91C1C")
  roc_data <- rbindlist(lapply(c(cp2_families, "SAPS-I", "SOFA"), function(model) {
    d <- oof[Model == model]
    cp2_roc(d$Outcome, d$Score)[, Model := model]
  }))
  roc_plot <- ggplot(roc_data, aes(FPR, TPR, color = Model)) +
    geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = 2) +
    geom_path(linewidth = 0.65) + scale_color_manual(values = colors) +
    labs(x = "False positive rate", y = "Sensitivity", color = NULL,
         title = "Figure 1. ROC curves from held-out fold predictions") +
    coord_equal() + theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, "cp2_roc.png"), roc_plot, width = 7, height = 4.8, dpi = 180)

  coefficient_data <- stability[abs(Coefficient) > 1e-8][1:min(18L, .N)]
  coefficient_data[, Feature := factor(Feature, levels = rev(Feature))]
  coefficient_plot <- ggplot(coefficient_data, aes(Coefficient, Feature)) +
    geom_vline(xintercept = 0, color = "grey65") +
    geom_segment(aes(x = 0, xend = Coefficient, yend = Feature), color = "#2563EB") +
    geom_point(color = "#2563EB", size = 2) +
    labs(x = "Log-odds coefficient per training-set standard deviation", y = NULL,
         title = "Figure 2. Largest elastic-net coefficients") +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
  ggsave(file.path(out_dir, "cp2_coefficients.png"), coefficient_plot, width = 8, height = 5.2, dpi = 180)

  prep <- fitted$preprocessing
  a <- cp2_bake(prep, x)
  base <- as.data.frame(a[rep(1L, 100L), , drop = FALSE])
  for (n in names(base)) base[[n]] <- 0
  smooth_results <- list()
  smooth_table <- summary(fitted$models[["GAM"]])$s.table
  for (variable in c("Age", "HR__mean", "Temp__mean")) {
    limits <- quantile(x[[variable]], c(0.025, 0.975), na.rm = TRUE)
    values <- seq(limits[[1L]], limits[[2L]], length.out = 100L)
    newdata <- base
    newdata[[variable]] <- (values - prep$centers[[variable]]) / prep$scales[[variable]]
    p <- predict(fitted$models[["GAM"]], newdata = newdata, type = "terms", se.fit = TRUE)
    term <- paste0("s(", variable, ")")
    smooth_results[[variable]] <- data.table(
      Variable = variable, X = values, Effect = p$fit[, term], SE = p$se.fit[, term],
      EDF = smooth_table[term, "edf"]
    )
  }
  smooth_data <- rbindlist(smooth_results)
  smooth_data[, Label := factor(Variable,
    levels = c("Age", "HR__mean", "Temp__mean"),
    labels = c("Age (years)", "Mean heart rate (beats/min)", "Mean temperature (degrees C)"))]
  smooth_plot <- ggplot(smooth_data, aes(X, Effect)) +
    geom_ribbon(aes(ymin = Effect - 1.96 * SE, ymax = Effect + 1.96 * SE), fill = "#DBEAFE") +
    geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
    geom_line(color = "#2563EB", linewidth = 0.7) +
    facet_wrap(~Label, ncol = 1, scales = "free") +
    labs(x = NULL, y = "Centered contribution to log odds",
         title = "Figure 3. GAM partial effects with approximate 95% bands") +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
  ggsave(file.path(out_dir, "cp2_gam_effects.png"), smooth_plot, width = 7, height = 7, dpi = 180)
  fwrite(smooth_data[, !"Label"], file.path(out_dir, "cp2_gam_plot_data.csv"))

  group_table <- data.frame(
    Group = c("Admission variables", "Last values", "Minimum and maximum", "Mean values", "Measurement counts"),
    Columns = c(9, 37, 74, 37, 37),
    Reason = c("Age, body size, sex, and ICU case mix",
               "Condition nearest the end of the window", "Most abnormal states during the window",
               "Typical measured state across 48 hours", "Monitoring intensity and never measured status")
  )
  catalog <- data.table(Feature = names(x))
  catalog[, Group := fcase(grepl("__count$", Feature), "Measurement counts",
    grepl("__last$", Feature), "Last values", grepl("__(min|max)$", Feature), "Minimum and maximum",
    grepl("__mean$", Feature), "Mean values", default = "Admission variables")]
  fwrite(catalog, file.path(out_dir, "cp2_predictor_catalog.csv"))

  page <- 0L
  page_start <- function(title) {
    page <<- page + 1L
    grid::grid.newpage()
    grid::grid.text(title, x = 0.08, y = 0.955, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 17, fontface = "bold", col = "#172554"))
    grid::grid.lines(x = c(0.08, 0.92), y = c(0.918, 0.918), gp = grid::gpar(col = "#CBD5E1"))
    grid::grid.text("BIOSTAT 707 - Checkpoint 2 - Felix Wang", x = 0.08, y = 0.035,
                    just = "left", gp = grid::gpar(fontsize = 8, col = "#475569"))
    grid::grid.text(as.character(page), x = 0.92, y = 0.035,
                    just = "right", gp = grid::gpar(fontsize = 8, col = "#475569"))
  }
  heading <- function(text, y) {
    grid::grid.text(text, x = 0.08, y = y, just = c("left", "top"),
                    gp = grid::gpar(fontsize = 11.5, fontface = "bold"))
  }
  paragraph <- function(text, y, size = 10, width = 99L) {
    lines <- strwrap(text, width = width)
    grid::grid.text(paste(lines, collapse = "\n"), x = 0.08, y = y, just = c("left", "top"),
                    gp = grid::gpar(fontsize = size, lineheight = 1.25))
    invisible(y - length(lines) * 0.018)
  }
  draw_table <- function(data, y, widths, size = 9) {
    t <- gridExtra::tableGrob(data, rows = NULL, theme = gridExtra::ttheme_minimal(
      base_size = size, padding = grid::unit(c(3.5, 4), "mm"),
      core = list(fg_params = list(hjust = 0, x = 0.04),
                  bg_params = list(fill = rep(c("white", "#F1F5F9"), length.out = nrow(data)), col = NA)),
      colhead = list(fg_params = list(fontface = "bold", hjust = 0, x = 0.04),
                     bg_params = list(fill = "#DBEAFE", col = NA))
    ))
    t$widths <- grid::unit(widths, "inches")
    grid::pushViewport(grid::viewport(x = 0.5, y = y))
    grid::grid.draw(t)
    grid::popViewport()
  }
  draw_plot <- function(plot, y, height) {
    print(plot, newpage = FALSE, vp = grid::viewport(x = 0.5, y = y, width = 0.86, height = height))
  }

  grDevices::pdf(file.path(out_dir, "checkpoint2.pdf"), width = 8.5, height = 11,
                  family = "Helvetica", onefile = TRUE, timestamp = FALSE,
                  author = "Felix Wang", title = "BIOSTAT 707 Checkpoint 2")
  on.exit(grDevices::dev.off(), add = TRUE)

  page_start("Feature engineering and supervised modeling")
  heading("AI-use statement", 0.88)
  paragraph(paste(
    "I used OpenAI Codex for coding, debugging, running checks, and drafting this report.",
    "Outputs were checked against the data, saved predictions, and figures. It saved time, but",
    "code that runs can still make statistical mistakes. I need to understand and review the analysis."
  ), 0.85)
  heading("2.1 Modeling table and prediction target", 0.72)
  paragraph(sprintf(paste(
    "I used set-a only: %s intensive care unit (ICU) stays, including %s in-hospital deaths (%.2f%%).",
    "Each row is one stay, and the script checks that RecordID is unique. The outcome is",
    "in-hospital death, joined by RecordID. Predictors summarize the first 48 ICU hours.",
    "This gives a prediction after the full observation window, not at admission."
  ), format(nrow(w), big.mark = ","), sum(w$InHospitalDeath), 100 * mean(w$InHospitalDeath)), 0.69)
  heading("Cleaning carried over from Checkpoint 1", 0.55)
  paragraph(paste(
    "The script rebuilds the long and wide tables from the original patient files on every run.",
    "It carries over the Checkpoint 1 cleaning rules. Missing codes of -1 are set to missing",
    "where appropriate. Original values remain in ValueRaw. Screening ranges include height",
    "100-250 cm, weight 20-300 kg, pH 6.5-8.0, and temperature 20-45 degrees C.",
    "These rules also remove zero height and weight values. Some extreme but non-negative",
    "laboratory and vital-sign values remain. This is a limitation rather than proof they are correct."
  ), 0.52)
  heading("What is excluded", 0.32)
  paragraph(paste(
    "RecordID is used only to identify stays and check splits. Hospital length of stay and",
    "survival time are excluded because they contain later outcome information.",
    "The Simplified Acute Physiology Score I (SAPS-I) and Sequential Organ Failure Assessment",
    "(SOFA) score are benchmarks, not inputs to the three models. No set-b or set-c files are read."
  ), 0.29)
  paragraph("Data source: PhysioNet/Computing in Cardiology Challenge 2012, version 1.0.0.", 0.12, size = 9)

  page_start("2.2 Predictors")
  paragraph(sprintf(paste(
    "Table 1 groups the %d input columns before training-fold zero-variance removal.",
    "ICU type is represented by indicator columns rather than treated as a continuous number.",
    "Gender also has an unknown category. The full column list is in cp2_predictor_catalog.csv."
  ), ncol(x)), 0.88)
  heading("Table 1. Predictor groups and clinical reasons", 0.77)
  draw_table(group_table, 0.63, c(1.6, 0.65, 4.7), size = 8.2)
  series_vars <- sort(sub("__count$", "", grep("__count$", names(x), value = TRUE)))
  heading("Measurements included", 0.47)
  paragraph(paste(paste(series_vars, collapse = ", "),
    ". Each has a last value, minimum, maximum, mean, and count."), 0.44, size = 9.5)
  heading("Missingness and preprocessing", 0.29)
  paragraph(paste(
    "No clinical variable is dropped for low coverage. Counts are transformed with log(1 + count),",
    "so zero means no valid measurement and high counts have less influence. First values are omitted",
    "to reduce redundancy. Missing numeric values are filled with the training-fold median.",
    "Columns with zero variance are removed using the training fold. Remaining columns are",
    "centered and scaled using that fold. The held-out fold never sets these values."
  ), 0.26)
  paragraph(paste("The generalized additive model (GAM) uses 12 clinically chosen continuous predictors, sex, ICU type, and",
    "lactate and ventilation measurement counts. This keeps the nonlinear model manageable."), 0.105, size = 9.5)

  page_start("2.3 Cross-validation and tuning")
  paragraph(paste(
    "I used five-fold stratified outer cross-validation (CV), with the same stay-level folds for",
    "every method. Each patient receives one out-of-fold prediction per method. In each outer",
    "training set, three-fold stratified inner CV chooses hyperparameters by mean Brier score",
    "(mean squared error of predicted probability). Preprocessing is fitted again inside every",
    "inner training fold. This prevents imputation and scaling from leaking across either level."
  ), 0.88)
  heading("Table 2. Prespecified tuning search", 0.72)
  tuning_table <- data.frame(
    Model = c("Elastic net", "GAM", "Random forest"),
    Search = c("alpha = 0.5; 20 lambdas: 0.0005 to 0.2",
               "12 cubic smooths, k = 4; gamma = 1 or 1.4",
               "300 trees; mtry = 14 or 48; min node = 10 or 30")
  )
  draw_table(tuning_table, 0.61, c(1.15, 5.8), size = 9)
  paragraph(paste(
    "Elastic net combines shrinkage with variable selection and is useful for correlated",
    "summaries of the same measurement. The generalized additive model (GAM) uses shrinkage",
    "smooths fitted by restricted maximum likelihood, with a fixed ridge penalty on linear terms.",
    "The forest uses median imputation, no early stopping, no oversampling, and no class weighting."
  ), 0.48)
  heading("Table 3. Outer-fold audit", 0.35)
  audit_print <- as.data.frame(audit[, .(Fold, TrainingStays, ValidationStays, ValidationDeaths, Overlap)])
  names(audit_print) <- c("Fold", "Training stays", "Held-out stays", "Held-out deaths", "Shared IDs")
  draw_table(audit_print, 0.23, c(0.55, 1.65, 1.65, 1.65, 1.45), size = 8.5)
  paragraph("All folds keep stays separate. Random seed: 7072026. Inner folds are shared across the three models.", 0.085, size = 9)

  page_start("2.7 Comparing models")
  paragraph(paste(
    "The area under the receiver operating characteristic (ROC) curve, or AUROC, measures ranking.",
    "A value of 0.5 is random ranking and 1 is perfect. Table 4 gives each fold, its mean, and standard deviation (SD).",
    "Missing SAPS-I or SOFA scores are filled with the outer training-set median for comparison",
    "on the same patients. These benchmarks are raw scores, not calibrated probabilities."
  ), 0.88)
  comparison_print <- as.data.frame(comparison)
  for (j in 2:ncol(comparison_print)) comparison_print[[j]] <- sprintf("%.3f", comparison_print[[j]])
  names(comparison_print) <- c("Method", paste0("F", 1:5), "Mean", "SD")
  heading("Table 4. Held-out AUROC by fold", 0.755)
  draw_table(comparison_print, 0.62, c(1.45, rep(0.7, 5), 0.85, 0.65), size = 8.5)
  draw_plot(roc_plot, 0.31, 0.40)
  best_model <- metrics[Model %chin% cp2_families][which.max(MeanAUROC)]
  benchmark_best <- metrics[Model %chin% c("SAPS-I", "SOFA")][which.max(MeanAUROC)]
  paragraph(sprintf(paste(
    "Figure 1 uses out-of-fold predictions for all five methods. %s has the highest model mean",
    "AUROC (%.3f). The best benchmark is %s (%.3f); the gap is %.3f, compared with fold SDs of",
    "%.3f and %.3f. The gain exceeds this fold variation, but this is not a significance test."
  ), best_model$Model, best_model$MeanAUROC, benchmark_best$Model, benchmark_best$MeanAUROC,
     best_model$MeanAUROC - benchmark_best$MeanAUROC, best_model$SDAUROC, benchmark_best$SDAUROC), 0.105, size = 9)

  page_start("2.4 Regularized logistic regression")
  selected_n <- sum(abs(stability$Coefficient) > 1e-8)
  paragraph(sprintf(paste(
    "The refit elastic-net model uses alpha = 0.5 and lambda = %.5f, selected by inner CV",
    "on all set-a. It retains %d nonzero coefficients, listed in Appendix Table A1.",
    "Figure 2 shows the 18 largest on a common scale. Positive values increase predicted log odds,",
    "conditional on the other inputs; they do not identify treatment effects."
  ), fitted$lambda, selected_n), 0.88)
  draw_plot(coefficient_plot, 0.56, 0.48)
  heading("Table 5. Selection stability for the largest coefficients", 0.305)
  st <- as.data.frame(stability[1:6, .(Feature, Coefficient, FoldsSelected)])
  st$Coefficient <- sprintf("%.3f", st$Coefficient)
  names(st) <- c("Predictor", "Coefficient", "Folds selected / 5")
  draw_table(st, 0.19, c(4.1, 1.1, 1.75), size = 8.5)
  stable_n <- stability[abs(Coefficient) > 1e-8 & FoldsSelected == 5L, .N]
  paragraph(sprintf(paste("%d of the %d full-data selected inputs appear in all five outer-fold fits.",
    "Correlated summaries can replace each other, so unstable selection limits interpretation."), stable_n, selected_n), 0.09, size = 8.8)

  page_start("2.5 Nonlinear effects")
  paragraph(paste(
    "Figure 3 shows centered smooth contributions to log odds from the refit GAM, with approximate",
    "95% confidence bands. The x ranges cover the 2.5th to 97.5th observed percentiles.",
    "The plots describe associations while the other inputs are held fixed; zero is a centered",
    "reference, not an absolute risk of zero. All outer-fold predictions used fold-specific fits."
  ), 0.88)
  draw_plot(smooth_plot, 0.48, 0.64)
  descriptions <- vapply(names(smooth_results), function(variable) {
    z <- smooth_results[[variable]]
    low <- which.min(z$Effect)
    high <- which.max(z$Effect)
    shape <- if (low > 10 && low < 90 && min(z$Effect[c(1, 100)]) > min(z$Effect) + 0.1) {
      "U-shaped"
    } else if (high > 10 && high < 90 && max(z$Effect[c(1, 100)]) < max(z$Effect) - 0.1) {
      "inverted U-shaped"
    } else if (tail(z$Effect, 1) > z$Effect[[1]]) "mostly increasing" else "mostly decreasing"
    label <- c(Age = "Age", HR__mean = "Heart rate", Temp__mean = "Temperature")[[variable]]
    sprintf("%s is %s (effective df %.2f).", label, shape, z$EDF[[1]])
  }, character(1))
  paragraph(paste(paste(descriptions, collapse = " "),
    "Higher age and fast heart rate are plausible severity markers. Higher risk at low temperature",
    "is also plausible in critical illness. Effective degrees of freedom above 1 suggest nonlinearity.",
    "These are associations, with uncertainty at extreme values."),
    0.13, size = 9)

  page_start("Final model and limitations")
  chosen_metrics <- metrics[Model == chosen]
  heading("Choice and refit", 0.88)
  paragraph(sprintf(paste(
    "I chose %s, with mean outer-fold AUROC %.3f (SD %.3f). My rule first finds the largest",
    "mean AUROC, then accepts methods within one standard error of that mean (cutoff %.3f).",
    "Among eligible methods, it favors elastic net, then GAM, then random forest for simplicity.",
    "This is a model-selection heuristic, not a statistical equivalence test."
  ), chosen, chosen_metrics$MeanAUROC, chosen_metrics$SDAUROC, cutoff), 0.85)
  paragraph(paste(
    "After choosing the family, I repeated the three-fold parameter search on all set-a and refit",
    "the model on all 4,000 stays. The saved cp2_final_model.rds contains the model, feature order,",
    "imputation medians, scaling values, and retained columns. Checkpoint 3 can load this frozen",
    "pipeline. cp2_oof_predictions.csv saves held-out predictions for all five methods."
  ), 0.69)
  heading("What the results do and do not establish", 0.54)
  paragraph(paste(
    "Inner tuning avoids choosing parameters with the outer validation outcomes. However, the",
    "same outer results choose the final model family, so its reported winning performance may",
    "still be optimistic. Folds are from one data source and do not test transport to a new hospital.",
    "Median imputation is simple but may lose detail. Some remaining extreme values may be",
    "measurement errors. Counts and ventilation can reflect care decisions as well as illness.",
    "A predictor being important does not mean that intervening on it would reduce mortality.",
    "Calibration and clinical usefulness will need separate evaluation in Checkpoint 3."
  ), 0.51)
  heading("Validation and reproducibility", 0.315)
  paragraph(paste(
    "Automated checks cover unique stays, disjoint training and validation IDs, one held-out",
    "prediction per method and stay, probability bounds, and identical predictions after saving",
    "and reloading the final pipeline. The checkpoint2 pixi task rebuilds all derived files and",
    "this PDF. Model and table outputs remain local; the PDF is committed with the source code."
  ), 0.285)
  heading("References", 0.15)
  paragraph(paste(
    "Silva I et al. (2012). Predicting in-hospital mortality of ICU patients. Computing in",
    "Cardiology 39:245-248. Data: https://physionet.org/content/challenge-2012/1.0.0/.",
    "Course checklist: https://steverozen.github.io/biostat707_2026/checkpoint2-3-checklist.html."
  ), 0.12, size = 8.5, width = 110L)

  retained <- stability[abs(Coefficient) > 1e-8]
  chunks <- split(seq_len(nrow(retained)), ceiling(seq_len(nrow(retained)) / 28L))
  for (j in seq_along(chunks)) {
    page_start(sprintf("Appendix: retained predictors (%d/%d)", j, length(chunks)))
    paragraph(paste(
      "Table A1 lists every nonzero elastic-net coefficient, ordered by absolute size.",
      "The fold count shows how often each input was selected across the five outer fits.",
      "Predictors not listed have coefficient zero in the refit model."
    ), 0.88, size = 9)
    d <- as.data.frame(retained[chunks[[j]], .(Feature, Coefficient, FoldsSelected)])
    d$Coefficient <- sprintf("%.4f", d$Coefficient)
    names(d) <- c("Predictor", "Coefficient", "Folds selected / 5")
    draw_table(d, 0.46, c(4.1, 1.15, 1.7), size = 8.3)
  }
}
