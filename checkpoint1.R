# BIOSTAT 707 Checkpoint 1: cohort characterization and EDA for set-a.
# Run with: pixi run --locked checkpoint1
# The script reads data/ and writes all derived files to output/.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

here::i_am("checkpoint1.R")
root <- here::here()

data_dir <- file.path(root, "data")
set_a_dir <- file.path(data_dir, "set-a")
out_dir <- file.path(root, "output")
expected_n <- 4000L
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stop_with_message <- function(message) {
  stop(message, call. = FALSE)
}

check_input_data <- function() {
  outcome_path <- file.path(data_dir, "Outcomes-a.txt")
  files <- list.files(set_a_dir, pattern = "^[0-9]+\\.txt$", full.names = TRUE)

  if (!file.exists(outcome_path)) {
    stop_with_message("Missing data/Outcomes-a.txt.")
  }
  if (length(files) != expected_n) {
    stop_with_message(
      sprintf("data/set-a contains %d record files. Expected %d.", length(files), expected_n)
    )
  }
  files
}

read_one_record <- function(path) {
  record_id <- as.integer(tools::file_path_sans_ext(basename(path)))
  x <- fread(path, showProgress = FALSE)
  required <- c("Time", "Parameter", "Value")
  if (!identical(names(x), required)) {
    stop_with_message(sprintf("Unexpected columns in %s.", basename(path)))
  }

  id_row <- x[Parameter == "RecordID", Value]
  if (length(id_row) != 1L || as.integer(id_row) != record_id) {
    stop_with_message(sprintf("RecordID check failed in %s.", basename(path)))
  }

  x <- x[Parameter != "RecordID"]
  x[, RecordID := record_id]
  setcolorder(x, c("RecordID", "Time", "Parameter", "Value"))
  x
}

add_time_variables <- function(x) {
  time_parts <- tstrsplit(x$Time, ":", fixed = TRUE)
  x[, TimeMinutes := as.integer(time_parts[[1L]]) * 60L + as.integer(time_parts[[2L]])]
  x[, TimeHours := TimeMinutes / 60]
  x
}

clean_measurements <- function(x) {
  x[, ValueRaw := Value]
  x[, `:=`(is_missing_code = ValueRaw == -1, is_invalid = FALSE, invalid_reason = NA_character_)]
  x[is_missing_code == TRUE, invalid_reason := "-1 missing code"]
  x[!is_missing_code & ValueRaw < 0, `:=`(
    is_invalid = TRUE,
    invalid_reason = "negative value is not valid in this data set"
  )]

  range_rules <- data.table(
    Parameter = c(
      "Gender", "Height", "ICUType", "Weight", "FiO2", "GCS",
      "MechVent", "pH", "SaO2", "Temp"
    ),
    lower = c(0, 100, 1, 20, 0, 3, 0, 6.5, 0, 20),
    upper = c(1, 250, 4, 300, 1, 15, 1, 8.0, 100, 45)
  )

  x <- range_rules[x, on = "Parameter"]
  x[
    !is.na(lower) & !is_missing_code & (ValueRaw < lower | ValueRaw > upper),
    `:=`(is_invalid = TRUE, invalid_reason = "outside variable definition")
  ]
  x[, Value := ValueRaw]
  x[is_missing_code == TRUE | is_invalid == TRUE, Value := NA_real_]
  x[, c("lower", "upper") := NULL]
  setcolorder(
    x,
    c(
      "RecordID", "Time", "TimeMinutes", "TimeHours", "Parameter",
      "ValueRaw", "Value", "is_missing_code", "is_invalid", "invalid_reason"
    )
  )
  setorder(x, RecordID, TimeMinutes, Parameter)
  x
}

make_quality_summary <- function(long_data) {
  long_data[, .(
    n_rows = .N,
    n_missing_code = sum(is_missing_code),
    n_invalid = sum(is_invalid),
    raw_min = if (all(is.na(ValueRaw))) NA_real_ else min(ValueRaw, na.rm = TRUE),
    raw_max = if (all(is.na(ValueRaw))) NA_real_ else max(ValueRaw, na.rm = TRUE),
    cleaned_min = if (all(is.na(Value))) NA_real_ else min(Value, na.rm = TRUE),
    cleaned_max = if (all(is.na(Value))) NA_real_ else max(Value, na.rm = TRUE)
  ), by = Parameter][order(Parameter)]
}

make_static_table <- function(long_data) {
  static_variables <- c("Age", "Gender", "Height", "ICUType", "Weight")
  static_long <- long_data[Time == "00:00" & Parameter %chin% static_variables]
  static <- dcast(
    static_long,
    RecordID ~ Parameter,
    value.var = "Value",
    fun.aggregate = function(z) if (length(z) == 0L) NA_real_ else z[[1L]]
  )
  missing_columns <- setdiff(static_variables, names(static))
  for (column in missing_columns) static[, (column) := NA_real_]
  setcolorder(static, c("RecordID", static_variables))
  static
}

make_time_series_summary <- function(long_data) {
  descriptor_only <- c("Age", "Gender", "Height", "ICUType")
  time_series <- long_data[!Parameter %chin% descriptor_only]
  valid <- time_series[!is.na(Value)][order(RecordID, Parameter, TimeMinutes)]

  summary_long <- valid[, .(
    count = as.numeric(.N),
    first = first(Value),
    last = last(Value),
    min = min(Value),
    max = max(Value),
    mean = mean(Value)
  ), by = .(RecordID, Parameter)]

  summary_melted <- melt(
    summary_long,
    id.vars = c("RecordID", "Parameter"),
    variable.name = "Statistic",
    value.name = "SummaryValue"
  )
  summary_melted[, Feature := paste(Parameter, Statistic, sep = "__")]
  wide <- dcast(summary_melted, RecordID ~ Feature, value.var = "SummaryValue")

  count_columns <- grep("__count$", names(wide), value = TRUE)
  for (column in count_columns) set(wide, which(is.na(wide[[column]])), column, 0)
  wide
}

format_continuous <- function(x, digits = 1L) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return("NA")
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE, type = 2)
  sprintf(
    paste0("%.", digits, "f [ %.", digits, "f, %.", digits, "f ]"),
    q[[2L]], q[[1L]], q[[3L]]
  )
}

format_categorical <- function(x, value) {
  denominator <- sum(!is.na(x))
  numerator <- sum(x == value, na.rm = TRUE)
  if (denominator == 0L) return("0 (NA%)")
  sprintf("%d (%.1f%%)", numerator, 100 * numerator / denominator)
}

make_table1 <- function(cohort) {
  groups <- list(
    cohort,
    cohort[InHospitalDeath == 0],
    cohort[InHospitalDeath == 1]
  )
  output <- data.table(
    Variable = character(),
    LevelOrStatistic = character(),
    Overall = character(),
    Survived = character(),
    Died = character()
  )

  add_continuous <- function(variable, label, digits = 1L) {
    values <- vapply(groups, function(g) format_continuous(g[[variable]], digits), character(1))
    missing <- vapply(groups, function(g) sprintf("%d", sum(is.na(g[[variable]]))), character(1))
    output <<- rbind(
      output,
      data.table(
        Variable = label,
        LevelOrStatistic = "Median [Q1, Q3]",
        Overall = values[[1L]], Survived = values[[2L]], Died = values[[3L]]
      ),
      data.table(
        Variable = label,
        LevelOrStatistic = "Missing, n",
        Overall = missing[[1L]], Survived = missing[[2L]], Died = missing[[3L]]
      )
    )
  }

  add_category <- function(variable, label, value, level_label) {
    values <- vapply(groups, function(g) format_categorical(g[[variable]], value), character(1))
    output <<- rbind(
      output,
      data.table(
        Variable = label,
        LevelOrStatistic = level_label,
        Overall = values[[1L]], Survived = values[[2L]], Died = values[[3L]]
      )
    )
  }

  add_category_missing <- function(variable, label) {
    values <- vapply(groups, function(g) sprintf("%d", sum(is.na(g[[variable]]))), character(1))
    output <<- rbind(
      output,
      data.table(
        Variable = label,
        LevelOrStatistic = "Missing, n",
        Overall = values[[1L]], Survived = values[[2L]], Died = values[[3L]]
      )
    )
  }

  sample_sizes <- vapply(groups, nrow, integer(1))
  output <- rbind(
    output,
    data.table(
      Variable = "Cohort size",
      LevelOrStatistic = "n",
      Overall = as.character(sample_sizes[[1L]]),
      Survived = as.character(sample_sizes[[2L]]),
      Died = as.character(sample_sizes[[3L]])
    )
  )

  add_continuous("Age", "Age (years)")
  add_category("Gender", "Gender", 0, "Female")
  add_category("Gender", "Gender", 1, "Male")
  add_category_missing("Gender", "Gender")
  add_category("ICUType", "ICU type", 1, "Coronary care unit")
  add_category("ICUType", "ICU type", 2, "Cardiac surgery recovery unit")
  add_category("ICUType", "ICU type", 3, "Medical ICU")
  add_category("ICUType", "ICU type", 4, "Surgical ICU")
  add_category_missing("ICUType", "ICU type")
  add_continuous("Height", "Admission height (cm)")
  add_continuous("Weight", "Admission weight (kg)")
  add_continuous("SAPS_I", "SAPS-I score", 0L)
  add_continuous("SOFA", "SOFA score", 0L)
  output
}

make_outcome_summary <- function(cohort) {
  continuous_rows <- function(variable, label) {
    x <- cohort[[variable]]
    if (variable == "Survival") x[x == -1] <- NA_real_
    x_valid <- x[!is.na(x)]
    data.table(
      Outcome = label,
      Level = "Continuous",
      N = length(x_valid),
      Missing = sum(is.na(x)),
      Percent = NA_real_,
      Mean = mean(x_valid),
      SD = sd(x_valid),
      Median = median(x_valid),
      Q1 = quantile(x_valid, 0.25, names = FALSE, type = 2),
      Q3 = quantile(x_valid, 0.75, names = FALSE, type = 2),
      Min = min(x_valid),
      Max = max(x_valid)
    )
  }

  result <- rbindlist(list(
    continuous_rows("SAPS_I", "SAPS-I score"),
    continuous_rows("SOFA", "SOFA score"),
    continuous_rows("LengthOfStay", "Hospital length of stay (days)"),
    continuous_rows("Survival", "Survival time when death was recorded (days)")
  ), fill = TRUE)

  death_counts <- cohort[, .N, by = InHospitalDeath][order(InHospitalDeath)]
  death_counts[, `:=`(
    Outcome = "In-hospital death",
    Level = fifelse(InHospitalDeath == 1, "Died", "Survived"),
    Missing = sum(is.na(cohort$InHospitalDeath)),
    Percent = 100 * N / sum(N),
    Mean = NA_real_, SD = NA_real_, Median = NA_real_, Q1 = NA_real_,
    Q3 = NA_real_, Min = NA_real_, Max = NA_real_
  )]
  death_counts[, InHospitalDeath := NULL]
  setcolorder(death_counts, names(result))
  rbind(result, death_counts, fill = TRUE)
}

make_presence_table <- function(long_data, record_ids) {
  descriptor_only <- c("Age", "Gender", "Height", "ICUType")
  valid <- unique(
    long_data[!Parameter %chin% descriptor_only & !is.na(Value), .(RecordID, Parameter)]
  )
  valid[, Measured := 1L]
  presence <- dcast(valid, RecordID ~ Parameter, value.var = "Measured", fill = 0L)
  all_ids <- data.table(RecordID = record_ids)
  presence <- merge(all_ids, presence, by = "RecordID", all.x = TRUE)
  variable_columns <- setdiff(names(presence), "RecordID")
  for (column in variable_columns) set(presence, which(is.na(presence[[column]])), column, 0L)
  presence
}

make_missingness_outputs <- function(long_data, cohort, presence) {
  variables <- setdiff(names(presence), "RecordID")
  counts <- long_data[
    !Parameter %chin% c("Age", "Gender", "Height", "ICUType") & !is.na(Value),
    .(MeasurementCount = .N),
    by = .(RecordID, Parameter)
  ]

  coverage <- counts[, .(
    RecordsMeasured = uniqueN(RecordID),
    MedianMeasurementsAmongMeasured = as.numeric(median(MeasurementCount)),
    Q1MeasurementsAmongMeasured = as.numeric(
      quantile(MeasurementCount, 0.25, names = FALSE, type = 2)
    ),
    Q3MeasurementsAmongMeasured = as.numeric(
      quantile(MeasurementCount, 0.75, names = FALSE, type = 2)
    )
  ), by = Parameter]
  coverage[, `:=`(
    RecordsNeverMeasured = expected_n - RecordsMeasured,
    PercentRecordsMeasured = 100 * RecordsMeasured / expected_n
  )]
  setcolorder(
    coverage,
    c(
      "Parameter", "RecordsMeasured", "RecordsNeverMeasured", "PercentRecordsMeasured",
      "MedianMeasurementsAmongMeasured", "Q1MeasurementsAmongMeasured",
      "Q3MeasurementsAmongMeasured"
    )
  )
  setorder(coverage, PercentRecordsMeasured, Parameter)
  fwrite(coverage, file.path(out_dir, "variable_coverage.csv"), na = "")

  presence_with_outcome <- merge(
    presence,
    cohort[, .(RecordID, InHospitalDeath)],
    by = "RecordID",
    all.x = TRUE
  )
  presence_long <- melt(
    presence_with_outcome,
    id.vars = c("RecordID", "InHospitalDeath"),
    measure.vars = variables,
    variable.name = "Parameter",
    value.name = "Measured"
  )
  missingness_vs_death <- presence_long[, .(
    N = .N,
    Deaths = sum(InHospitalDeath == 1, na.rm = TRUE),
    DeathRate = mean(InHospitalDeath == 1, na.rm = TRUE)
  ), by = .(Parameter, Measured)]
  missingness_vs_death[, MeasurementStatus := fifelse(Measured == 1L, "Measured", "Never measured")]
  missingness_wide <- dcast(
    missingness_vs_death,
    Parameter ~ MeasurementStatus,
    value.var = c("N", "Deaths", "DeathRate")
  )
  if (all(c("DeathRate_Measured", "DeathRate_Never measured") %in% names(missingness_wide))) {
    missingness_wide[, DeathRateDifference := DeathRate_Measured - `DeathRate_Never measured`]
  }
  missingness_wide <- missingness_wide[order(-abs(DeathRateDifference))]
  fwrite(missingness_wide, file.path(out_dir, "missingness_vs_death.csv"), na = "")

  temporal <- unique(long_data[
    !Parameter %chin% c("Age", "Gender", "Height", "ICUType") &
      !is.na(Value) & TimeHours >= 0 & TimeHours <= 48,
    .(RecordID, Parameter, Hour = pmin(floor(TimeHours), 47L))
  ])[, .(RecordsMeasured = uniqueN(RecordID)), by = .(Parameter, Hour)]
  temporal <- merge(
    CJ(Parameter = variables, Hour = 0:47, unique = TRUE),
    temporal,
    by = c("Parameter", "Hour"),
    all.x = TRUE
  )
  temporal[is.na(RecordsMeasured), RecordsMeasured := 0L]
  temporal[, PercentRecordsMeasured := 100 * RecordsMeasured / expected_n]
  fwrite(temporal, file.path(out_dir, "temporal_measurement_coverage.csv"), na = "")

  variable_order <- coverage[order(PercentRecordsMeasured), Parameter]
  temporal[, Parameter := factor(Parameter, levels = variable_order)]
  temporal_plot <- ggplot(temporal, aes(x = Hour, y = Parameter, fill = PercentRecordsMeasured)) +
    geom_tile() +
    scale_fill_viridis_c(name = "% measured", limits = c(0, 100)) +
    scale_x_continuous(breaks = seq(0, 48, by = 6), expand = c(0, 0)) +
    labs(
      title = "Measurement availability during the first 48 ICU hours",
      x = "Hours after ICU admission",
      y = "Variable"
    ) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))
  ggsave(
    file.path(out_dir, "missingness_temporal.png"),
    temporal_plot,
    width = 11,
    height = 8,
    dpi = 180
  )

  plot_presence <- copy(presence_with_outcome)
  plot_presence[, VariablesMeasured := rowSums(.SD), .SDcols = variables]
  setorder(plot_presence, InHospitalDeath, VariablesMeasured, RecordID)
  plot_presence[, PatientOrder := seq_len(.N), by = InHospitalDeath]
  map_long <- melt(
    plot_presence,
    id.vars = c("RecordID", "InHospitalDeath", "PatientOrder"),
    measure.vars = variable_order,
    variable.name = "Parameter",
    value.name = "Measured"
  )
  map_long[, Parameter := factor(Parameter, levels = variable_order)]
  map_long[, DeathStatus := factor(
    InHospitalDeath,
    levels = c(0, 1),
    labels = c("Survived", "Died")
  )]
  missingness_plot <- ggplot(map_long, aes(x = Parameter, y = PatientOrder, fill = factor(Measured))) +
    geom_raster() +
    facet_grid(DeathStatus ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(
      values = c("0" = "#E5E7EB", "1" = "#2563EB"),
      labels = c("Not measured", "Measured"),
      name = NULL
    ) +
    labs(
      title = "Variable measurement map by ICU record",
      subtitle = "Records are ordered by death status and number of variables measured",
      x = "Variable",
      y = "ICU records"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )
  ggsave(
    file.path(out_dir, "missingness_map_wide.png"),
    missingness_plot,
    width = 12,
    height = 8,
    dpi = 180
  )

  coverage[, Parameter := factor(Parameter, levels = Parameter)]
  coverage_plot <- ggplot(coverage, aes(x = Parameter, y = PercentRecordsMeasured)) +
    geom_col(fill = "#2563EB") +
    coord_flip() +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(
      title = "Percentage of ICU records with each variable measured",
      x = "Variable",
      y = "Records measured (%)"
    ) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"))
  ggsave(
    file.path(out_dir, "variable_coverage.png"),
    coverage_plot,
    width = 9,
    height = 8,
    dpi = 180
  )
}

main <- function() {
  files <- check_input_data()
  message("Reading 4,000 set-a records...")
  raw_long <- rbindlist(lapply(files, read_one_record), use.names = TRUE)
  long_data <- clean_measurements(add_time_variables(raw_long))

  if (uniqueN(long_data$RecordID) != expected_n) {
    stop_with_message("The combined long table does not contain 4,000 records.")
  }
  if (any(long_data$TimeMinutes < 0 | long_data$TimeMinutes > 48 * 60, na.rm = TRUE)) {
    stop_with_message("A time value is outside the 0 to 48 hour window.")
  }

  fwrite(long_data, file.path(out_dir, "set-a_long.csv"), na = "")
  fwrite(make_quality_summary(long_data), file.path(out_dir, "measurement_range_check.csv"), na = "")

  outcomes <- fread(file.path(data_dir, "Outcomes-a.txt"), showProgress = FALSE)
  setnames(
    outcomes,
    old = c("SAPS-I", "Length_of_stay", "In-hospital_death"),
    new = c("SAPS_I", "LengthOfStay", "InHospitalDeath")
  )
  if (nrow(outcomes) != expected_n || uniqueN(outcomes$RecordID) != expected_n) {
    stop_with_message("Outcomes-a.txt does not contain 4,000 unique records.")
  }
  for (column in c("SAPS_I", "SOFA", "LengthOfStay")) {
    set(outcomes, which(outcomes[[column]] == -1), column, NA_real_)
  }

  static <- make_static_table(long_data)
  time_series_wide <- make_time_series_summary(long_data)
  cohort <- merge(static, outcomes, by = "RecordID", all = TRUE)
  wide <- merge(cohort, time_series_wide, by = "RecordID", all.x = TRUE)
  count_columns <- grep("__count$", names(wide), value = TRUE)
  for (column in count_columns) set(wide, which(is.na(wide[[column]])), column, 0)
  base_columns <- names(cohort)
  summary_variables <- sort(unique(sub("__.*$", "", setdiff(names(wide), base_columns))))
  statistic_order <- c("count", "first", "last", "min", "max", "mean")
  summary_columns <- unlist(
    lapply(summary_variables, function(variable) paste(variable, statistic_order, sep = "__")),
    use.names = FALSE
  )
  setcolorder(wide, c(base_columns, summary_columns))
  setorder(wide, RecordID)

  if (nrow(wide) != expected_n || uniqueN(wide$RecordID) != expected_n) {
    stop_with_message("The wide table does not contain one row for each of 4,000 records.")
  }
  fwrite(wide, file.path(out_dir, "set-a_wide.csv"), na = "")
  fwrite(make_table1(cohort), file.path(out_dir, "table1.csv"), na = "")
  fwrite(make_outcome_summary(cohort), file.path(out_dir, "outcome_summary.csv"), na = "")

  presence <- make_presence_table(long_data, cohort$RecordID)
  make_missingness_outputs(long_data, cohort, presence)

  validation <- data.table(
    Check = c(
      "Record files",
      "Unique records in long table",
      "Rows in wide table",
      "Unique records in wide table",
      "Outcome rows",
      "Time values outside 0-48 hours"
    ),
    Value = c(
      length(files),
      uniqueN(long_data$RecordID),
      nrow(wide),
      uniqueN(wide$RecordID),
      nrow(outcomes),
      sum(long_data$TimeMinutes < 0 | long_data$TimeMinutes > 48 * 60, na.rm = TRUE)
    ),
    Expected = c(expected_n, expected_n, expected_n, expected_n, expected_n, 0L)
  )
  validation[, Passed := Value == Expected]
  fwrite(validation, file.path(out_dir, "validation_checks.csv"), na = "")
  if (!all(validation$Passed)) stop_with_message("One or more validation checks failed.")

  message("Checkpoint 1 analysis completed. Files are in output/.")
}

main()
