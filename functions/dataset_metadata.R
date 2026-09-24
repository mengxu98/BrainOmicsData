normalize_missing_metadata <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  missing_values <- c(
    "", "NA", "N/A", "NULL", "None", "Unknown",
    "Unknown/not reported", "not reported", "not applicable"
  )
  x[
    is.na(x) |
      tolower(x) %in% tolower(missing_values)
  ] <- NA_character_
  x
}

read_geo_series_matrix_sample_metadata <- function(file) {
  if (!file.exists(file)) {
    stop("GEO series-matrix file is missing: ", file)
  }
  connection <- if (grepl("[.]gz$", file, ignore.case = TRUE)) {
    gzfile(file, "rt")
  } else {
    file(file, "rt")
  }
  lines <- tryCatch(
    readLines(connection, warn = FALSE),
    finally = close(connection)
  )
  sample_lines <- lines[startsWith(lines, "!Sample_")]
  if (length(sample_lines) == 0L) {
    stop("GEO series matrix contains no !Sample_ metadata rows")
  }
  parse_line <- function(line) {
    scan(
      text = line,
      what = character(),
      sep = "\t",
      quote = '"',
      quiet = TRUE,
      comment.char = ""
    )
  }
  parsed <- lapply(sample_lines, parse_line)
  keys <- vapply(parsed, `[[`, character(1), 1L)
  values <- lapply(parsed, function(value) value[-1L])
  accession_index <- which(keys == "!Sample_geo_accession")
  title_index <- which(keys == "!Sample_title")
  if (length(accession_index) != 1L || length(title_index) != 1L) {
    stop("GEO series matrix must contain one title and accession row")
  }
  n <- length(values[[accession_index]])
  if (n == 0L || any(lengths(values) != n)) {
    stop("GEO series-matrix sample metadata has inconsistent widths")
  }

  result <- data.frame(
    Sample_geo_accession = values[[accession_index]],
    Sample_title = values[[title_index]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ordinary <- which(
    !keys %in% c(
      "!Sample_geo_accession", "!Sample_title",
      "!Sample_characteristics_ch1"
    )
  )
  ordinary_occurrence <- integer(0)
  for (index in ordinary) {
    base_name <- sub("^!", "", keys[[index]])
    occurrence <- sum(keys[ordinary[ordinary <= index]] == keys[[index]])
    column <- if (occurrence == 1L) {
      base_name
    } else {
      paste0(base_name, "_", occurrence)
    }
    result[[column]] <- values[[index]]
  }

  characteristics <- which(keys == "!Sample_characteristics_ch1")
  for (index in characteristics) {
    characteristic <- values[[index]]
    characteristic_key <- trimws(sub(":.*$", "", characteristic))
    if (length(unique(characteristic_key)) != 1L) {
      stop("one GEO characteristics row contains multiple field names")
    }
    column <- paste0(
      "characteristic_",
      gsub("[^A-Za-z0-9]+", "_", tolower(characteristic_key[[1L]]))
    )
    if (column %in% names(result)) {
      stop("duplicated GEO characteristic field: ", column)
    }
    result[[column]] <- trimws(sub("^[^:]+:[[:space:]]*", "", characteristic))
  }
  if (anyNA(result$Sample_geo_accession) ||
    anyNA(result$Sample_title) ||
    anyDuplicated(result$Sample_geo_accession) ||
    anyDuplicated(result$Sample_title)) {
    stop("GEO series matrix contains incomplete or duplicated sample IDs")
  }
  result
}

bind_rows_by_name <- function(rows) {
  keep <- !vapply(rows, is.null, logical(1))
  rows <- rows[keep]
  if (length(rows) == 0L) {
    return(data.frame())
  }
  if (any(!vapply(rows, is.data.frame, logical(1)))) {
    stop("all row-binding inputs must be data frames")
  }

  column_order <- unique(unlist(
    lapply(rows, names),
    use.names = FALSE
  ))
  aligned <- lapply(rows, function(x) {
    missing_columns <- setdiff(column_order, names(x))
    for (column in missing_columns) {
      x[[column]] <- NA
    }
    x[, column_order, drop = FALSE]
  })
  result <- do.call(rbind, aligned)
  rownames(result) <- NULL
  result
}

metadata_count_values <- function(x, label = "count") {
  values <- if (is.numeric(x)) {
    as.numeric(x)
  } else {
    text <- trimws(as.character(x))
    invalid_text <- is.na(text) | !grepl("^[0-9]+$", text)
    if (any(invalid_text)) {
      stop(label, " contains missing or non-integer text values")
    }
    as.numeric(text)
  }
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0) ||
    any(values != floor(values))) {
    stop(label, " contains invalid non-negative whole-number counts")
  }
  values
}

formal_source_access_record <- function(dataset, source_access_file) {
  if (!file.exists(source_access_file)) {
    stop("source access registry is missing: ", source_access_file)
  }
  registry <- read.delim(
    source_access_file,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c(
    "dataset", "source_accession", "source_repository",
    "publication_doi", "verified_title", "journal",
    "publication_year", "source_url", "repository_record_url"
  )
  missing <- setdiff(required, names(registry))
  if (length(missing) > 0L) {
    stop(
      "source access registry is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  selected <- registry[registry$dataset == dataset, , drop = FALSE]
  if (nrow(selected) != 1L) {
    stop(dataset, " must have exactly one source access registry row")
  }
  for (column in c(
    "source_accession", "source_repository", "publication_doi",
    "verified_title", "journal", "publication_year",
    "repository_record_url"
  )) {
    value <- normalize_missing_metadata(selected[[column]])
    if (length(value) != 1L || is.na(value)) {
      stop(dataset, " source registry lacks ", column)
    }
  }
  if (!grepl("^10[.]", selected$publication_doi[[1L]])) {
    stop(dataset, " source registry has an invalid publication DOI")
  }
  selected
}

add_source_publication_metadata <- function(meta, source_record) {
  if (!is.data.frame(meta) || !is.data.frame(source_record) ||
    nrow(source_record) != 1L) {
    stop("source publication metadata requires one validated registry row")
  }
  field_map <- c(
    Source_Accession = "source_accession",
    Source_Repository = "source_repository",
    Source_Publication_DOI = "publication_doi",
    Source_Publication_Title = "verified_title",
    Source_Publication_Journal = "journal",
    Source_Publication_Year = "publication_year",
    Source_Repository_Record_URL = "repository_record_url"
  )
  missing <- setdiff(unname(field_map), names(source_record))
  if (length(missing) > 0L) {
    stop("source registry row lacks: ", paste(missing, collapse = ", "))
  }
  for (target in names(field_map)) {
    value <- normalize_missing_metadata(source_record[[field_map[[target]]]])
    if (length(value) != 1L || is.na(value)) {
      stop("source registry row has no value for ", field_map[[target]])
    }
    meta[[target]] <- rep(value, nrow(meta))
  }
  meta
}

restore_deterministic_metadata_provenance <- function(
  metadata,
  previous_metadata
) {
  if ("Source_Cell_Index" %in% names(previous_metadata)) {
    source_cell_index <- suppressWarnings(as.integer(
      previous_metadata$Source_Cell_Index
    ))
    expected <- seq_len(nrow(previous_metadata))
    if (anyNA(source_cell_index) ||
      !identical(source_cell_index, expected) ||
      nrow(metadata) != nrow(previous_metadata)) {
      stop(
        "existing Source_Cell_Index is not the complete deterministic ",
        "source-cell order"
      )
    }
    metadata$Source_Cell_Index <- source_cell_index
  }
  metadata
}

parse_dataset_age <- function(age, source_basis = NULL) {
  raw <- normalize_missing_metadata(age)
  if (is.null(source_basis)) {
    source_basis <- rep(NA_character_, length(raw))
  } else if (length(source_basis) == 1L) {
    source_basis <- rep(source_basis, length(raw))
  } else if (length(source_basis) != length(raw)) {
    stop("source_basis must have length 1 or match age")
  }
  source_basis <- normalize_missing_metadata(source_basis)

  normalize_source_basis <- function(value) {
    if (is.na(value)) {
      return(NA_character_)
    }
    text <- tolower(trimws(value))
    if (grepl("gestational|last menstrual|^gw$|^ga$", text)) {
      return("gestational")
    }
    if (grepl(
      "post[- ]?(conception|fertilization)|^pcw$",
      text,
      perl = TRUE
    )) {
      return("postconceptional")
    }
    if (grepl("postnatal|age at death|chronological", text)) {
      return("postnatal")
    }
    stop("unsupported age source basis: ", value)
  }
  source_basis <- unname(vapply(
    source_basis,
    normalize_source_basis,
    character(1)
  ))

  age_parse_result <- function(value = NA_real_,
                               lower = value,
                               upper = value,
                               unit = NA_character_,
                               explicit = FALSE,
                               method = "unparsed",
                               formula = NA_character_,
                               confidence = "unknown",
                               source_unit = NA_character_,
                               source_basis = NA_character_,
                               conversion_applied = FALSE) {
    list(
      age_value = value,
      age_lower = lower,
      age_upper = upper,
      age_unit = unit,
      age_parse_explicit = explicit,
      age_representative_method = method,
      age_parser_formula = formula,
      age_parser_confidence = confidence,
      age_source_unit = source_unit,
      age_source_basis = source_basis,
      age_conversion_applied = conversion_applied
    )
  }

  first_numeric_match <- function(text, pattern, group = 1L) {
    matched <- regexec(pattern, text, perl = TRUE)
    pieces <- regmatches(text, matched)[[1L]]
    if (length(pieces) <= group) {
      return(NA_real_)
    }
    suppressWarnings(as.numeric(pieces[[group + 1L]]))
  }

  numeric_range_match <- function(text, pattern) {
    matched <- regexec(pattern, text, perl = TRUE)
    pieces <- regmatches(text, matched)[[1L]]
    if (length(pieces) != 3L) {
      return(c(NA_real_, NA_real_))
    }
    suppressWarnings(as.numeric(pieces[2:3]))
  }

  week_basis_from_text <- function(text, supplied_basis) {
    if (grepl(
      "pcw|post[- ]?(conception|fertilization)",
      text,
      perl = TRUE
    )) {
      return("postconceptional")
    }
    if (grepl(
      paste0(
        "gestational|",
        "(?:^|[^[:alpha:]])(?:gw|ga)",
        "(?=[0-9[:space:]+-]|$)"
      ),
      text,
      perl = TRUE
    )) {
      return("gestational")
    }
    supplied_basis
  }

  week_basis_is_explicit <- function(text) {
    grepl(
      paste0(
        "pcw|post[- ]?(conception|fertilization)|gestational|",
        "(?:^|[^[:alpha:]])(?:gw|ga)",
        "(?=[0-9[:space:]+-]|$)"
      ),
      text,
      perl = TRUE
    )
  }

  week_source_unit <- function(basis) {
    if (identical(basis, "gestational")) {
      return("gestational weeks")
    }
    if (identical(basis, "postconceptional")) {
      return("postconceptional weeks")
    }
    "weeks (basis unspecified)"
  }

  canonicalize_week_values <- function(values, basis) {
    if (is.na(basis)) {
      return(rep(NA_real_, length(values)))
    }
    if (identical(basis, "gestational")) {
      values <- values - 2
    }
    values[values < 0] <- NA_real_
    values
  }

  parse_one <- function(value, supplied_basis) {
    if (is.na(value)) {
      return(age_parse_result())
    }
    text <- tolower(trimws(value))

    if (grepl("newborn|neonatal", text) &&
      !grepl("[0-9]", text)) {
      return(age_parse_result(
        value = 0,
        unit = "years",
        explicit = TRUE,
        method = "lexical birth anchor",
        formula = "newborn/neonatal = 0 years",
        confidence = "medium",
        source_unit = "lexical developmental stage",
        source_basis = "postnatal"
      ))
    }

    week_range <- numeric_range_match(
      text,
      paste0(
        "(?:\\b(?:gw|ga|pcw)\\s*)?",
        "([0-9]+(?:\\.[0-9]+)?)\\s*(?:-|to)\\s*",
        "([0-9]+(?:\\.[0-9]+)?)\\s*",
        "(?:gw|ga|pcw|gestational\\s+weeks?|fetal\\s+weeks?|",
        "weeks?\\s+post[- ]?(?:conception|fertilization)|weeks?|w)\\b"
      )
    )
    if (anyNA(week_range)) {
      week_range <- numeric_range_match(
        text,
        paste0(
          "\\b(?:gw|ga|pcw)\\s*",
          "([0-9]+(?:\\.[0-9]+)?)\\s*(?:-|to)\\s*",
          "([0-9]+(?:\\.[0-9]+)?)\\b"
        )
      )
    }
    if (!anyNA(week_range)) {
      week_basis <- week_basis_from_text(text, supplied_basis)
      canonical_range <- canonicalize_week_values(
        sort(week_range),
        week_basis
      )
      if (anyNA(canonical_range)) {
        return(age_parse_result(
          explicit = FALSE,
          method = "unparsed week range with unspecified basis",
          confidence = "requires source-specific week basis",
          source_unit = week_source_unit(week_basis),
          source_basis = week_basis
        ))
      }
      return(age_parse_result(
        value = mean(canonical_range),
        lower = canonical_range[[1L]],
        upper = canonical_range[[2L]],
        unit = "PCW",
        explicit = week_basis_is_explicit(text),
        method = "midpoint of reported range",
        formula = if (identical(week_basis, "gestational")) {
          "reported gestational-week bounds - 2 weeks = PCW; midpoint used"
        } else {
          "(reported lower PCW + reported upper PCW) / 2"
        },
        confidence = if (week_basis_is_explicit(text)) {
          "reported range with explicit week basis"
        } else {
          "reported range with dataset-supplied week basis"
        },
        source_unit = week_source_unit(week_basis),
        source_basis = week_basis,
        conversion_applied = identical(week_basis, "gestational")
      ))
    }

    prefixed_week <- first_numeric_match(
      text,
      "\\b(?:gw|ga|pcw)\\s*([0-9]+(?:\\.[0-9]+)?)",
      group = 1L
    )
    suffixed_week <- first_numeric_match(
      text,
      paste0(
        "([0-9]+(?:\\.[0-9]+)?)(?:st|nd|rd|th)?\\s*",
        "(?:gw|ga|pcw|gestational\\s+weeks?|fetal\\s+weeks?|",
        "weeks?\\s+post[- ]?(?:conception|fertilization)|weeks?|w)"
      ),
      group = 1L
    )
    week_value <- if (!is.na(prefixed_week)) {
      prefixed_week
    } else {
      suffixed_week
    }
    if (!is.na(week_value)) {
      week_basis <- week_basis_from_text(text, supplied_basis)
      extra_days <- first_numeric_match(
        text,
        "\\+\\s*([0-9]+(?:\\.[0-9]+)?)\\s*(?:d|days?)\\b",
        group = 1L
      )
      if (!is.na(extra_days)) {
        week_value <- week_value + extra_days / 7
      }
      canonical_week <- canonicalize_week_values(
        week_value,
        week_basis
      )
      if (is.na(canonical_week)) {
        return(age_parse_result(
          explicit = FALSE,
          method = "unparsed week value with unspecified basis",
          confidence = "requires source-specific week basis",
          source_unit = week_source_unit(week_basis),
          source_basis = week_basis
        ))
      }
      return(age_parse_result(
        value = canonical_week,
        unit = "PCW",
        explicit = week_basis_is_explicit(text),
        method = "reported point age",
        formula = if (identical(week_basis, "gestational")) {
          if (is.na(extra_days)) {
            "reported gestational weeks - 2 = PCW"
          } else {
            paste(
              "reported gestational weeks + reported days / 7",
              "- 2 = PCW"
            )
          }
        } else if (is.na(extra_days)) {
          "reported PCW"
        } else {
          "reported PCW + reported days / 7"
        },
        confidence = if (week_basis_is_explicit(text)) {
          "high"
        } else {
          "dataset-supplied week basis"
        },
        source_unit = week_source_unit(week_basis),
        source_basis = week_basis,
        conversion_applied = identical(week_basis, "gestational")
      ))
    }

    composite_match <- regexec(
      paste0(
        "([0-9]+(?:\\.[0-9]+)?)\\s*years?\\s+",
        "([0-9]+(?:\\.[0-9]+)?)\\s*days?"
      ),
      text,
      perl = TRUE
    )
    composite <- regmatches(text, composite_match)[[1L]]
    if (length(composite) == 3L) {
      return(age_parse_result(
        value =
          as.numeric(composite[[2L]]) +
            as.numeric(composite[[3L]]) / 365,
        unit = "years",
        explicit = TRUE,
        method = "reported point age",
        formula = "reported years + reported days / 365",
        confidence = "high",
        source_unit = "years and days",
        source_basis = "postnatal"
      ))
    }

    month_range <- numeric_range_match(
      text,
      paste0(
        "([0-9]+(?:\\.[0-9]+)?)\\s*(?:-|to)\\s*",
        "([0-9]+(?:\\.[0-9]+)?)\\s*months?\\b"
      )
    )
    if (!anyNA(month_range)) {
      month_range <- sort(month_range) / 12
      return(age_parse_result(
        value = mean(month_range),
        lower = month_range[[1L]],
        upper = month_range[[2L]],
        unit = "years",
        explicit = TRUE,
        method = "midpoint of reported range",
        formula = "(reported lower months + reported upper months) / 24",
        confidence = "reported range",
        source_unit = "months",
        source_basis = "postnatal",
        conversion_applied = TRUE
      ))
    }

    month_value <- first_numeric_match(
      text,
      "([0-9]+(?:\\.[0-9]+)?)[- ]*months?\\b",
      group = 1L
    )
    if (!is.na(month_value)) {
      return(age_parse_result(
        value = month_value / 12,
        unit = "years",
        explicit = TRUE,
        method = "reported point age",
        formula = "reported months / 12",
        confidence = "high",
        source_unit = "months",
        source_basis = "postnatal",
        conversion_applied = TRUE
      ))
    }

    year_range <- numeric_range_match(
      text,
      paste0(
        "([0-9]+(?:\\.[0-9]+)?)\\s*(?:-|to)\\s*",
        "([0-9]+(?:\\.[0-9]+)?)\\s*(?:years?|yrs?|yr)\\b"
      )
    )
    if (!anyNA(year_range)) {
      year_range <- sort(year_range)
      return(age_parse_result(
        value = mean(year_range),
        lower = year_range[[1L]],
        upper = year_range[[2L]],
        unit = "years",
        explicit = TRUE,
        method = "midpoint of reported range",
        formula = "(reported lower years + reported upper years) / 2",
        confidence = "reported range",
        source_unit = "years",
        source_basis = "postnatal"
      ))
    }

    year_open_lower <- first_numeric_match(
      text,
      "^([0-9]+(?:\\.[0-9]+)?)\\s*\\+\\s*(?:years?|yrs?|yr)$",
      group = 1L
    )
    if (!is.na(year_open_lower)) {
      return(age_parse_result(
        value = year_open_lower,
        lower = year_open_lower,
        upper = Inf,
        unit = "years",
        explicit = TRUE,
        method = "reported open-ended lower bound",
        formula = paste(
          "reported lower bound retained; representative value equals",
          "the lower bound"
        ),
        confidence = "open-ended reported range",
        source_unit = "years",
        source_basis = "postnatal"
      ))
    }

    year_value <- first_numeric_match(
      text,
      "([0-9]+(?:\\.[0-9]+)?)[- ]*(?:years?|yrs?|yr)\\b",
      group = 1L
    )
    if (!is.na(year_value)) {
      return(age_parse_result(
        value = year_value,
        unit = "years",
        explicit = TRUE,
        method = "reported point age",
        formula = "reported years",
        confidence = "high",
        source_unit = "years",
        source_basis = "postnatal"
      ))
    }

    day_range <- numeric_range_match(
      text,
      paste0(
        "([0-9]+(?:\\.[0-9]+)?)\\s*(?:-|to)\\s*",
        "([0-9]+(?:\\.[0-9]+)?)\\s*(?:days?|d)\\b"
      )
    )
    if (!anyNA(day_range)) {
      day_range <- sort(day_range) / 365
      return(age_parse_result(
        value = mean(day_range),
        lower = day_range[[1L]],
        upper = day_range[[2L]],
        unit = "years",
        explicit = TRUE,
        method = "midpoint of reported range",
        formula = "(reported lower days + reported upper days) / 730",
        confidence = "reported range",
        source_unit = "days",
        source_basis = "postnatal",
        conversion_applied = TRUE
      ))
    }

    day_value <- first_numeric_match(
      text,
      "([0-9]+(?:\\.[0-9]+)?)[- ]*(?:days?|d)\\b",
      group = 1L
    )
    if (!is.na(day_value)) {
      return(age_parse_result(
        value = day_value / 365,
        unit = "years",
        explicit = TRUE,
        method = "reported point age",
        formula = "reported days / 365",
        confidence = "high",
        source_unit = "days",
        source_basis = "postnatal",
        conversion_applied = TRUE
      ))
    }

    if (grepl("^[0-9]+(?:\\.[0-9]+)?$", text, perl = TRUE)) {
      return(age_parse_result(
        value = as.numeric(text),
        unit = "years",
        explicit = FALSE,
        method = "numeric value with inferred unit",
        formula = "bare numeric value interpreted as years",
        confidence = "unit inferred",
        source_unit = "unitless numeric",
        source_basis = "postnatal"
      ))
    }

    age_parse_result()
  }

  unique_input <- unique(data.frame(
    age = raw,
    source_basis = source_basis,
    stringsAsFactors = FALSE
  ))
  parsed_unique <- Map(
    parse_one,
    unique_input$age,
    unique_input$source_basis
  )
  input_key <- paste(raw, source_basis, sep = "\r")
  unique_key <- paste(
    unique_input$age,
    unique_input$source_basis,
    sep = "\r"
  )
  idx <- match(input_key, unique_key)
  extract_parsed <- function(name, type) {
    unname(vapply(parsed_unique, `[[`, type, name)[idx])
  }

  data.frame(
    age_raw = raw,
    age_value = extract_parsed("age_value", numeric(1)),
    age_lower = extract_parsed("age_lower", numeric(1)),
    age_upper = extract_parsed("age_upper", numeric(1)),
    age_unit = extract_parsed("age_unit", character(1)),
    age_parse_explicit = extract_parsed(
      "age_parse_explicit",
      logical(1)
    ),
    age_representative_method = extract_parsed(
      "age_representative_method",
      character(1)
    ),
    age_parser_formula = extract_parsed(
      "age_parser_formula",
      character(1)
    ),
    age_parser_confidence = extract_parsed(
      "age_parser_confidence",
      character(1)
    ),
    age_source_unit = extract_parsed(
      "age_source_unit",
      character(1)
    ),
    age_source_basis = extract_parsed(
      "age_source_basis",
      character(1)
    ),
    age_conversion_applied = extract_parsed(
      "age_conversion_applied",
      logical(1)
    ),
    stringsAsFactors = FALSE
  )
}

dataset_age_interval <- function(age_value, age_unit) {
  stage <- rep(NA_character_, length(age_value))
  intervals <- list(
    c("PCW", 4, 8, "S1"),
    c("PCW", 8, 10, "S2"),
    c("PCW", 10, 13, "S3"),
    c("PCW", 13, 16, "S4"),
    c("PCW", 16, 19, "S5"),
    c("PCW", 19, 24, "S6"),
    c("PCW", 24, 40, "S7"),
    c("years", 0, 0.5, "S8"),
    c("years", 0.5, 1, "S9"),
    c("years", 1, 6, "S10"),
    c("years", 6, 12, "S11"),
    c("years", 12, 20, "S12"),
    c("years", 20, 40, "S13"),
    c("years", 40, 60, "S14"),
    c("years", 60, Inf, "S15")
  )
  for (interval in intervals) {
    selected <- !is.na(age_value) &
      !is.na(age_unit) &
      age_unit == interval[[1L]] &
      age_value >= as.numeric(interval[[2L]]) &
      age_value < as.numeric(interval[[3L]])
    stage[selected] <- interval[[4L]]
  }
  stage
}

apply_curated_age_interval_decisions <- function(meta) {
  if (!is.data.frame(meta) || !"Dataset" %in% names(meta)) {
    return(meta)
  }
  n <- nrow(meta)
  if (!"Stage" %in% names(meta)) {
    meta$Stage <- rep(NA_character_, n)
  }
  if (!"Age_Interval_Assignment_Status" %in% names(meta)) {
    meta$Age_Interval_Assignment_Status <- rep(NA_character_, n)
  }
  if (!"Age_Interval_Curation_Evidence" %in% names(meta)) {
    meta$Age_Interval_Curation_Evidence <- rep(NA_character_, n)
  }
  if (!"Age_Interval_Curation_Confidence" %in% names(meta)) {
    meta$Age_Interval_Curation_Confidence <- rep(NA_character_, n)
  }
  if (!"Continuous_Age_Eligibility" %in% names(meta)) {
    meta$Continuous_Age_Eligibility <- rep(NA_character_, n)
  }

  age_input <- rep(NA_character_, n)
  if ("Age_Harmonization_Input" %in% names(meta)) {
    age_input <- normalize_missing_metadata(meta$Age_Harmonization_Input)
  }
  if ("Age" %in% names(meta)) {
    fallback <- normalize_missing_metadata(meta$Age)
    age_input[is.na(age_input)] <- fallback[is.na(age_input)]
  }

  unresolved_stage <- is.na(normalize_missing_metadata(meta$Stage))
  gse67835_age_match <- !is.na(age_input) & grepl(
    "^[[:space:]]*14[[:space:]]*-[[:space:]]*16[[:space:]]*PCW[[:space:]]*$",
    age_input,
    ignore.case = TRUE
  )
  gse67835_slim_match <- unresolved_stage
  if ("BrainRegion" %in% names(meta)) {
    gse67835_slim_match <- gse67835_slim_match &
      meta$BrainRegion == "Cerebral cortex"
  }
  gse67835_prenatal <- meta$Dataset == "GSE67835" &
    (gse67835_age_match | gse67835_slim_match)
  meta$Stage[gse67835_prenatal] <- "S4"
  meta$Age_Interval_Assignment_Status[gse67835_prenatal] <- paste(
    "curated study-level range assignment: the 14-16 PCW prenatal",
    "GSE67835 specimens are represented by the range midpoint (15 PCW)",
    "and assigned to S4"
  )
  meta$Age_Interval_Curation_Evidence[gse67835_prenatal] <- paste(
    "Darmanis et al. 2015 reports 16-18 gestational weeks; conversion",
    "to 14-16 postconception weeks is retained, and the documented",
    "descriptive interval decision assigns the cohort to S4"
  )
  meta$Age_Interval_Curation_Confidence[gse67835_prenatal] <-
    "moderate: source range is published; interval uses a group-level midpoint"
  meta$Continuous_Age_Eligibility[gse67835_prenatal] <-
    "not eligible: exact donor age is not available"

  gse144136 <- meta$Dataset == "GSE144136" &
    is.na(normalize_missing_metadata(meta$Stage))
  meta$Stage[gse144136] <- "S13"

  if (!"AgeIntervalID" %in% names(meta)) {
    meta$AgeIntervalID <- rep(NA_character_, n)
  }
  curated_stage <- meta$Dataset %in% c("GSE67835", "GSE144136") &
    !is.na(normalize_missing_metadata(meta$Stage))
  meta$AgeIntervalID[curated_stage] <- meta$Stage[curated_stage]

  meta
}

add_age_schema <- function(meta) {
  if (!"Age" %in% names(meta)) {
    return(meta)
  }
  age_input <- normalize_missing_metadata(meta$Age)
  use_harmonized_input <- rep(FALSE, nrow(meta))
  if ("Age_Harmonization_Input" %in% names(meta)) {
    harmonized_input <- normalize_missing_metadata(
      meta$Age_Harmonization_Input
    )
    use_harmonized_input <- !is.na(harmonized_input)
    age_input[use_harmonized_input] <-
      harmonized_input[use_harmonized_input]
  }
  # Restore the published control-group mean in older processed metadata.
  # This is a descriptive cohort representative, never an individual age.
  if (!"Age_Source_Raw" %in% names(meta)) {
    meta$Age_Source_Raw <- normalize_missing_metadata(meta$Age)
  }
  if ("Dataset" %in% names(meta)) {
    legacy_control <- meta$Dataset %in% "GSE144136" & age_input %in% "38 years"
    age_input[legacy_control] <- "38.71 years"
    meta$Age[legacy_control] <- "38.71 years"
    if ("Age_Harmonization_Input" %in% names(meta)) {
      meta$Age_Harmonization_Input[legacy_control] <- "38.71 years"
    }
  }
  source_basis <- rep(NA_character_, nrow(meta))
  if ("Age_Source_Basis" %in% names(meta)) {
    use_source_basis <- !use_harmonized_input & !is.na(age_input)
    source_basis[use_source_basis] <-
      meta$Age_Source_Basis[use_source_basis]
  }
  parsed <- parse_dataset_age(
    age_input,
    source_basis = source_basis
  )

  if (!"Age_Source_Raw" %in% names(meta)) {
    meta$Age_Source_Raw <- normalize_missing_metadata(meta$Age)
  }
  if (!"Age_Harmonization_Input" %in% names(meta)) {
    meta$Age_Harmonization_Input <- age_input
  } else {
    meta$Age_Harmonization_Input[!use_harmonized_input] <-
      age_input[!use_harmonized_input]
  }
  if (!"Age_Source_Unit" %in% names(meta)) {
    meta$Age_Source_Unit <- parsed$age_source_unit
  }
  if (!"Age_Source_Basis" %in% names(meta)) {
    meta$Age_Source_Basis <- parsed$age_source_basis
  }
  if (!"Age_Conversion_Formula" %in% names(meta)) {
    meta$Age_Conversion_Formula <- parsed$age_parser_formula
  }
  if (!"Age_Conversion_Confidence" %in% names(meta)) {
    meta$Age_Conversion_Confidence <- parsed$age_parser_confidence
  }
  if (!"Age_Conversion_Applied" %in% names(meta)) {
    meta$Age_Conversion_Applied <- parsed$age_conversion_applied
  }
  meta$Unit <- ifelse(parsed$age_unit == "PCW", "PCW", "Years")
  meta$Unit[is.na(parsed$age_unit)] <- NA_character_
  meta$Age_num <- parsed$age_value
  meta$Age_Lower <- parsed$age_lower
  meta$Age_Upper <- parsed$age_upper
  meta$Age_Canonical_Value <- parsed$age_value
  meta$Age_Canonical_Lower <- parsed$age_lower
  meta$Age_Canonical_Upper <- parsed$age_upper
  meta$Age_Canonical_Unit <- parsed$age_unit
  meta$Age_Parse_Explicit <- parsed$age_parse_explicit
  meta$Age_Representative_Method <-
    parsed$age_representative_method
  meta$Age_Parser_Formula <- parsed$age_parser_formula
  meta$Age_Parser_Confidence <- parsed$age_parser_confidence
  representative_stage <- dataset_age_interval(
    age_value = parsed$age_value,
    age_unit = parsed$age_unit
  )
  lower_stage <- dataset_age_interval(
    age_value = parsed$age_lower,
    age_unit = parsed$age_unit
  )
  upper_stage <- dataset_age_interval(
    age_value = parsed$age_upper,
    age_unit = parsed$age_unit
  )
  is_range <- !is.na(parsed$age_lower) &
    !is.na(parsed$age_upper) &
    parsed$age_lower != parsed$age_upper
  final_open_ended_range <- is_range &
    is.infinite(parsed$age_upper) &
    lower_stage == "S15"
  range_with_one_interval <- is_range &
    !is.na(lower_stage) &
    (lower_stage == upper_stage | final_open_ended_range)
  range_crosses_interval <- is_range & !range_with_one_interval
  meta$Stage <- representative_stage
  meta$Stage[range_with_one_interval] <-
    lower_stage[range_with_one_interval]
  meta$Stage[range_crosses_interval] <- NA_character_
  meta$Age_Interval_Assignment_Status <- ifelse(
    is.na(parsed$age_value),
    "unparsed",
    ifelse(
      range_crosses_interval,
      "ambiguous reported range crosses interval boundary",
      ifelse(
        is.na(meta$Stage),
        "outside defined intervals",
        ifelse(
          is_range,
          "reported range contained within one interval",
          "representative point assigned"
        )
      )
    )
  )
  meta$Age_Interval_Boundary_Convention <-
    "lower inclusive; upper exclusive"
  meta <- apply_curated_age_interval_decisions(meta)
  if ("Dataset" %in% names(meta)) {
    group_mean <- meta$Dataset %in% "GSE144136" &
      age_input %in% c("38.71 years", "41.06 years")
    meta$Age_Representative_Method[group_mean] <- "published diagnosis-group mean"
    meta$Age_Parser_Formula[group_mean] <- "published group mean retained without donor imputation"
    meta$Age_Parser_Confidence[group_mean] <- "low: group-level representative"
    meta$Age_Conversion_Formula[group_mean] <- "published diagnosis-group mean used only as a descriptive stage representative"
    meta$Age_Conversion_Confidence[group_mean] <- "low: group mean; exact donor ages unavailable in public sources"
    meta$Age_Interval_Curation_Confidence[group_mean] <- "low: descriptive stage based on a published group mean"
    meta$Continuous_Age_Eligibility[group_mean] <- "not eligible: exact donor age is not available"
  }
  meta
}

format_canonical_age_number <- function(value) {
  format(
    value,
    scientific = FALSE,
    trim = TRUE,
    digits = 10
  )
}

canonical_age_text_from_parsed <- function(parsed) {
  n <- length(parsed$age_value)
  result <- rep(NA_character_, n)
  point <- !is.na(parsed$age_value) &
    !is.na(parsed$age_lower) &
    !is.na(parsed$age_upper) &
    parsed$age_lower == parsed$age_upper
  ranged <- !is.na(parsed$age_value) & !point
  result[point] <- paste(
    format_canonical_age_number(parsed$age_value[point]),
    parsed$age_unit[point]
  )
  result[ranged] <- paste0(
    format_canonical_age_number(parsed$age_lower[ranged]),
    "-",
    format_canonical_age_number(parsed$age_upper[ranged]),
    " ",
    parsed$age_unit[ranged]
  )
  result
}

standardize_source_age_metadata <- function(
  meta,
  source_age,
  source_basis = NULL,
  source_unit = NULL,
  source_reference = NA_character_
) {
  if (!is.data.frame(meta) || length(source_age) != nrow(meta)) {
    stop("source age must contain one value per metadata row")
  }
  if (is.null(source_basis)) {
    source_basis <- rep(NA_character_, nrow(meta))
  } else if (length(source_basis) == 1L) {
    source_basis <- rep(source_basis, nrow(meta))
  }
  if (length(source_basis) != nrow(meta)) {
    stop("source age basis must contain one value per metadata row")
  }
  parsed <- parse_dataset_age(source_age, source_basis = source_basis)
  canonical_age <- canonical_age_text_from_parsed(parsed)
  if (is.null(source_unit)) {
    source_unit <- parsed$age_source_unit
  } else if (length(source_unit) == 1L) {
    source_unit <- rep(source_unit, nrow(meta))
  }
  if (length(source_unit) != nrow(meta)) {
    stop("source age unit must contain one value per metadata row")
  }
  if (length(source_reference) == 1L) {
    source_reference <- rep(source_reference, nrow(meta))
  }
  if (length(source_reference) != nrow(meta)) {
    stop("source age reference must contain one value per metadata row")
  }

  parsed_value <- !is.na(parsed$age_value)
  source_unit_normalized <- normalize_missing_metadata(source_unit)
  canonical_unit_normalized <- normalize_missing_metadata(parsed$age_unit)
  unit_changed <- parsed_value &
    !is.na(source_unit_normalized) &
    source_unit_normalized != canonical_unit_normalized

  meta$Age_Source_Raw <- normalize_missing_metadata(source_age)
  meta$Age_Source_Unit <- source_unit_normalized
  meta$Age_Source_Basis <- parsed$age_source_basis
  meta$Age_Source_Reference <- normalize_missing_metadata(source_reference)
  meta$Age_Harmonization_Input <- canonical_age
  meta$Age_Conversion_Formula <- parsed$age_parser_formula
  meta$Age_Conversion_Confidence <- parsed$age_parser_confidence
  meta$Age_Conversion_Applied <- parsed_value &
    (parsed$age_conversion_applied | unit_changed)
  meta$Age <- canonical_age
  meta
}

apply_existing_dataset_age_rules <- function(meta, dataset) {
  if (!identical(dataset, "Nowakowski_et_al_2017")) {
    return(meta)
  }
  provenance_fields <- c(
    "Age_Source_Raw", "Age_Source_Basis", "Age_Harmonization_Input",
    "Age_Conversion_Formula", "Age_Conversion_Confidence",
    "Age_Conversion_Applied"
  )
  if (all(provenance_fields %in% names(meta)) &&
    any(!is.na(normalize_missing_metadata(
      meta$Age_Harmonization_Input
    )))) {
    return(meta)
  }
  if (!"Age" %in% names(meta)) {
    stop("Nowakowski metadata is missing the retained Age field")
  }

  legacy_age <- normalize_missing_metadata(meta$Age)
  recoverable <- !is.na(legacy_age) & grepl(
    "\\bPCW\\s*$",
    legacy_age,
    ignore.case = TRUE,
    perl = TRUE
  )
  recovered_source_age <- legacy_age
  recovered_source_age[recoverable] <- sub(
    "\\s*PCW\\s*$",
    "w",
    legacy_age[recoverable],
    ignore.case = TRUE,
    perl = TRUE
  )
  parsed <- parse_dataset_age(
    recovered_source_age,
    source_basis = ifelse(recoverable, "gestational", NA_character_)
  )
  converted <- recoverable &
    !is.na(parsed$age_value) &
    parsed$age_unit == "PCW"

  canonical_age <- rep(NA_character_, nrow(meta))
  point_age <- converted &
    !is.na(parsed$age_lower) &
    !is.na(parsed$age_upper) &
    parsed$age_lower == parsed$age_upper
  range_age <- converted & !point_age
  canonical_age[point_age] <- paste(
    format_canonical_age_number(parsed$age_value[point_age]),
    "PCW"
  )
  canonical_age[range_age] <- paste0(
    format_canonical_age_number(parsed$age_lower[range_age]),
    "-",
    format_canonical_age_number(parsed$age_upper[range_age]),
    " PCW"
  )

  meta$Age_Legacy_Processed <- legacy_age
  meta$Age_Source_Raw <- ifelse(
    recoverable,
    recovered_source_age,
    legacy_age
  )
  meta$Age_Source_Unit <- ifelse(
    recoverable,
    "gestational weeks",
    NA_character_
  )
  meta$Age_Source_Reference <- ifelse(
    recoverable,
    paste(
      "Nowakowski et al. 2017, DOI 10.1126/science.aap8809;",
      "legacy donor_age gestational-week label recovered"
    ),
    NA_character_
  )
  meta$Age_Source_Basis <- ifelse(
    recoverable,
    "gestational age",
    NA_character_
  )
  meta$Age_Source_Recovery_Method <- ifelse(
    recoverable,
    paste(
      "reversed the legacy deterministic substitution",
      "gsub('w', ' PCW', donor_age)"
    ),
    "source age not recoverable from retained metadata"
  )
  meta$Age_Harmonization_Input <- ifelse(
    converted,
    canonical_age,
    legacy_age
  )
  meta$Age_Conversion_Formula <- ifelse(
    converted,
    parsed$age_parser_formula,
    NA_character_
  )
  meta$Age_Conversion_Confidence <- ifelse(
    converted,
    paste(
      "high: source gestational-week semantics and reversible",
      "legacy text substitution"
    ),
    "unresolved"
  )
  meta$Age_Conversion_Applied <- converted
  meta$Age[converted] <- canonical_age[converted]
  meta
}

metadata_source_value <- function(meta, source, default = NA_character_) {
  if (is.null(source) || length(source) == 0L || all(is.na(source))) {
    return(rep(default, nrow(meta)))
  }
  if (length(source) != 1L) {
    stop("metadata source must identify exactly one column")
  }
  if (!source %in% names(meta)) {
    stop("metadata source column is missing: ", source)
  }
  meta[[source]]
}

metadata_constant_value <- function(constants, name, n) {
  value <- constants[[name]]
  if (is.null(value)) {
    return(rep(NA_character_, n))
  }
  rep(as.character(value), n)
}

canonical_dataset_id <- function(dataset, level, value) {
  value <- normalize_missing_metadata(value)
  ifelse(
    is.na(value),
    NA_character_,
    paste(dataset, level, value, sep = ":")
  )
}

build_dataset_metadata <- function(
  raw_meta,
  dataset,
  field_map,
  constants = list(),
  default_analysis_role = "pending",
  default_analysis_include = TRUE,
  default_exclusion_reason = NA_character_
) {
  stopifnot(is.data.frame(raw_meta), length(dataset) == 1L)

  n <- nrow(raw_meta)
  cell_id <- metadata_source_value(raw_meta, field_map$cell_id)
  donor_raw <- metadata_source_value(raw_meta, field_map$donor_id)
  sample_raw <- metadata_source_value(raw_meta, field_map$sample_id)
  specimen_raw <- metadata_source_value(raw_meta, field_map$specimen_id)
  library_raw <- metadata_source_value(raw_meta, field_map$library_id)
  source_record_raw <- metadata_source_value(
    raw_meta,
    field_map$source_record_id
  )
  technical_batch_raw <- metadata_source_value(
    raw_meta,
    field_map$technical_batch_id
  )

  cell_id <- normalize_missing_metadata(cell_id)
  donor_raw <- normalize_missing_metadata(donor_raw)
  sample_raw <- normalize_missing_metadata(sample_raw)
  specimen_raw <- normalize_missing_metadata(specimen_raw)
  library_raw <- normalize_missing_metadata(library_raw)
  source_record_raw <- normalize_missing_metadata(source_record_raw)
  technical_batch_raw <- normalize_missing_metadata(technical_batch_raw)

  specimen_fallback <- is.na(specimen_raw)
  specimen_raw[specimen_fallback] <- sample_raw[specimen_fallback]
  specimen_fallback <- is.na(specimen_raw)
  specimen_raw[specimen_fallback] <- donor_raw[specimen_fallback]

  library_verification <- metadata_constant_value(
    constants,
    "library_id_verification_status",
    n
  )
  technical_batch_verification <- metadata_constant_value(
    constants,
    "technical_batch_verification_status",
    n
  )
  verified_library <- grepl(
    "^verified([ :]|$)",
    tolower(library_verification)
  )
  verified_technical_batch <- grepl(
    "^verified([ :]|$)",
    tolower(technical_batch_verification)
  )
  if (any(!is.na(library_raw) & !verified_library)) {
    stop(
      dataset,
      " library_id is populated without source-verified library semantics"
    )
  }
  if (any(!is.na(technical_batch_raw) & !verified_technical_batch)) {
    stop(
      dataset,
      paste(
        "technical_batch_id is populated without a paper- or",
        "repository-verified technical-batch definition"
      )
    )
  }

  get_field <- function(name) {
    source <- field_map[[name]]
    value <- metadata_source_value(raw_meta, source)
    if (all(is.na(value))) {
      value <- metadata_constant_value(constants, name, n)
    }
    normalize_missing_metadata(value)
  }
  get_field_source_column <- function(name, value) {
    source <- field_map[[name]]
    if (is.null(source) || length(source) == 0L || all(is.na(source))) {
      return(rep(NA_character_, n))
    }
    if (length(source) != 1L) {
      stop("metadata field source must identify exactly one column")
    }
    ifelse(is.na(value), NA_character_, as.character(source))
  }
  species <- get_field("species")
  species_raw <- metadata_source_value(
    raw_meta,
    field_map$species_raw
  )
  if (all(is.na(species_raw))) {
    species_raw <- species
  }
  species_raw <- normalize_missing_metadata(species_raw)
  sex_value <- get_field("sex")
  sex_assignment_method <- metadata_constant_value(
    constants,
    "sex_assignment_method",
    n
  )
  sex_assignment_method[
    !is.na(sex_value) & is.na(normalize_missing_metadata(
      sex_assignment_method
    ))
  ] <- "source reported"
  sex_assignment_method[is.na(sex_value)] <- "not reported"
  cell_type_label <- get_field("cell_type")
  cell_type_original_source <- field_map$cell_type_original
  if (is.null(cell_type_original_source) ||
    length(cell_type_original_source) == 0L ||
    all(is.na(cell_type_original_source))) {
    cell_type_original_source <- field_map$cell_type
  }
  cell_type_original_label <- normalize_missing_metadata(
    metadata_source_value(raw_meta, cell_type_original_source)
  )
  if (all(is.na(cell_type_original_label))) {
    cell_type_original_label <- cell_type_label
  }
  cell_type_level_1 <- get_field("cell_type_level_1")
  cell_type_level_2 <- get_field("cell_type_level_2")
  cell_type_level_3 <- get_field("cell_type_level_3")
  source_cluster_id <- get_field("cell_type_cluster_id")

  result <- data.frame(
    Cells = cell_id,
    Dataset = rep(dataset, n),
    Original_Cell_ID = cell_id,
    Original_Donor_ID = donor_raw,
    Original_Sample_ID = sample_raw,
    Original_Specimen_ID = specimen_raw,
    Original_Library_ID = library_raw,
    Original_Source_Sample_ID = sample_raw,
    Original_Source_Record_ID = source_record_raw,
    Original_Technical_Batch_ID = technical_batch_raw,
    Donor_ID = canonical_dataset_id(dataset, "donor", donor_raw),
    Specimen_ID = canonical_dataset_id(dataset, "specimen", specimen_raw),
    Library_ID = canonical_dataset_id(dataset, "library", library_raw),
    Source_Record_ID = canonical_dataset_id(
      dataset,
      "source_record",
      source_record_raw
    ),
    Technical_Batch_ID = canonical_dataset_id(
      dataset,
      "technical_batch",
      technical_batch_raw
    ),
    Global_Donor_ID = canonical_dataset_id(dataset, "donor", donor_raw),
    Donor_ID_Biological_Meaning = metadata_constant_value(
      constants,
      "donor_id_semantics",
      n
    ),
    Specimen_ID_Biological_Meaning = metadata_constant_value(
      constants,
      "specimen_id_semantics",
      n
    ),
    Library_ID_Biological_Meaning = metadata_constant_value(
      constants,
      "library_id_semantics",
      n
    ),
    Source_Record_ID_Biological_Meaning = metadata_constant_value(
      constants,
      "source_record_id_semantics",
      n
    ),
    Technical_Batch_ID_Biological_Meaning = metadata_constant_value(
      constants,
      "technical_batch_semantics",
      n
    ),
    Donor_ID_Verification_Status = ifelse(
      is.na(donor_raw),
      "unknown",
      "source reported"
    ),
    Specimen_ID_Verification_Status = ifelse(
      is.na(specimen_raw),
      "unknown",
      metadata_constant_value(
        constants,
        "specimen_id_verification_status",
        n
      )
    ),
    Library_ID_Verification_Status = ifelse(
      is.na(library_raw),
      "unknown",
      library_verification
    ),
    Library_ID_Evidence = metadata_constant_value(
      constants,
      "library_id_evidence",
      n
    ),
    Library_ID_Source_Column = get_field_source_column(
      "library_id",
      library_raw
    ),
    Library_ID_Derivation_Rule = metadata_constant_value(
      constants,
      "library_id_derivation_rule",
      n
    ),
    Technical_Batch_Verification_Status = ifelse(
      is.na(technical_batch_raw),
      "not available",
      technical_batch_verification
    ),
    Technical_Batch_Evidence = metadata_constant_value(
      constants,
      "technical_batch_evidence",
      n
    ),
    Technical_Batch_Source_Column = get_field_source_column(
      "technical_batch_id",
      technical_batch_raw
    ),
    Technical_Batch_Derivation_Rule = metadata_constant_value(
      constants,
      "technical_batch_derivation_rule",
      n
    ),
    Technical_Batch_Use_Status = ifelse(
      is.na(technical_batch_raw),
      "not available",
      metadata_constant_value(
        constants,
        "technical_batch_use_status",
        n
      )
    ),
    Technical_Batch_Use_Evidence = metadata_constant_value(
      constants,
      "technical_batch_use_evidence",
      n
    ),
    Species = species,
    Species_raw = species_raw,
    Age = get_field("age"),
    Sex = sex_value,
    Sex_Source_Raw = sex_value,
    Sex_Source_Standardized = ifelse(
      standardize_reported_sex(sex_value) %in% c("Female", "Male"),
      standardize_reported_sex(sex_value),
      NA_character_
    ),
    Sex_Source_Column = get_field_source_column("sex", sex_value),
    Sex_Assignment_Method = sex_assignment_method,
    BrainRegion = get_field("brain_region"),
    brain_region_source_label = get_field("brain_region_source"),
    brain_region_ontology_label = get_field(
      "brain_region_ontology_label"
    ),
    brain_region_ontology_id = get_field("brain_region_ontology_id"),
    brain_region_ontology_source = metadata_constant_value(
      constants,
      "brain_region_ontology_source",
      n
    ),
    Diagnosis_raw = get_field("diagnosis"),
    CellType_raw = cell_type_label,
    source_cell_type_original_label = cell_type_original_label,
    source_cell_type_label = cell_type_label,
    source_cell_type_level_1 = cell_type_level_1,
    source_cell_type_level_2 = cell_type_level_2,
    source_cell_type_level_3 = cell_type_level_3,
    source_cluster_id = source_cluster_id,
    source_cell_type_original_label_source_column = ifelse(
      is.na(cell_type_original_label),
      NA_character_,
      as.character(cell_type_original_source)
    ),
    source_cell_type_original_label_semantics = metadata_constant_value(
      constants,
      "cell_type_original_label_semantics",
      n
    ),
    source_cell_type_label_semantics = metadata_constant_value(
      constants,
      "cell_type_label_semantics",
      n
    ),
    source_cell_type_label_source_column = get_field_source_column(
      "cell_type",
      cell_type_label
    ),
    source_cell_type_level_1_source_column = get_field_source_column(
      "cell_type_level_1",
      cell_type_level_1
    ),
    source_cell_type_level_2_source_column = get_field_source_column(
      "cell_type_level_2",
      cell_type_level_2
    ),
    source_cell_type_level_3_source_column = get_field_source_column(
      "cell_type_level_3",
      cell_type_level_3
    ),
    source_cluster_id_source_column = get_field_source_column(
      "cell_type_cluster_id",
      source_cluster_id
    ),
    source_cell_type_ontology_label = get_field(
      "cell_type_ontology_label"
    ),
    source_cell_type_ontology_id = get_field(
      "cell_type_ontology_id"
    ),
    source_cell_type_source_column = metadata_constant_value(
      constants,
      "cell_type_source_column",
      n
    ),
    source_cell_type_source_file = metadata_constant_value(
      constants,
      "cell_type_source_file",
      n
    ),
    Technology = get_field("technology"),
    Sequence = get_field("modality"),
    Assay_Type = get_field("assay"),
    Sequencing_Platform = get_field("sequencing_platform"),
    Library_Chemistry = get_field("chemistry"),
    Analysis_Role = rep(default_analysis_role, n),
    Analysis_Include = rep(isTRUE(default_analysis_include), n),
    Exclusion_Reason = rep(
      normalize_missing_metadata(default_exclusion_reason),
      n
    ),
    Duplicate_Group = rep(NA_character_, n),
    Duplicate_Evidence = rep(NA_character_, n),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Backward-compatible aliases used by the existing integration scripts.
  result$Sample <- result$Donor_ID
  result$Sample_ID <- result$Specimen_ID

  original_columns <- setdiff(names(raw_meta), names(result))
  if (length(original_columns) > 0L) {
    result <- cbind(
      result,
      raw_meta[, original_columns, drop = FALSE]
    )
  }
  rownames(result) <- result$Cells
  result
}

join_external_metadata <- function(
  matrix_meta,
  external_meta,
  key_column = 1L,
  matrix_key_column = "Original_Cell_ID",
  match_column = "External_Metadata_Matched",
  require_all_matrix_cells = TRUE,
  require_all_metadata_cells = TRUE
) {
  if (length(matrix_key_column) != 1L ||
    !matrix_key_column %in% names(matrix_meta)) {
    stop(
      "matrix metadata key column is missing: ",
      matrix_key_column
    )
  }
  if (length(match_column) != 1L ||
    is.na(match_column) ||
    match_column == "" ||
    match_column %in% names(matrix_meta)) {
    stop("external metadata match column is invalid: ", match_column)
  }
  if (is.numeric(key_column)) {
    if (length(key_column) != 1L ||
      key_column < 1L ||
      key_column > ncol(external_meta)) {
      stop("external metadata key column index is invalid")
    }
    key_index <- as.integer(key_column)
    key_name <- names(external_meta)[[key_index]]
  } else {
    if (length(key_column) != 1L ||
      !key_column %in% names(external_meta)) {
      stop("external metadata key column is missing: ", key_column)
    }
    key_name <- key_column
    key_index <- match(key_name, names(external_meta))
  }

  matrix_ids <- normalize_missing_metadata(
    matrix_meta[[matrix_key_column]]
  )
  external_ids <- normalize_missing_metadata(external_meta[[key_index]])
  if (anyNA(external_ids) || any(external_ids == "")) {
    stop("external metadata contains missing cell IDs")
  }
  if (anyDuplicated(external_ids)) {
    stop("external metadata contains duplicated cell IDs")
  }

  idx <- match(matrix_ids, external_ids)
  matrix_without_metadata <- sum(is.na(idx))
  metadata_without_matrix <- sum(!external_ids %in% matrix_ids)
  if (isTRUE(require_all_matrix_cells) && matrix_without_metadata > 0L) {
    stop(
      "external metadata is missing ",
      matrix_without_metadata,
      " matrix cells"
    )
  }
  if (isTRUE(require_all_metadata_cells) && metadata_without_matrix > 0L) {
    stop(
      metadata_without_matrix,
      " external metadata cells are absent from the matrix"
    )
  }

  external_column_indices <- setdiff(
    seq_along(external_meta),
    key_index
  )
  external_columns <- names(external_meta)[external_column_indices]
  if (anyNA(external_columns) || any(external_columns == "")) {
    stop("external metadata contains unnamed non-key columns")
  }
  collisions <- intersect(names(matrix_meta), external_columns)
  if (length(collisions) > 0L) {
    stop(
      "external metadata column collision: ",
      paste(collisions, collapse = ", ")
    )
  }
  aligned <- external_meta[
    idx,
    external_column_indices,
    drop = FALSE
  ]
  rownames(aligned) <- NULL
  result <- cbind(matrix_meta, aligned)
  result[[match_column]] <- !is.na(idx)
  attr(result, "external_metadata_audit") <- data.frame(
    matrix_key_column = matrix_key_column,
    external_key_column = key_name,
    match_column = match_column,
    matrix_cells = length(matrix_ids),
    external_metadata_cells = length(external_ids),
    matched_cells = sum(!is.na(idx)),
    matrix_cells_without_metadata = matrix_without_metadata,
    metadata_cells_without_matrix = metadata_without_matrix,
    stringsAsFactors = FALSE
  )
  result
}

apply_donor_crosswalk <- function(meta, crosswalk_file) {
  if (!file.exists(crosswalk_file)) {
    return(meta)
  }
  crosswalk <- read.delim(
    crosswalk_file,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c(
    "Dataset", "Original_Donor_ID", "Global_Donor_ID",
    "Duplicate_Group", "Duplicate_Evidence",
    "Primary_Analysis_Include", "Analysis_Role",
    "Exclusion_Reason"
  )
  missing <- setdiff(required, names(crosswalk))
  if (length(missing) > 0L) {
    stop(
      "donor crosswalk is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  key <- paste(meta$Dataset, meta$Original_Donor_ID, sep = "\r")
  crosswalk_key <- paste(
    crosswalk$Dataset,
    crosswalk$Original_Donor_ID,
    sep = "\r"
  )
  idx <- match(key, crosswalk_key)
  matched <- !is.na(idx)
  meta$Global_Donor_ID[matched] <- normalize_missing_metadata(
    crosswalk$Global_Donor_ID[idx[matched]]
  )
  meta$Duplicate_Group[matched] <- normalize_missing_metadata(
    crosswalk$Duplicate_Group[idx[matched]]
  )
  meta$Duplicate_Evidence[matched] <- normalize_missing_metadata(
    crosswalk$Duplicate_Evidence[idx[matched]]
  )
  if (any(matched)) {
    include <- tolower(crosswalk$Primary_Analysis_Include[idx[matched]])
    if (any(!include %in% c("true", "false"))) {
      stop(
        "Primary_Analysis_Include must be true or false in donor crosswalk"
      )
    }
    meta$Analysis_Include[matched] <- include == "true"
    meta$Analysis_Role[matched] <- crosswalk$Analysis_Role[idx[matched]]
    meta$Exclusion_Reason[matched] <- normalize_missing_metadata(
      crosswalk$Exclusion_Reason[idx[matched]]
    )
  }
  meta
}

validate_dataset_metadata <- function(
  meta,
  matrix_cells = NULL,
  expected_cells = NULL
) {
  required <- c(
    "Cells", "Dataset", "Original_Cell_ID", "Original_Donor_ID",
    "Original_Specimen_ID", "Original_Library_ID",
    "Original_Source_Record_ID", "Original_Technical_Batch_ID",
    "Donor_ID", "Specimen_ID", "Library_ID", "Source_Record_ID",
    "Technical_Batch_ID", "Global_Donor_ID", "Analysis_Role",
    "Analysis_Include", "Exclusion_Reason",
    "Donor_ID_Biological_Meaning",
    "Specimen_ID_Biological_Meaning",
    "Library_ID_Biological_Meaning",
    "Source_Record_ID_Biological_Meaning",
    "Technical_Batch_ID_Biological_Meaning",
    "Donor_ID_Verification_Status",
    "Specimen_ID_Verification_Status",
    "Library_ID_Verification_Status",
    "Library_ID_Evidence", "Library_ID_Source_Column",
    "Library_ID_Derivation_Rule",
    "Technical_Batch_Verification_Status",
    "Technical_Batch_Evidence", "Technical_Batch_Source_Column",
    "Technical_Batch_Derivation_Rule", "Technical_Batch_Use_Status",
    "Technical_Batch_Use_Evidence", "Assay_Type",
    "Sequencing_Platform", "Library_Chemistry",
    "Sex_Source_Raw", "Sex_Source_Standardized", "Sex_Source_Column",
    "Sex_Assignment_Method",
    "sex_raw", "sex_standardized", "sex_provenance",
    "sex_standardization_status", "Sex_Donor_Conflict",
    "Sex_Specimen_Conflict",
    "CellType_raw", "source_cell_type_original_label",
    "source_cell_type_original_label_source_column",
    "source_cell_type_original_label_semantics",
    "source_cell_type_label", "source_cell_type_label_semantics",
    "source_cell_type_level_1", "source_cell_type_level_2",
    "source_cell_type_level_3", "source_cluster_id",
    "source_cell_type_label_source_column",
    "source_cell_type_level_1_source_column",
    "source_cell_type_level_2_source_column",
    "source_cell_type_level_3_source_column",
    "source_cluster_id_source_column",
    "source_cell_type_source_column", "source_cell_type_source_file",
    "source_cell_type_annotation_status",
    "source_cell_type_annotation_scope",
    "source_cell_type_mapping_method", "source_cell_type_available"
  )
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0L) {
    stop(
      "canonical metadata is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  if (anyNA(meta$Cells) || any(meta$Cells == "")) {
    stop("canonical metadata contains missing cell IDs")
  }
  if (anyDuplicated(meta$Cells)) {
    stop("canonical metadata contains duplicated cell IDs")
  }
  if (!is.null(expected_cells) && nrow(meta) != expected_cells) {
    stop(
      "full metadata cell count mismatch: expected ",
      expected_cells,
      ", observed ",
      nrow(meta)
    )
  }
  if (!is.null(matrix_cells) && !identical(meta$Cells, matrix_cells)) {
    stop("canonical metadata cell order does not match matrix columns")
  }
  if (any(!meta$Analysis_Include & is.na(meta$Exclusion_Reason))) {
    stop("excluded cells must have an explicit Exclusion_Reason")
  }
  if (any(meta$Analysis_Include & meta$Analysis_Role %in%
    c("recorded_not_integrated", "reference_excluded", "holdout_excluded"))) {
    stop("included cells carry an excluded analysis role")
  }
  library_present <- !is.na(normalize_missing_metadata(meta$Library_ID))
  library_verified <- grepl(
    "^verified([ :]|$)",
    tolower(as.character(meta$Library_ID_Verification_Status))
  )
  if (any(library_present & !library_verified)) {
    stop("canonical Library_ID contains an unverified source record")
  }
  library_evidence <- !is.na(normalize_missing_metadata(
    meta$Library_ID_Evidence
  ))
  if (any(library_present & !library_evidence)) {
    stop("canonical Library_ID lacks source evidence")
  }
  library_origin <-
    !is.na(normalize_missing_metadata(meta$Library_ID_Source_Column)) |
      !is.na(normalize_missing_metadata(meta$Library_ID_Derivation_Rule))
  if (any(library_present & !library_origin)) {
    stop("canonical Library_ID lacks a source column or derivation rule")
  }
  technical_batch_present <- !is.na(normalize_missing_metadata(
    meta$Technical_Batch_ID
  ))
  technical_batch_verified <- grepl(
    "^verified([ :]|$)",
    tolower(as.character(meta$Technical_Batch_Verification_Status))
  )
  if (any(technical_batch_present & !technical_batch_verified)) {
    stop("canonical Technical_Batch_ID contains an unverified batch label")
  }
  technical_batch_evidence <- !is.na(normalize_missing_metadata(
    meta$Technical_Batch_Evidence
  ))
  if (any(technical_batch_present & !technical_batch_evidence)) {
    stop("canonical Technical_Batch_ID lacks source evidence")
  }
  technical_batch_origin <-
    !is.na(normalize_missing_metadata(meta$Technical_Batch_Source_Column)) |
      !is.na(normalize_missing_metadata(meta$Technical_Batch_Derivation_Rule))
  if (any(technical_batch_present & !technical_batch_origin)) {
    stop("canonical Technical_Batch_ID lacks a source column or derivation rule")
  }
  technical_use_status <- tolower(as.character(
    meta$Technical_Batch_Use_Status
  ))
  eligible_technical <- grepl("^eligible([ :]|$)", technical_use_status)
  if (any(eligible_technical & !technical_batch_present)) {
    stop("technical-batch correction eligibility was asserted without a batch")
  }
  technical_use_evidence <- !is.na(normalize_missing_metadata(
    meta$Technical_Batch_Use_Evidence
  ))
  if (any(eligible_technical & !technical_use_evidence)) {
    stop("eligible technical batch lacks correction-use evidence")
  }
  annotation_available <- as.logical(meta$source_cell_type_available)
  exact_annotation_source <- Reduce(
    `|`,
    lapply(
      c(
        "source_cell_type_original_label_source_column",
        "source_cell_type_label_source_column",
        "source_cell_type_level_1_source_column",
        "source_cell_type_level_2_source_column",
        "source_cell_type_level_3_source_column",
        "source_cluster_id_source_column"
      ),
      function(column) {
        !is.na(normalize_missing_metadata(meta[[column]]))
      }
    )
  )
  if (any(annotation_available & !exact_annotation_source)) {
    stop("source cell-type annotation lacks its exact source column")
  }
  original_annotation_present <- !is.na(normalize_missing_metadata(
    meta$source_cell_type_original_label
  ))
  original_annotation_source_present <- !is.na(
    normalize_missing_metadata(
      meta$source_cell_type_original_label_source_column
    )
  )
  if (any(original_annotation_present & !original_annotation_source_present)) {
    stop("original-study cell-type annotation lacks its exact source column")
  }
  invisible(TRUE)
}

count_unique_metadata_rows <- function(meta, columns) {
  columns <- intersect(columns, names(meta))
  if (length(columns) == 0L) {
    return(data.frame(Cells = nrow(meta)))
  }
  values <- meta[, columns, drop = FALSE]
  if (requireNamespace("data.table", quietly = TRUE)) {
    grouped <- data.table::as.data.table(values)[
      ,
      .(Cells = .N),
      by = columns
    ]
    return(as.data.frame(grouped, stringsAsFactors = FALSE))
  }
  encoded <- lapply(values, function(column) {
    column <- as.character(column)
    column[is.na(column)] <- "<NA>"
    column
  })
  key <- do.call(paste, c(encoded, sep = "\r"))
  unique_key <- unique(key)
  first <- match(unique_key, key)
  counts <- tabulate(match(key, unique_key), nbins = length(unique_key))
  result <- values[first, , drop = FALSE]
  rownames(result) <- NULL
  result$Cells <- counts
  result
}

canonical_metadata_crosswalks <- function(meta) {
  crosswalks <- list(
    age = count_unique_metadata_rows(
      meta,
      c(
        "Dataset", "Age_Source_Raw", "Age_Source_Unit",
        "Age_Source_Basis", "Age_Source_Reference",
        "Age_Legacy_Processed",
        "Age_Source_Recovery_Method", "Age_Harmonization_Input",
        "Age_Conversion_Formula", "Age_Conversion_Confidence",
        "Age_Conversion_Applied", "Age", "Age_num",
        "Age_Lower", "Age_Upper", "Age_Canonical_Value",
        "Age_Canonical_Lower", "Age_Canonical_Upper",
        "Age_Canonical_Unit", "Unit", "Age_Representative_Method",
        "Age_Parser_Formula", "Age_Parser_Confidence", "Stage",
        "Age_Interval_Assignment_Status",
        "Age_Interval_Curation_Evidence",
        "Age_Interval_Curation_Confidence",
        "Continuous_Age_Eligibility",
        "Age_Interval_Boundary_Convention", "AgeIntervalID",
        "AgeInterval", "AgeRange"
      )
    ),
    region = count_unique_metadata_rows(
      meta,
      c(
        "Dataset", "brain_region_raw", "brain_region_source_label",
        "brain_region_source_column", "brain_region_standardized",
        "brain_region_harmonized", "brain_region_anatomical_level",
        "brain_region_laterality", "brain_region_ontology_label",
        "brain_region_mapping_method",
        "brain_region_mapping_confidence", "brain_region_ontology_id",
        "brain_region_ontology_source",
        "brain_region_ontology_mapping_status"
      )
    ),
    source_cell_type = count_unique_metadata_rows(
      meta,
      c(
        "Dataset", "Source_Accession", "Source_Publication_DOI",
        "Source_Publication_Title",
        "source_cell_type_original_label",
        "source_cell_type_original_label_source_column",
        "source_cell_type_original_label_semantics",
        "source_cell_type_label", "source_cell_type_label_semantics",
        "source_cell_type_level_1", "source_cell_type_level_2",
        "source_cell_type_level_3", "source_cluster_id",
        "source_cell_type_label_source_column",
        "source_cell_type_level_1_source_column",
        "source_cell_type_level_2_source_column",
        "source_cell_type_level_3_source_column",
        "source_cluster_id_source_column",
        "source_cell_type_source_column",
        "source_cell_type_source_file",
        "source_cell_type_annotation_status",
        "source_cell_type_annotation_scope",
        "source_cell_type_mapping_method",
        "source_cell_type_confidence",
        "source_cell_type_ontology_label",
        "source_cell_type_ontology_id",
        "source_cell_type_doublet_flag",
        "source_cell_type_available"
      )
    ),
    technical_batch = count_unique_metadata_rows(
      meta,
      c(
        "Dataset", "Original_Technical_Batch_ID",
        "Technical_Batch_ID",
        "Technical_Batch_ID_Biological_Meaning",
        "Technical_Batch_Verification_Status",
        "Technical_Batch_Evidence", "Technical_Batch_Source_Column",
        "Technical_Batch_Derivation_Rule", "Technical_Batch_Use_Status",
        "Technical_Batch_Use_Evidence", "Assay_Type",
        "sequencing_modality_raw",
        "sequencing_modality_standardized", "Technology",
        "Sequencing_Platform", "Library_Chemistry"
      )
    ),
    biological_unit = count_unique_metadata_rows(
      meta,
      c(
        "Dataset", "Original_Donor_ID",
        "Original_Source_Sample_ID", "Original_Source_Record_ID",
        "Global_Donor_ID", "Specimen_ID", "Library_ID",
        "Source_Record_ID", "Technical_Batch_ID",
        "Donor_ID_Biological_Meaning",
        "Specimen_ID_Biological_Meaning",
        "Library_ID_Biological_Meaning",
        "Library_ID_Evidence",
        "Library_ID_Source_Column", "Library_ID_Derivation_Rule",
        "Source_Record_ID_Biological_Meaning",
        "Technical_Batch_ID_Biological_Meaning",
        "Donor_ID_Verification_Status",
        "Specimen_ID_Verification_Status",
        "Library_ID_Verification_Status",
        "Technical_Batch_Verification_Status",
        "Technical_Batch_Evidence", "Technical_Batch_Source_Column",
        "Technical_Batch_Derivation_Rule", "Technical_Batch_Use_Status",
        "Technical_Batch_Use_Evidence", "Donor_ID_Derivation_Rule",
        "Sex_Source_Raw", "Sex_Source_Standardized", "Sex_Source_Column",
        "Sex_Assignment_Method", "Sex", "sex_standardized",
        "sex_provenance", "sex_standardization_status",
        "Sex_Donor_Conflict", "Sex_Specimen_Conflict",
        "Duplicate_Group", "Duplicate_Evidence", "Analysis_Include",
        "Analysis_Role", "Exclusion_Reason"
      )
    )
  )
  source("functions/canonical_donors.R")
  crosswalks$biological_unit <- add_canonical_donor_identity(
    crosswalks$biological_unit
  )
  for (name in names(crosswalks)) {
    if (sum(crosswalks[[name]]$Cells) != nrow(meta)) {
      stop(name, " metadata crosswalk does not cover every cell")
    }
  }
  crosswalks
}

metadata_single_value <- function(meta, column) {
  if (!column %in% names(meta)) {
    return(NA_character_)
  }
  value <- unique(stats::na.omit(normalize_missing_metadata(
    meta[[column]]
  )))
  if (length(value) == 0L) {
    return(NA_character_)
  }
  paste(sort(as.character(value)), collapse = ";")
}

metadata_distinct_count <- function(meta, column, selected = NULL) {
  if (!column %in% names(meta)) {
    return(0L)
  }
  if (is.null(selected)) {
    selected <- rep(TRUE, nrow(meta))
  }
  value <- normalize_missing_metadata(meta[[column]][selected])
  length(unique(stats::na.omit(value)))
}

metadata_value_examples <- function(
  meta,
  column,
  selected = NULL,
  limit = 30L
) {
  if (!column %in% names(meta)) {
    return(NA_character_)
  }
  if (is.null(selected)) {
    selected <- rep(TRUE, nrow(meta))
  }
  value <- sort(unique(stats::na.omit(normalize_missing_metadata(
    meta[[column]][selected]
  ))))
  if (length(value) == 0L) {
    return(NA_character_)
  }
  value <- gsub("[[:cntrl:]]+", " ", as.character(value))
  truncated <- length(value) > limit
  value <- utils::head(value, limit)
  paste0(
    paste(value, collapse = " | "),
    if (truncated) " | <additional values omitted>" else ""
  )
}

technical_batch_design_audit <- function(meta) {
  validate_dataset_metadata(meta)
  dataset <- metadata_single_value(meta, "Dataset")
  if (is.na(dataset) || grepl(";", dataset, fixed = TRUE)) {
    stop("technical-batch design audit requires exactly one dataset")
  }

  batch <- normalize_missing_metadata(meta$Technical_Batch_ID)
  known <- !is.na(batch)
  analysis_include <- if ("Analysis_Include" %in% names(meta)) {
    value <- as.logical(meta$Analysis_Include)
    !is.na(value) & value
  } else {
    rep(TRUE, nrow(meta))
  }
  verification <- normalize_missing_metadata(
    meta$Technical_Batch_Verification_Status
  )
  verified <- known & grepl(
    "paper-verified|repository-verified|verified",
    tolower(verification)
  )
  use_status <- normalize_missing_metadata(meta$Technical_Batch_Use_Status)
  eligible <- known & grepl(
    "^eligible([ :]|$)",
    tolower(use_status)
  )
  levels <- sort(unique(stats::na.omit(batch)))
  analysis_levels <- sort(unique(stats::na.omit(
    batch[analysis_include]
  )))

  make_row <- function(record_type, selected, batch_id = NA_character_) {
    # Missing batch identifiers do not belong to any known batch level.
    # Comparisons such as `batch == level` otherwise retain NA values and
    # propagate them into the level cell counts through sum().
    selected <- !is.na(selected) & selected
    donors <- metadata_distinct_count(meta, "Global_Donor_ID", selected)
    specimens <- metadata_distinct_count(meta, "Specimen_ID", selected)
    libraries <- metadata_distinct_count(meta, "Library_ID", selected)
    source_records <- metadata_distinct_count(
      meta,
      "Source_Record_ID",
      selected
    )
    independent_units <- max(donors, specimens)
    analysis_selected <- selected & analysis_include
    analysis_donors <- metadata_distinct_count(
      meta,
      "Global_Donor_ID",
      analysis_selected
    )
    analysis_specimens <- metadata_distinct_count(
      meta,
      "Specimen_ID",
      analysis_selected
    )
    analysis_libraries <- metadata_distinct_count(
      meta,
      "Library_ID",
      analysis_selected
    )
    analysis_independent_units <- max(
      analysis_donors,
      analysis_specimens
    )
    data.frame(
      Dataset = dataset,
      Source_Accession = metadata_single_value(meta, "Source_Accession"),
      Source_Publication_DOI = metadata_single_value(
        meta,
        "Source_Publication_DOI"
      ),
      Source_Publication_Title = metadata_single_value(
        meta,
        "Source_Publication_Title"
      ),
      Record_Type = record_type,
      Technical_Batch_ID = batch_id,
      Technical_Batch_ID_Biological_Meaning = metadata_single_value(
        meta[selected, , drop = FALSE],
        "Technical_Batch_ID_Biological_Meaning"
      ),
      Technical_Batch_Verification_Status = metadata_single_value(
        meta[selected, , drop = FALSE],
        "Technical_Batch_Verification_Status"
      ),
      Technical_Batch_Evidence = metadata_single_value(
        meta[selected, , drop = FALSE],
        "Technical_Batch_Evidence"
      ),
      Technical_Batch_Source_Column = metadata_single_value(
        meta[selected, , drop = FALSE],
        "Technical_Batch_Source_Column"
      ),
      Technical_Batch_Derivation_Rule = metadata_single_value(
        meta[selected, , drop = FALSE],
        "Technical_Batch_Derivation_Rule"
      ),
      Technical_Batch_Use_Status = metadata_single_value(
        meta[selected, , drop = FALSE],
        "Technical_Batch_Use_Status"
      ),
      Technical_Batch_Use_Evidence = metadata_single_value(
        meta[selected, , drop = FALSE],
        "Technical_Batch_Use_Evidence"
      ),
      Cells = sum(selected),
      Analysis_Cells = sum(analysis_selected),
      Donors = donors,
      Specimens = specimens,
      Libraries = libraries,
      Source_Records = source_records,
      Independent_Donor_Or_Specimen_Units = independent_units,
      Spans_Multiple_Donors = donors >= 2L,
      Spans_Multiple_Specimens = specimens >= 2L,
      Spans_Multiple_Donors_Or_Specimens = independent_units >= 2L,
      Analysis_Donors = analysis_donors,
      Analysis_Specimens = analysis_specimens,
      Analysis_Libraries = analysis_libraries,
      Analysis_Independent_Donor_Or_Specimen_Units =
        analysis_independent_units,
      Analysis_Spans_Multiple_Donors = analysis_donors >= 2L,
      Analysis_Spans_Multiple_Specimens = analysis_specimens >= 2L,
      Analysis_Spans_Multiple_Donors_Or_Specimens =
        analysis_independent_units >= 2L,
      Age_Intervals = metadata_distinct_count(meta, "Stage", selected),
      Age_Interval_Values = metadata_value_examples(
        meta,
        "Stage",
        selected
      ),
      Brain_Regions = metadata_distinct_count(
        meta,
        "brain_region_standardized",
        selected
      ),
      Brain_Region_Values = metadata_value_examples(
        meta,
        "brain_region_standardized",
        selected
      ),
      Reported_Sex_Values = metadata_distinct_count(
        meta,
        "sex_standardized",
        selected
      ),
      Sex_Values = metadata_value_examples(
        meta,
        "sex_standardized",
        selected
      ),
      Original_Cell_Type_Labels = metadata_distinct_count(
        meta,
        "source_cell_type_original_label",
        selected
      ),
      Original_Cell_Type_Examples = metadata_value_examples(
        meta,
        "source_cell_type_original_label",
        selected
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  level_rows <- lapply(levels, function(level) {
    make_row("batch_level", batch == level, level)
  })
  level_table <- if (length(level_rows) > 0L) {
    do.call(rbind, level_rows)
  } else {
    NULL
  }
  all_levels_cross_units <- length(levels) > 0L &&
    all(level_table$Spans_Multiple_Donors_Or_Specimens)
  full_coverage <- sum(known) == nrow(meta)
  analysis_full_coverage <- sum(known & analysis_include) ==
    sum(analysis_include)
  analysis_level_rows <- level_table$Technical_Batch_ID %in% analysis_levels
  analysis_all_levels_cross_units <- length(analysis_levels) > 0L &&
    all(level_table$Analysis_Spans_Multiple_Donors_Or_Specimens[
      analysis_level_rows
    ])
  all_known_verified <- any(known) && all(verified[known])
  source_declares_eligible <- any(eligible & analysis_include) &&
    all(eligible[known & analysis_include])
  design_gate_passed <- length(analysis_levels) >= 2L && full_coverage &&
    analysis_full_coverage && all_known_verified &&
    analysis_all_levels_cross_units

  summary <- make_row("dataset_summary", rep(TRUE, nrow(meta)))
  summary$Cells <- nrow(meta)
  summary$Technical_Batch_Known_Cells <- sum(known)
  summary$Technical_Batch_Missing_Cells <- sum(!known)
  summary$Technical_Batch_Coverage_Fraction <- if (nrow(meta) > 0L) {
    sum(known) / nrow(meta)
  } else {
    NA_real_
  }
  summary$Technical_Batch_Levels <- length(levels)
  summary$Analysis_Technical_Batch_Known_Cells <-
    sum(known & analysis_include)
  summary$Analysis_Technical_Batch_Missing_Cells <-
    sum(!known & analysis_include)
  summary$Analysis_Technical_Batch_Coverage_Fraction <- if (
    sum(analysis_include) > 0L
  ) {
    sum(known & analysis_include) / sum(analysis_include)
  } else {
    NA_real_
  }
  summary$Analysis_Technical_Batch_Levels <- length(analysis_levels)
  summary$Verified_Technical_Batch_Cells <- sum(verified)
  summary$Eligible_Technical_Batch_Cells <- sum(eligible)
  summary$Batch_Levels_Spanning_Multiple_Donors <- if (length(levels)) {
    sum(level_table$Spans_Multiple_Donors)
  } else {
    0L
  }
  summary$Batch_Levels_Spanning_Multiple_Specimens <- if (length(levels)) {
    sum(level_table$Spans_Multiple_Specimens)
  } else {
    0L
  }
  summary$All_Batch_Levels_Span_Multiple_Donors_Or_Specimens <-
    all_levels_cross_units
  summary$All_Analysis_Batch_Levels_Span_Multiple_Donors_Or_Specimens <-
    analysis_all_levels_cross_units
  summary$Technical_Batch_Design_Gate_Passed <- design_gate_passed
  summary$Declared_Eligibility_Consistent_With_Design_Gate <-
    !source_declares_eligible || design_gate_passed
  summary$Technical_Batch_Metadata_Completeness_Status <- if (!any(known)) {
    "not available"
  } else if (full_coverage) {
    "complete for all cells"
  } else {
    "incomplete; missing batch values remain"
  }
  summary$Biological_Confounding_Screen_Status <- if (!any(known)) {
    "not assessable without a verified technical-batch field"
  } else {
    paste(
      "requires review of batch-level donor, specimen, age, region, sex",
      "and original cell-type distributions"
    )
  }

  if (is.null(level_table)) {
    return(summary)
  }
  for (column in setdiff(names(summary), names(level_table))) {
    level_table[[column]] <- NA
  }
  level_table <- level_table[, names(summary), drop = FALSE]
  rbind(summary, level_table)
}

source_technical_field_inventory <- function(meta, source_meta = NULL) {
  validate_dataset_metadata(meta)
  if (is.null(source_meta)) {
    source_meta <- meta
  }
  if (!is.data.frame(source_meta) || nrow(source_meta) != nrow(meta)) {
    stop("source metadata and canonical metadata row counts differ")
  }
  source_id_column <- intersect(
    c("Source_Object_Cell_ID", "Original_Cell_ID", "Cells"),
    names(source_meta)
  )
  alignment_status <- "row count matched; source cell ID unavailable"
  if (length(source_id_column) > 0L) {
    source_ids <- as.character(source_meta[[source_id_column[[1L]]]])
    if (!identical(source_ids, as.character(meta$Cells))) {
      stop("source metadata and canonical metadata cell order differ")
    }
    alignment_status <- paste(
      "exact cell order verified using",
      source_id_column[[1L]]
    )
  }

  generated_pattern <- paste0(
    "^(Original_|Technical_Batch_|Integration_Batch_|Source_Accession$|",
    "Source_Publication_|Source_Repository|Library_ID$|Specimen_ID$|",
    "Global_Donor_ID$|Donor_ID$|Assay_Type$|Sequencing_Platform$|",
    "Library_Chemistry$|Technology$|Sequence$|sequencing_modality_)"
  )
  candidate_pattern <- paste0(
    "batch|library|lane|flow[._ -]*cell|gem[._ -]*well|multiplex|",
    "chemistry|platform|instrument|sequenc|assay|technology|",
    "(^|[._ -])run([._ -]|$)|(^|[._ -])pool([._ -]|$)"
  )
  generated <- grepl(generated_pattern, names(source_meta), ignore.case = TRUE)
  candidate <- grepl(candidate_pattern, names(source_meta), ignore.case = TRUE)
  candidate_columns <- names(source_meta)[candidate & !generated]
  exact_source <- unique(stats::na.omit(normalize_missing_metadata(
    meta$Technical_Batch_Source_Column
  )))
  missing_exact_source <- setdiff(exact_source, names(source_meta))
  if (length(missing_exact_source) > 0L) {
    stop(
      "audited technical-batch source column is absent from source metadata: ",
      paste(missing_exact_source, collapse = ", ")
    )
  }
  derivation_rule <- unique(stats::na.omit(normalize_missing_metadata(
    meta$Technical_Batch_Derivation_Rule
  )))
  direct_identity_status <- NA_character_
  if (length(exact_source) == 1L && length(derivation_rule) == 0L) {
    source_value <- normalize_missing_metadata(source_meta[[exact_source]])
    audited_value <- normalize_missing_metadata(
      meta$Original_Technical_Batch_ID
    )
    if (!identical(source_value, audited_value)) {
      stop(
        "direct technical-batch source values differ from the audited ",
        "original technical-batch values"
      )
    }
    direct_identity_status <-
      "exact direct-source value identity verified for every row"
  } else if (length(exact_source) > 0L && length(derivation_rule) > 0L) {
    direct_identity_status <-
      "derived field; use the recorded derivation rule instead of identity"
  }
  candidate_columns <- unique(c(candidate_columns, exact_source))
  candidate_columns <- candidate_columns[candidate_columns %in% names(source_meta)]

  base_row <- function(record_type, column = NA_character_) {
    value <- if (!is.na(column)) {
      normalize_missing_metadata(source_meta[[column]])
    } else {
      rep(NA_character_, nrow(source_meta))
    }
    exact <- !is.na(column) && column %in% exact_source
    data.frame(
      Dataset = metadata_single_value(meta, "Dataset"),
      Source_Accession = metadata_single_value(meta, "Source_Accession"),
      Source_Publication_DOI = metadata_single_value(
        meta,
        "Source_Publication_DOI"
      ),
      Source_Publication_Title = metadata_single_value(
        meta,
        "Source_Publication_Title"
      ),
      Record_Type = record_type,
      Source_Field = column,
      Source_Metadata_Columns = ncol(source_meta),
      Source_Metadata_Column_Names = paste(
        names(source_meta),
        collapse = ";"
      ),
      Source_Metadata_Rows = nrow(source_meta),
      Nonmissing_Rows = if (is.na(column)) 0L else sum(!is.na(value)),
      Missing_Rows = if (is.na(column)) nrow(source_meta) else sum(is.na(value)),
      Coverage_Fraction = if (is.na(column) || nrow(source_meta) == 0L) {
        NA_real_
      } else {
        sum(!is.na(value)) / nrow(source_meta)
      },
      Distinct_Nonmissing_Values = if (is.na(column)) {
        0L
      } else {
        length(unique(stats::na.omit(value)))
      },
      Example_Values = if (is.na(column)) {
        NA_character_
      } else {
        metadata_value_examples(source_meta, column, limit = 12L)
      },
      Exact_Audited_Technical_Batch_Source = exact,
      Technical_Batch_Verification_Status = if (exact) {
        metadata_single_value(meta, "Technical_Batch_Verification_Status")
      } else {
        "not semantically verified"
      },
      Technical_Batch_Evidence = if (exact) {
        metadata_single_value(meta, "Technical_Batch_Evidence")
      } else {
        NA_character_
      },
      Technical_Batch_Derivation_Rule = if (exact) {
        metadata_single_value(meta, "Technical_Batch_Derivation_Rule")
      } else {
        NA_character_
      },
      Source_To_Audited_Batch_Value_Status = if (exact) {
        direct_identity_status
      } else {
        NA_character_
      },
      Semantic_Admission_Status = if (exact) {
        paste(
          "linked to the explicitly audited technical-batch source field;",
          "eligibility is determined separately by the design audit"
        )
      } else if (identical(record_type, "candidate_field")) {
        paste(
          "candidate name only; do not use without paper or repository",
          "evidence and a passing design audit"
        )
      } else {
        "no candidate technical source field detected"
      },
      Source_Cell_Order_Alignment = alignment_status,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  summary <- base_row("dataset_summary")
  summary$Candidate_Source_Fields <- length(candidate_columns)
  if (length(candidate_columns) == 0L) {
    return(summary)
  }
  rows <- lapply(candidate_columns, function(column) {
    base_row("candidate_field", column)
  })
  inventory <- do.call(rbind, rows)
  inventory$Candidate_Source_Fields <- NA_integer_
  rbind(summary, inventory)
}

canonical_metadata_audits <- function(meta, source_meta = NULL) {
  list(
    technical_batch_design = technical_batch_design_audit(meta),
    source_technical_field_inventory = source_technical_field_inventory(
      meta,
      source_meta = source_meta
    )
  )
}

write_metadata_table_atomic <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(file, ".tmp.", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  gzip_output <- grepl("[.]gz$", file, ignore.case = TRUE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(
      x,
      temporary,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "NA",
      compress = if (gzip_output) "gzip" else "none"
    )
  } else if (gzip_output) {
    connection <- gzfile(temporary, "wt")
    tryCatch(
      write.table(
        x,
        connection,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
      ),
      finally = close(connection)
    )
  } else {
    write.table(
      x,
      temporary,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "NA"
    )
  }
  if (file.exists(file)) {
    unlink(file)
  }
  if (!file.rename(temporary, file)) {
    stop("failed to atomically write metadata table: ", file)
  }
  invisible(file)
}

write_canonical_metadata_crosswalks <- function(
  meta,
  output_dir,
  source_meta = NULL
) {
  crosswalks <- canonical_metadata_crosswalks(meta)
  files <- c(
    age = "age_crosswalk.tsv",
    region = "region_crosswalk.tsv",
    source_cell_type = "source_cell_type_crosswalk.tsv",
    technical_batch = "technical_batch_crosswalk.tsv",
    biological_unit = "donor_specimen_library_crosswalk.tsv.gz"
  )
  for (name in names(files)) {
    write_metadata_table_atomic(
      crosswalks[[name]],
      file.path(output_dir, files[[name]])
    )
  }
  audits <- canonical_metadata_audits(meta, source_meta = source_meta)
  audit_files <- c(
    technical_batch_design = "technical_batch_design_audit.tsv",
    source_technical_field_inventory =
      "source_technical_field_inventory.tsv"
  )
  for (name in names(audit_files)) {
    write_metadata_table_atomic(
      audits[[name]],
      file.path(output_dir, audit_files[[name]])
    )
  }
  unname(c(files, audit_files))
}

dataset_metadata_audit <- function(meta) {
  count_true <- function(column, default = FALSE) {
    if (!column %in% names(meta)) {
      return(sum(rep(default, nrow(meta))))
    }
    value <- as.logical(meta[[column]])
    value[is.na(value)] <- FALSE
    sum(value)
  }
  count_reported <- function(column) {
    if (!column %in% names(meta)) {
      return(0L)
    }
    sum(!is.na(normalize_missing_metadata(meta[[column]])))
  }
  data.frame(
    Dataset = unique(meta$Dataset),
    full_cells = nrow(meta),
    included_cells = sum(meta$Analysis_Include),
    excluded_cells = sum(!meta$Analysis_Include),
    reported_donors = length(unique(stats::na.omit(
      meta$Global_Donor_ID
    ))),
    reported_specimens = length(unique(stats::na.omit(
      meta$Specimen_ID
    ))),
    reported_libraries = length(unique(stats::na.omit(
      meta$Library_ID
    ))),
    verified_technical_batches = length(unique(stats::na.omit(
      meta$Technical_Batch_ID
    ))),
    correction_eligible_technical_batch_cells = sum(grepl(
      "^eligible([ :]|$)",
      tolower(as.character(meta$Technical_Batch_Use_Status))
    )),
    source_age_reported_cells = count_reported("Age_Source_Raw"),
    canonical_age_parsed_cells = if ("Age_num" %in% names(meta)) {
      sum(!is.na(meta$Age_num))
    } else {
      0L
    },
    age_interval_assigned_cells = count_reported("Stage"),
    brain_region_reported_cells = count_reported(
      "brain_region_source_label"
    ),
    brain_region_ontology_mapped_cells = count_reported(
      "brain_region_ontology_id"
    ),
    source_cell_type_available_cells = count_true(
      "source_cell_type_available"
    ),
    source_cell_type_original_label_cells = count_reported(
      "source_cell_type_original_label"
    ),
    source_cell_type_unavailable_cells = nrow(meta) - count_true(
      "source_cell_type_available"
    ),
    source_cell_type_level_1_cells = count_reported(
      "source_cell_type_level_1"
    ),
    source_cell_type_level_2_cells = count_reported(
      "source_cell_type_level_2"
    ),
    source_cell_type_level_3_cells = count_reported(
      "source_cell_type_level_3"
    ),
    source_cluster_id_cells = count_reported("source_cluster_id"),
    source_cell_type_exact_source_cells = sum(Reduce(
      `|`,
      lapply(
        c(
          "source_cell_type_original_label_source_column",
          "source_cell_type_label_source_column",
          "source_cell_type_level_1_source_column",
          "source_cell_type_level_2_source_column",
          "source_cell_type_level_3_source_column",
          "source_cluster_id_source_column"
        ),
        function(column) {
          !is.na(normalize_missing_metadata(meta[[column]]))
        }
      )
    )),
    reported_sex_cells = count_reported("Sex_Source_Raw"),
    standardized_sex_cells = count_reported("Sex"),
    donor_sex_conflict_cells = if ("Sex_Donor_Conflict" %in% names(meta)) {
      sum(as.logical(meta$Sex_Donor_Conflict), na.rm = TRUE)
    } else {
      0L
    },
    specimen_sex_conflict_cells = if (
      "Sex_Specimen_Conflict" %in% names(meta)
    ) {
      sum(as.logical(meta$Sex_Specimen_Conflict), na.rm = TRUE)
    } else {
      0L
    },
    full_cell_fraction = 1,
    metadata_schema_version = brainomics_metadata_schema_version(),
    stringsAsFactors = FALSE
  )
}
