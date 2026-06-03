#!/usr/bin/env Rscript

# EC metadata -> activity mapping for LPMO analyses.
#
# Usage:
#   Rscript scripts/ec_activity_mapping.R \
#     --input metadata.csv \
#     --output activity_mapping.tsv \
#     --id-col protein_id \
#     --ec-col ec_number \
#     --family-col family

parse_args <- function(args) {
  opts <- list(
    input = NULL,
    output = NULL,
    id_col = "protein_id",
    ec_col = "ec_number",
    family_col = "family"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--input" && i < length(args)) {
      opts$input <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--output" && i < length(args)) {
      opts$output <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--id-col" && i < length(args)) {
      opts$id_col <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--ec-col" && i < length(args)) {
      opts$ec_col <- args[[i + 1]]
      i <- i + 2
    } else if (key == "--family-col" && i < length(args)) {
      opts$family_col <- args[[i + 1]]
      i <- i + 2
    } else {
      stop(sprintf("Unknown or incomplete argument: %s", key), call. = FALSE)
    }
  }

  if (is.null(opts$input) || is.null(opts$output)) {
    stop(
      "Missing required arguments. Required: --input <file> --output <file>",
      call. = FALSE
    )
  }

  opts
}

pick_separator <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.tsv$", lower) || grepl("\\.txt$", lower)) {
    return("\t")
  }
  ","
}

read_table_auto <- function(path) {
  sep <- pick_separator(path)
  read.table(
    path,
    sep = sep,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "\"",
    comment.char = ""
  )
}

write_table_auto <- function(df, path) {
  sep <- pick_separator(path)
  write.table(df, file = path, sep = sep, row.names = FALSE, quote = TRUE)
}

normalize_ec <- function(ec_raw) {
  ec <- trimws(as.character(ec_raw))
  if (is.na(ec) || ec == "") {
    return(character(0))
  }
  matches <- gregexpr("1\\.14\\.99\\.[0-9-]+", ec)
  tokens <- regmatches(ec, matches)[[1]]
  if (length(tokens) == 1 && identical(tokens, "")) {
    return(character(0))
  }
  unique(tokens[nzchar(tokens)])
}

map_ec_token <- function(ec_token, family_raw) {
  family <- toupper(trimws(as.character(family_raw)))

  if (identical(ec_token, "1.14.99.54")) {
    return(list(
      substrate_class = "cellulose",
      regio_class = "C1",
      activity_label = "cellulose_C1_hydroxylating",
      mapping_rule = "exact_1.14.99.54"
    ))
  }

  if (identical(ec_token, "1.14.99.56")) {
    return(list(
      substrate_class = "cellulose",
      regio_class = "C4",
      activity_label = "cellulose_C4_dehydrogenating",
      mapping_rule = "exact_1.14.99.56"
    ))
  }

  if (identical(ec_token, "1.14.99.53")) {
    return(list(
      substrate_class = "chitin",
      regio_class = "C1",
      activity_label = "chitin_C1_hydroxylating",
      mapping_rule = "exact_1.14.99.53"
    ))
  }

  if (identical(ec_token, "1.14.99.55")) {
    return(list(
      substrate_class = "starch",
      regio_class = "C1",
      activity_label = "starch_C1_hydroxylating",
      mapping_rule = "exact_1.14.99.55"
    ))
  }

  if (identical(ec_token, "1.14.99.-")) {
    if (grepl("AA10", family, fixed = TRUE)) {
      return(list(
        substrate_class = "unknown",
        regio_class = "C4",
        activity_label = "unknown_substrate_C4",
        mapping_rule = "special_1.14.99.-_AA10"
      ))
    }
    if (grepl("AA14", family, fixed = TRUE)) {
      return(list(
        substrate_class = "unknown",
        regio_class = "C1",
        activity_label = "unknown_substrate_C1",
        mapping_rule = "special_1.14.99.-_AA14"
      ))
    }
    if (grepl("AA17", family, fixed = TRUE)) {
      return(list(
        substrate_class = "unknown",
        regio_class = "C4",
        activity_label = "unknown_substrate_C4",
        mapping_rule = "special_1.14.99.-_AA17"
      ))
    }
    return(list(
      substrate_class = "unknown",
      regio_class = "unknown",
      activity_label = "unknown",
      mapping_rule = "special_1.14.99.-_non_AA10_AA14_AA17"
    ))
  }

  NULL
}

map_ec_activity <- function(ec_raw, family_raw) {
  tokens <- normalize_ec(ec_raw)
  if (length(tokens) == 0) {
    return(list(
      substrate_class = "unknown",
      regio_class = "unknown",
      activity_label = "unknown",
      mapping_rule = "missing_ec"
    ))
  }

  mapped_items <- list()
  for (token in tokens) {
    mapped_item <- map_ec_token(token, family_raw)
    if (!is.null(mapped_item)) {
      mapped_items[[length(mapped_items) + 1]] <- mapped_item
    }
  }

  if (length(mapped_items) > 0) {
    substrate_values <- unique(vapply(mapped_items, function(x) x$substrate_class, character(1)))
    substrate_values <- substrate_values[substrate_values %in% c("chitin", "cellulose", "starch")]
    regio_values <- unique(tolower(vapply(mapped_items, function(x) x$regio_class, character(1))))
    has_c1 <- "c1" %in% regio_values
    has_c4 <- "c4" %in% regio_values
    regio_class <- "unknown"
    if (has_c1 && has_c4) {
      regio_class <- "C1/C4"
    } else if (has_c1) {
      regio_class <- "C1"
    } else if (has_c4) {
      regio_class <- "C4"
    }
    return(list(
      substrate_class = if (length(substrate_values) > 0) paste(sort(substrate_values), collapse = "+") else "unknown",
      regio_class = regio_class,
      activity_label = paste(sort(unique(vapply(mapped_items, function(x) x$activity_label, character(1)))), collapse = "+"),
      mapping_rule = paste(sort(unique(vapply(mapped_items, function(x) x$mapping_rule, character(1)))), collapse = "+")
    ))
  }

  list(
    substrate_class = "unknown",
    regio_class = "unknown",
    activity_label = "unknown",
    mapping_rule = "unmapped_ec"
  )
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- parse_args(args)

  metadata <- read_table_auto(opts$input)

  required_cols <- c(opts$id_col, opts$ec_col)
  missing_cols <- setdiff(required_cols, names(metadata))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!(opts$family_col %in% names(metadata))) {
    metadata[[opts$family_col]] <- ""
  }

  mapped_list <- mapply(
    FUN = map_ec_activity,
    ec_raw = metadata[[opts$ec_col]],
    family_raw = metadata[[opts$family_col]],
    SIMPLIFY = FALSE
  )

  metadata$mapped_substrate_class <- vapply(mapped_list, function(x) x$substrate_class, character(1))
  metadata$mapped_regio_class <- vapply(mapped_list, function(x) x$regio_class, character(1))
  metadata$mapped_activity_label <- vapply(mapped_list, function(x) x$activity_label, character(1))
  metadata$mapping_rule <- vapply(mapped_list, function(x) x$mapping_rule, character(1))

  write_table_auto(metadata, opts$output)

  message(sprintf("Wrote %d mapped rows to %s", nrow(metadata), opts$output))
}

main()
