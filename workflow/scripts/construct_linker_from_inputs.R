library(openxlsx, quietly = TRUE)
library(stringr, quietly = TRUE)

#' Combine subject ID, analyte ID, and sequencing index ID
#' into a '_'-delimited identifier for exported data preparation.
#'
#' @details
#' If any of the required columns have an NA entry for a particular subject,
#' the entire resultant ID is set to NA.
#'
#' @param df data.frame; input data minimally containing columns 'subject',
#' 'analyte', and 'index'
#' @return character vector; constructed output IDs for each row in input
#' data.frame, or NA if any of the required inputs for the row were NA
construct.output.stems <- function(df) {
  stopifnot(is.data.frame(df))
  stopifnot(
    "input data frame contains column 'subject'" = "subject" %in% colnames(df),
    "input data frame contains column 'analyte'" = "analyte" %in% colnames(df),
    "input data frame contains column 'index'" = "index" %in% colnames(df)
  )
  ## try to construct IDs of the format SUBJECTID_LSID_SQID
  output.stems <- paste(df$subject, df$analyte, df$index, sep = "_")
  output.stems[is.na(df$subject) | is.na(df$analyte) | is.na(df$index)] <- NA
  output.stems
}

#' Search for freetext ID remappings and apply them as needed
#'
#' The idea here is that there are evidently some freetext
#' notes along the lines of
#' "SAMPLE SWAP - this is actually PMGRC-\d+-\d+-\d" that
#' only exist in freetext in rows that are labeled by
#' _unmapped_ subject ID. For downstream analysis, these ID
#' mappings must be detected and applied. To make things
#' even more confusing, there doesn't appear to be any sense
#' to which columns of the logbook this information is reported in.
#'
#' The logic in this initial implementation will likely need to be
#' made more complex when different upstream freetext field types
#' are detected.
#'
#' @param df data.frame, input logbook sheet from read.xlsx
#' @param vec character vector, vector of actual subject IDs from
#' current sheet
#' @return character vector, mapped version of input character vector
#' with ID relabelings applied
apply.id.mappings <- function(df, vec) {
  stopifnot(is.data.frame(df))
  stopifnot(is.character(vec))
  stopifnot("input vector should represent subject ID column of data frame" = nrow(df) == length(vec))
  known.patterns <- c("^SAMPLE SWAP.*this is actually (PMGRC-\\d+-\\d+-\\d).*$" = "\\1")
  for (i in seq_len(ncol(df))) {
    new.ids <- rep(NA, length(vec))
    for (j in seq_len(length(known.patterns))) {
      has.pattern <- stringr::str_detect(df[, i], names(known.patterns)[j])
      has.pattern[is.na(has.pattern)] <- FALSE
      new.ids[has.pattern] <- stringr::str_replace(
        df[has.pattern, i],
        names(known.patterns)[j],
        known.patterns[j]
      )
    }
    vec[!is.na(new.ids)] <- new.ids[!is.na(new.ids)]
  }
  vec
}

#' Parse legacy manual logbook data, attempt to minimally standardize,
#' and convert into some sort of useful format.
#'
#' @description
#' For legacy support, this function supports
#' input in the form of a manually annotated Excel document created by lab
#' staff. The format largely defies description, and varies from tab to tab.
#' For any future applications, this functionality should _never_ be used;
#' rather, input information should be provided via the data type-specific
#' linker files.
#'
#' @param input.fn character; name of spreadsheet containing logbook data.
#' In practical use, this was a local copy pulled from a google sheet,
#' to make individual runs of the workflow comparatively resilient to the
#' mercurial changes of the file upstream.
#' @return data.frame; formatted linker data
parse.logbook <- function(input.fn) {
  stopifnot(
    is.character(input.fn),
    length(input.fn) == 1
  )
  sheet.names <- openxlsx::getSheetNames(input.fn)
  subject.id <- c()
  jira.tickets <- c()
  index.id <- c()
  analyte.id <- c()
  project.id <- c()
  sex.id <- c()
  for (sheet.name in sheet.names) {
    df <- openxlsx::read.xlsx(input.fn, sheet = sheet.name, check.names = FALSE)
    colnames(df) <- tolower(colnames(df))
    if (colnames(df)[1] == "pmgrc.id") {
      ## resolve chaos
      subject.col <- 1
      jira.col <- which(stringr::str_detect(colnames(df), "jira\\.ticket\\.for\\.batches\\.in\\.flight"))
      index.col <- which(stringr::str_detect(colnames(df), "sq\\.id"))
      analyte.col <- which(stringr::str_detect(colnames(df), "ls\\.id"))
      project.col <- which(stringr::str_detect(colnames(df), "ru\\.id"))
      sex.col <- which(stringr::str_detect(colnames(df), "biological\\.sex"))
      if (length(jira.col) == 0) {
        jira.col <- which(stringr::str_detect(colnames(df), "^jira\\.ticket\\(s\\)\\.\\("))
      }
      if (length(jira.col) != 1) {
        next
      }
      ## "fix": upstream logbook has sporadic freetext annotations that indicate when a
      ## subject needs to have an ID update applied to it. for some IDs, this information
      ## is not recorded but rather just applied without annotation. It's not clear what
      ## governs the difference between those two situations
      df[, subject.col] <- apply.id.mappings(df, df[, subject.col])
      for (row.num in seq_len(nrow(df))) {
        subject.id <- c(subject.id, df[row.num, subject.col])
        jira.tickets <- c(jira.tickets, df[row.num, jira.col])
        index.id <- c(index.id, df[row.num, index.col])
        analyte.id <- c(analyte.id, df[row.num, analyte.col])
        project.id <- c(project.id, df[row.num, project.col])
        if (length(sex.col) == 1) {
          sex.id <- c(sex.id, df[row.num, sex.col])
        } else {
          sex.id <- c(sex.id, "Unknown")
        }
      }
    }
  }
  res <- data.frame(
    subject = subject.id,
    jira = jira.tickets,
    project = project.id,
    index = index.id,
    analyte = analyte.id,
    sex = sex.id
  )
  ## remove entries with useless content
  res <- res[!(is.na(res[, 1]) & is.na(res[, 2]) & is.na(res[, 3]) & is.na(res[, 4]) & is.na(res[, 5])), ]
  ## construct output stems for deliverables, if sufficient information is present to fit the accepted format
  output.stems <- construct.output.stems(res)
  res$external <- output.stems
  res
}

#' Contextually overwrite or append linker information to an existing
#' data frame from a simple two column id->value linker file
#'
#' @details
#' If entirely new subjects are encountered in the linker file,
#' those subjects will be appended to the existing data frame
#' with relevant linker information in the appropriate column
#' and NA everywhere else
#'
#' @param df data.frame existing linker information, to which new
#' information will be overwritten or appended
#' @param linker.fn character; name of simple linker file with new
#' information to be added to the existing data frame
#' @param target.colname character; name of column in existing data
#' frame to which linker information will be saved. column does not
#' necessarily need to already exist, although it is intended to
#' @return data.frame; modified version of input data frame with
#' new information overwritten or appended
add.linker.data <- function(df, linker.fn, target.colname) {
  stopifnot(is.data.frame(df))
  stopifnot(
    is.character(linker.fn),
    length(linker.fn) == 1
  )
  stopifnot(is.character(target.colname))
  linker.df <- read.table(linker.fn, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
  rownames(linker.df) <- linker.df[, 1]
  ## handle the situation where people not yet present are being added here
  new.subjects <- linker.df[!(linker.df[, 1] %in% df$subject), 1]
  if (length(new.subjects) > 0) {
    new.df <- data.frame(
      subject = new.subjects,
      jira = NA,
      project = NA,
      index = new.subjects,
      analyte = NA,
      sex = NA,
      external = NA
    )
    df <- rbind(df, new.df)
  }
  df[df[, "subject"] %in% linker.df[, 1], target.colname] <-
    linker.df[df[df[, "subject"] %in% linker.df[, 1], "subject"], 2]
  df
}

#' Run primary logic of this script, wrapped such that sourcing
#' this file will not cause actual code execution.
#'
#' @param logbook.fn character or NULL; name of input legacy
#' logbook file. Can be NULL, in which case it is effectively ignored
#' @param sex.linker.fn character or NULL; name of input sex linker
#' file. Can be NULL, in which case it is effectively ignored
#' @param external.id.linker.fn character or NULL; name of input external
#' ID linker file. Can be NULL, in which case it is effectively ignored
#' @param out.fn character; name of file to which to write output linker data
run.construct.linker <- function(logbook.fn,
                                 sex.linker.fn,
                                 external.id.linker.fn,
                                 out.fn) {
  stopifnot(
    is.character(logbook.fn) || is.null(logbook.fn),
    length(logbook.fn) <= 1
  )
  stopifnot(
    is.character(sex.linker.fn) || is.null(sex.linker.fn),
    length(sex.linker.fn) <= 1
  )
  stopifnot(
    is.character(external.id.linker.fn) || is.null(external.id.linker.fn),
    length(external.id.linker.fn) <= 1
  )
  stopifnot(
    is.character(out.fn),
    length(out.fn) == 1
  )
  df <- data.frame(
    subject = c("A"),
    jira = c("A"),
    project = c("A"),
    index = c("A"),
    analyte = c("A"),
    sex = c("A")
  )
  df <- df[-1, ]
  if (!is.null(logbook.fn)) {
    df <- parse.logbook(logbook.fn)
  }
  if (!is.null(sex.linker.fn)) {
    df <- add.linker.data(df, sex.linker.fn, "sex")
  }
  if (!is.null(external.id.linker.fn)) {
    df <- add.linker.data(df, external.id.linker.fn, "external")
  }
  ## deal with possibility that no external ID has been provided
  df[is.na(df[, "external"]), "external"] <- df[is.na(df[, "external"]), "subject"]

  ## deal with apparent internal newlines
  for (index in seq_len(ncol(df))) {
    df[, index] <- stringr::str_replace_all(df[, index], "\n", " ")
  }

  write.table(df, out.fn, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

if (exists("snakemake")) {
  run.construct.linker(
    snakemake@params[["logbook"]],
    snakemake@params[["sex_linker"]],
    snakemake@params[["external_id_linker"]],
    snakemake@output[["linker"]]
  )
}
