library(openxlsx)
library(stringr)

construct.output.stems <- function(df) {
  ## try to construct IDs of the format PMGRCID_LSID_SQID
  output.stems <- paste(df$pmgrc, df$ls, df$sq, sep = "_")
  output.stems[is.na(df$pmgrc) | is.na(df$ls) | is.na(df$sq)] <- NA
  output.stems
}

run.construct.linker <- function(input.fn, output.fn) {
  sheet.names <- openxlsx::getSheetNames(input.fn)
  pmgrc.id <- c()
  jira.tickets <- c()
  sq.id <- c()
  ls.id <- c()
  ru.id <- c()
  sex.id <- c()
  for (sheet.name in sheet.names) {
    df <- openxlsx::read.xlsx(input.fn, sheet = sheet.name, check.names = FALSE)
    colnames(df) <- tolower(colnames(df))
    if (colnames(df)[1] == "pmgrc.id") {
      ## resolve chaos
      pmgrc.col <- 1
      jira.col <- which(stringr::str_detect(colnames(df), "jira\\.ticket\\.for\\.batches\\.in\\.flight"))
      sq.col <- which(stringr::str_detect(colnames(df), "sq\\.id"))
      ls.col <- which(stringr::str_detect(colnames(df), "ls\\.id"))
      ru.col <- which(stringr::str_detect(colnames(df), "ru\\.id"))
      sex.col <- which(stringr::str_detect(colnames(df), "biological\\.sex"))
      if (length(jira.col) != 1) {
        next
      }
      for (row.num in seq_len(nrow(df))) {
        pmgrc.id <- c(pmgrc.id, df[row.num, pmgrc.col])
        jira.tickets <- c(jira.tickets, df[row.num, jira.col])
        sq.id <- c(sq.id, df[row.num, sq.col])
        ls.id <- c(ls.id, df[row.num, ls.col])
        ru.id <- c(ru.id, df[row.num, ru.col])
        if (length(sex.col) == 1) {
          sex.id <- c(sex.id, df[row.num, sex.col])
        } else {
          sex.id <- c(sex.id, "Unknown")
        }
      }
    }
  }
  res <- data.frame(
    pmgrc = pmgrc.id,
    jira = jira.tickets,
    ru = ru.id,
    sq = sq.id,
    ls = ls.id,
    sex = sex.id
  )
  ## remove entries with useless content
  res <- res[!(is.na(res[, 1]) & is.na(res[, 2]) & is.na(res[, 3]) & is.na(res[, 4]) & is.na(res[, 5])), ]
  ## construct output stems for deliverables, if sufficient information is present to fit the accepted format
  output.stems <- construct.output.stems(res)
  res$output <- output.stems
  write.table(res, output.fn, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

if (exists("snakemake")) {
  run.construct.linker(
    snakemake@input[["logbook"]],
    snakemake@output[["linker"]]
  )
}
