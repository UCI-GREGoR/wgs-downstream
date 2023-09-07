library(openxlsx)

run.create.manifest <- function(data.model.tsv, linker.fn, excluded.samples.fn, project.ids, sample.ids, out.fn) {
  data.model <- read.table(data.model.tsv,
    header = TRUE, stringsAsFactors = FALSE, sep = "\t",
    comment.char = "", quote = ""
  )
  linker <- read.table(linker.fn,
    header = TRUE, stringsAsFactors = FALSE, sep = "\t",
    comment.char = "", quote = ""
  )
  excluded.samples <- c()
  if (file.info(excluded.samples.fn)$size > 0) {
    excluded.samples <- read.table(excluded.samples.fn,
      header = FALSE,
      stringsAsFactors = FALSE, sep = "\t",
      comment.char = "", quote = ""
    )[, 1]
  }
  stopifnot(length(project.ids) == length(sample.ids))
  samplename <- c()
  samplestatus <- c()
  samplejson <- c()
  for (i in seq_len(length(project.ids))) {
    sample.index <- which(((linker$project == project.ids[i]) | is.na(linker$project)) & linker$index == sample.ids[i])
    ## some input bams are low-depth NA24385
    if (length(sample.index) != 1) {
      next
    }
    sample.id <- linker$subject[sample.index]
    project.id <- project.ids[i]
    if (sample.id %in% data.model[, "participant_id"] && !(sample.id %in% excluded.samples)) {
      affected.status <- data.model[data.model[, "participant_id"] == sample.id, "affected_status"]
      samplename <- c(samplename, sample.id)
      samplestatus <- c(
        samplestatus,
        ifelse(grepl("unaffected", affected.status, ignore.case = TRUE),
          "control", "case"
        )
      )
      samplejson <- c(samplejson, paste(getwd(), "/results/expansionhunter_denovo/profiles/", project.id,
        "/", sample.id, ".str_profile.json",
        sep = ""
      ))
    }
  }
  out.df <- data.frame(
    sample = samplename,
    status = samplestatus,
    json = samplejson
  )
  write.table(out.df, out.fn, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}

if (exists("snakemake")) {
  run.create.manifest(
    snakemake@input[["affected_status"]],
    snakemake@input[["linker"]],
    snakemake@input[["exclusions"]],
    snakemake@params[["projectids"]],
    snakemake@params[["sampleids"]],
    snakemake@output[["tsv"]]
  )
}
