run.create.manifest <- function(data.model.tsv, sex.manifest.fn, excluded.samples.fn, sample.ids, out.fn) {
  data.model <- read.table(data.model.tsv,
    header = TRUE, stringsAsFactors = FALSE, sep = "\t",
    comment.char = "", quote = ""
  )
  sex.manifest <- read.table(sex.manifest.fn,
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
  samplename <- c()
  samplestatus <- c()
  samplejson <- c()
  for (sample.id in sample.ids) {
    if (sample.id %in% data.model[, "participant_id"] && !(sample.id %in% excluded.samples)) {
      affected.status <- data.model[data.model[, "participant_id"] == sample.id, "affected_status"]
      samplename <- c(samplename, sample.id)
      samplestatus <- c(
        samplestatus,
        ifelse(grepl("unaffected", affected.status, ignore.case = TRUE),
          "control", "case"
        )
      )
      samplejson <- c(samplejson, paste(getwd(), "/results/expansionhunter_denovo/profiles/",
        sample.id, ".str_profile.json",
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
    snakemake@input[["sex_manifest"]],
    snakemake@input[["exclusions"]],
    snakemake@params[["sampleids"]],
    snakemake@output[["tsv"]]
  )
}
