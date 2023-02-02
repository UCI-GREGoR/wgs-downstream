library(openxlsx)

run.create.manifest <- function(data.model.xlsx, linker.fn, project.ids, sample.ids, out.fn) {
  data.model <- openxlsx::read.xlsx(data.model.xlsx, sheet = "participant", check.names = FALSE)
  linker <- read.table(linker.fn,
    header = TRUE, stringsAsFactors = FALSE, sep = "\t",
    comment.char = "", quote = ""
  )
  stopifnot(length(project.ids) == length(sample.ids))
  samplename <- c()
  samplestatus <- c()
  samplejson <- c()
  for (i in seq_len(length(project.ids))) {
    pmgrc.index <- which(linker$ru == project.ids[i] & linker$sq == sample.ids[i])
    ## some input bams are low-depth NA24385
    if (length(pmgrc.index) != 1) {
      next
    }
    pmgrc.id <- linker$pmgrc[pmgrc.index]
    ru.id <- project.ids[i]
    if (pmgrc.id %in% data.model[, "participant_id"]) {
      affected.status <- data.model[data.model[, "participant_id"] == pmgrc.id, "affected_status"]
      samplename <- c(samplename, pmgrc.id)
      samplestatus <- c(
        samplestatus,
        ifelse(grepl("unaffected", affected.status, ignore.case = TRUE),
          "control", "case"
        )
      )
      samplejson <- c(samplejson, paste(getwd(), "/results/profiles/", ru.id,
        "/", pmgrc.id, ".str_profile.json",
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
    snakemake@input[["data_model"]],
    snakemake@input[["linker"]],
    snakemake@params[["projectids"]],
    snakemake@params[["sampleids"]],
    snakemake@output[["tsv"]]
  )
}
