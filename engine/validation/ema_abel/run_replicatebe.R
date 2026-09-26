# Reference ABEL results with replicateBE (EMA Method A) for the bioequivalence.py comparison.
# Usage: Rscript run_replicatebe.R <folder>. Writes rdsNN.csv (data) and replicatebe_method_a.csv (results).
suppressMessages(library(replicateBE))
args <- commandArgs(trailingOnly = TRUE); dir <- if (length(args)) args[1] else "."
sets <- c("rds01", "rds02", "rds03", "rds04", "rds05", "rds06", "rds07", "rds08", "rds09", "rds10", "rds11", "rds12",
          "rds13", "rds14", "rds15", "rds16", "rds17", "rds18", "rds19", "rds20", "rds21", "rds22", "rds23", "rds24",
          "rds25", "rds26", "rds27", "rds28", "rds29", "rds30")
out <- list()
for (s in sets) {
  d <- get(s, envir = asNamespace("replicateBE"))
  d <- d[, c("subject", "period", "sequence", "treatment", "PK")]
  write.csv(data.frame(subject = d$subject, sequence = d$sequence, period = d$period, treatment = d$treatment,
                       value = d$PK), file.path(dir, paste0(s, ".csv")), row.names = FALSE, na = "")
  r <- tryCatch(method.A(data = get(s, envir = asNamespace("replicateBE")), print = FALSE, details = TRUE, verbose = FALSE),
                error = function(e) NULL)
  if (!is.null(r)) {
    g <- function(k) if (is.null(r[[k]])) NA else r[[k]]
    out[[s]] <- data.frame(set = s, design = g('Design'), n = g('n'), missing = sum(is.na(d$PK)),
                           CVwR = g('CVwR(%)'), EL.lo = if (is.na(g('L(%)'))) g('BE.lo(%)') else g('L(%)'), EL.hi = if (is.na(g('U(%)'))) g('BE.hi(%)') else g('U(%)'), PE = g('PE(%)'),
                           CL.lo = g('CL.lo(%)'), CL.hi = g('CL.hi(%)'), BE = g('BE'))
  }
}
# Incomplete datasets: the same method on complete subjects only (bioequivalence.py accepts complete designs only)
for (s in sets) {
  d <- get(s, envir = asNamespace("replicateBE"))
  d <- d[!is.na(d$PK), ]
  np <- nlevels(droplevels(factor(d$period)))
  keep <- names(which(table(droplevels(factor(d$subject))) == np))
  if (length(keep) == nlevels(droplevels(factor(d$subject)))) next
  dc <- droplevels(d[as.character(d$subject) %in% keep, ])
  nm <- paste0(s, "_complete")
  write.csv(data.frame(subject = dc$subject, sequence = dc$sequence, period = dc$period, treatment = dc$treatment,
                       value = dc$PK), file.path(dir, paste0(nm, ".csv")), row.names = FALSE, na = "")
  r <- tryCatch(method.A(data = dc, print = FALSE, details = TRUE, verbose = FALSE), error = function(e) NULL)
  if (!is.null(r)) {
    g <- function(k) if (is.null(r[[k]])) NA else r[[k]]
    out[[nm]] <- data.frame(set = nm, design = g('Design'), n = g('n'), missing = 0,
                            CVwR = g('CVwR(%)'), EL.lo = if (is.na(g('L(%)'))) g('BE.lo(%)') else g('L(%)'),
                            EL.hi = if (is.na(g('U(%)'))) g('BE.hi(%)') else g('U(%)'), PE = g('PE(%)'),
                            CL.lo = g('CL.lo(%)'), CL.hi = g('CL.hi(%)'), BE = g('BE'))
  }
}
write.csv(do.call(rbind, out), file.path(dir, "replicatebe_method_a.csv"), row.names = FALSE)
cat("replicateBE", as.character(packageVersion("replicateBE")), "\n")
