# Reference NCA with PKNCA for the nca.py comparison (see README.md in this folder).
# Usage: Rscript run_pknca.R <folder>   (reads single_dose.csv, steady_state.csv; writes pknca_*.csv)
suppressMessages(library(PKNCA))
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "."

# ICH M13A 2.2.2.2: BLQ = 0 in PK calculations. PKNCA's default drops BLQ values
# between quantifiable ones (conc.blq middle = "drop"); keep them as zero instead.
PKNCA.options(conc.blq = list(first = "keep", middle = "keep", last = "keep"))

run <- function(file, intervals, method) {
  d <- read.csv(file.path(dir, file), stringsAsFactors = FALSE)
  conc <- PKNCAconc(d, conc ~ time | id)
  dose <- PKNCAdose(data.frame(id = unique(d$id), time = 0, dose = 1), dose ~ time | id)
  res <- pk.nca(PKNCAdata(conc, dose, intervals = intervals, options = list(auc.method = method)))
  out <- as.data.frame(res$result)
  out$method <- method
  out[, c("id", "start", "end", "PPTESTCD", "PPORRES", "method")]
}

single <- data.frame(start = 0, end = Inf, cmax = TRUE, tmax = TRUE, tlast = TRUE, clast.obs = TRUE,
                     auclast = TRUE, lambda.z = TRUE, half.life = TRUE, lambda.z.n.points = TRUE,
                     adj.r.squared = TRUE, aucinf.obs = TRUE, aucinf.pred = TRUE)
partial <- data.frame(start = c(0, 0), end = c(2, 24), aucint.all = TRUE, aucint.last = TRUE)
ss <- data.frame(start = 0, end = 12, cmax = TRUE, tmax = TRUE, cmin = TRUE, ctrough = TRUE, aucint.last = TRUE,
                 cav.int.last = TRUE, swing = TRUE)

for (m in c("linear", "lin up/log down")) {
  tag <- if (m == "linear") "linear" else "linlog"
  write.csv(rbind(run("single_dose.csv", single, m), run("single_dose.csv", partial, m)),
            file.path(dir, paste0("pknca_single_", tag, ".csv")), row.names = FALSE)
  write.csv(run("steady_state.csv", ss, m), file.path(dir, paste0("pknca_ss_", tag, ".csv")), row.names = FALSE)
}
cat("PKNCA", as.character(packageVersion("PKNCA")), "R", R.version$major, R.version$minor, "\n")

# PKNCA's default BLQ handling (middle BLQ dropped), for comparison with the ICH M13A zero rule above.
PKNCA.options(default = TRUE)
write.csv(run("single_dose.csv", single, "linear"), file.path(dir, "pknca_single_linear_default_blq.csv"), row.names = FALSE)
