# Prepare the bundled phq9_nhanes dataset.
#
# Source -----------------------------------------------------------------
# NHANES Mental Health - Depression Screener (file DPQ_L), September 2024
# release:
#   https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Questionnaire&CycleBeginYear=2024
# The SAS transport file DPQ_L.xpt was downloaded, read into R and saved
# unchanged as data-raw/nhanes/DPQ_L_2024.rds (6337 rows: SEQN, DPQ010 to
# DPQ100). NHANES public-use data are in the public domain
# (https://www.cdc.gov/nchs/policy/data-release-policy.html).
#
# Selection ----------------------------------------------------------------
# 1. The nine PHQ-9 items DPQ010-DPQ090. DPQ100 (how difficult the
#    problems have made things) is not part of the PHQ-9 sum score and is
#    left out.
# 2. Codes 7 (refused) and 9 (don't know) are set to missing.
# 3. Respondents with responses 0-3 on all nine items are kept (5455 of
#    6337: 818 have no PHQ-9 data, 64 have at least one missing, refused
#    or don't know item).
# 4. A simple random sample of N_SAMPLE of them, with a fixed seed. There is
#    no selection on scores, so the share of zero sum scores is that of the
#    complete responders (about 24%).
#
# NHANES survey weights and design are ignored, so the data do not give
# population estimates.
#
# To regenerate data/phq9_nhanes.rda, run from the package root:
#   Rscript data-raw/phq9_nhanes.R

N_SAMPLE <- 1000L
SEED     <- 20241106L

raw <- as.data.frame(readRDS("data-raw/nhanes/DPQ_L_2024.rds"))
items_raw <- sprintf("DPQ0%d0", 1:9)
stopifnot(nrow(raw) == 6337L, all(c("SEQN", items_raw) %in% names(raw)))

resp <- as.matrix(raw[, items_raw])
resp[resp %in% c(7, 9)] <- NA
complete <- rowSums(is.na(resp)) == 0L
stopifnot(sum(complete) == 5455L, all(resp[complete, ] %in% 0:3))

pool <- data.frame(resp[complete, ], seqn = as.integer(raw$SEQN[complete]))
names(pool) <- c(paste0("q", 1:9), "seqn")

set.seed(SEED)
phq9_nhanes <- pool[sort(sample(nrow(pool), N_SAMPLE)), ]
for (i in 1:9) {
  phq9_nhanes[[paste0("q", i)]] <- as.integer(phq9_nhanes[[paste0("q", i)]])
}
rownames(phq9_nhanes) <- NULL

# Report the share of zero sum scores in the pool and in the sample
s_pool <- rowSums(pool[, 1:9])
s_samp <- rowSums(phq9_nhanes[, 1:9])
message(sprintf(
  "Zero sum scores: %.1f%% of %d complete responders, %.1f%% of the sample.",
  100 * mean(s_pool == 0), nrow(pool), 100 * mean(s_samp == 0)))

save(phq9_nhanes, file = "data/phq9_nhanes.rda", compress = "bzip2",
     version = 3)

loaded <- new.env()
load("data/phq9_nhanes.rda", envir = loaded)
stopifnot(identical(loaded$phq9_nhanes, phq9_nhanes))
message("data/phq9_nhanes.rda written (", nrow(phq9_nhanes), " rows, ",
        ncol(phq9_nhanes), " columns).")
