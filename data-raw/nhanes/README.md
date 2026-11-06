# NHANES source files

- `DPQ_L_2024.rds`: NHANES Mental Health - Depression Screener (DPQ_L.xpt,
  September 2024 release), read into R and saved unchanged as .rds.
  Downloaded from
  https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Questionnaire&CycleBeginYear=2024 6337 rows: SEQN, DPQ010 to DPQ100.

`data-raw/phq9_nhanes.R` reads it and writes `data/phq9_nhanes.rda`. The
script header documents the selection (complete responders, seed, sample
size, no selection on scores, survey weights ignored).

NHANES public-use data are in the public domain
(https://www.cdc.gov/nchs/policy/data-release-policy.html).
