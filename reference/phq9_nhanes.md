# PHQ-9 Depression Screener (NHANES Random Sample)

A simple random sample of 1000 respondents to the Patient Health
Questionnaire 9-item (PHQ-9) depression screener in the U.S. National
Health and Nutrition Examination Survey (NHANES), September 2024 release
(file DPQ_L). Respondents were sampled from the 5455 who answered all
nine items with a valid response (0 to 3). There is no selection on
scores, so the data keep the many zero responses of a general population
sample: 21.6% have a sum score of zero (24.3% among all complete
responders). Intended for the hurdle partial credit model, see
[`hurdle_acat`](https://pgmj.github.io/easyRaschBayes/reference/hurdle_acat.md).
For a sample with higher symptom levels, see
[`phq9`](https://pgmj.github.io/easyRaschBayes/reference/phq9.md).

## Usage

``` r
phq9_nhanes
```

## Format

A data frame with 1000 rows and 10 variables:

- q1:

  Little interest or pleasure in doing things. Integer 0–3.

- q2:

  Feeling down, depressed, or hopeless. Integer 0–3.

- q3:

  Trouble falling/staying asleep, or sleeping too much. Integer 0–3.

- q4:

  Feeling tired or having little energy. Integer 0–3.

- q5:

  Poor appetite or overeating. Integer 0–3.

- q6:

  Feeling bad about yourself, or that you are a failure or have let
  yourself or your family down. Integer 0–3.

- q7:

  Trouble concentrating on things, such as reading the newspaper or
  watching television. Integer 0–3.

- q8:

  Moving or speaking so slowly that other people could have noticed, or
  the opposite, being so fidgety or restless that you have been moving
  around a lot more than usual. Integer 0–3.

- q9:

  Thoughts that you would be better off dead, or of hurting yourself in
  some way. Integer 0–3.

- seqn:

  NHANES respondent sequence number (SEQN), for linking to other NHANES
  files of the same cycle.

Each item uses a four-point ordinal response scale, scored 0
(`"Not at all"`), 1 (`"Several days"`), 2 (`"More than half the days"`)
and 3 (`"Nearly every day"`). The highest categories are sparse for some
items, for example 1% of respondents chose 3 on `q9`.

## Source

U.S. Centers for Disease Control and Prevention, National Center for
Health Statistics. *National Health and Nutrition Examination Survey*,
September 2024 release, file DPQ_L.
<https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component=Questionnaire&CycleBeginYear=2024>.
NHANES data are released to the public domain by the U.S. federal
government (<https://www.cdc.gov/nchs/policy/data-release-policy.html>).

## Details

Codes 7 (refused) and 9 (don't know) were treated as missing, and
respondents with any missing item were left out before sampling. The
NHANES survey weights and design are ignored, so the data do not give
population estimates. The script `data-raw/phq9_nhanes.R` in the package
source rebuilds the dataset from the NHANES file.

## References

Kroenke, K., Spitzer, R. L., & Williams, J. B. W. (2001). The PHQ-9:
Validity of a brief depression severity measure. *Journal of General
Internal Medicine*, 16(9), 606–613.
[doi:10.1046/j.1525-1497.2001.016009606.x](https://doi.org/10.1046/j.1525-1497.2001.016009606.x)

## Examples

``` r
data(phq9_nhanes)
table(rowSums(phq9_nhanes[, 1:9]) == 0)
#> 
#> FALSE  TRUE 
#>   784   216 
```
