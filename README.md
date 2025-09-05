# Archived analysis code for: **Benning, John W., Jedidiah Carlson, Olivia S. Smith, Ruth G. Shaw, and Arbel Harpak. "Confounding Fuels Misinterpretation in Human Genetics." Proceedings of the Royal Society B (2025).**

**Abstract**: *The scientific literature has seen a resurgence of interest in genetic influences on human behavior and socioeconomic outcomes. Such studies face the central difficulty of distinguishing possible causal influences, in particular genetic and non-genetic ones. When confounding between possible influences is not rigorously addressed, it invites over- and misinterpretation of data. We illustrate the breadth of this problem through a discussion of the literature and a reanalysis of two examples. Clark (2023) suggested that patterns of similarity in social status between relatives indicate that social status is largely determined by one’s DNA. We show that the paper’s conclusions are based on the conflation of genetic and non-genetic transmission (for example, of wealth) within families. Song & Zhang (2024) posited that genetic variants underlying bisexual behavior are maintained in the population because they also affect risk-taking behavior and thereby confer an evolutionary fitness advantage through increased sexual promiscuity. In this case, too, we show that possible explanations cannot be distinguished, but only one is chosen and presented as a conclusion. We discuss how issues of confounding apply more broadly to studies that claim to establish genetic underpinnings to human behavior and societal outcomes.*

This repository preserves the R workflows that accompany our article. The documents `Clark2023.Rmd` and `SZ2024.R` will reproduce the statistical analyses and figures reported in the paper. These repository files are licensed under CC BY.

------------------------------------------------------------------------

> **Note**   No raw data are committed to GitHub because they are covered by the original authors’ CC BY license. Follow the steps in the next section to obtain them.

------------------------------------------------------------------------

## Workflow

Data can be downloaded from the archived data sets accompanying Clark 2023 and Song and Zhang 2024 (see below). Analyses based on data from Clark 2023 data are done in R using `Clark2023.rmd` and analyses based on data from Song and Zhang 2024 data are done in R using `SZ2024.R`.

## Data sources

### Clark 2023 analyses

`clark_reportedCors.csv` contains the reported correlations in Clark (2023) Table 2 and is found in this repository.

All other data come from the supplementary datasets accompanying:

> Clark, Gregory. "The inheritance of social status: England, 1600 to 2022." Proceedings of the National Academy of Sciences 120, no. 27 (2023): e2300926120.

Supplementary data can be found at <https://www.pnas.org/doi/10.1073/pnas.2300926120#supplementary-materials>. These datasets are licensed as CC-BY.

The supplementary data contain four Excel workbooks (Datasets 1–4) from which we extracted the following `.csv` files:

| `csv` file | Source in Clark 2023 Supplementary data | Contents |
|----|----|----|
| `clark_wealth.csv` | Dataset 4, *Figure 3 - ded occ lwealth* tab | Maternal / paternal wealth and status |
| `clark_lit.csv` | Dataset 4, *Figure 3 - Literacy* tab | Literacy data |
| `clark_fatherSon.csv` | Dataset 4, *Figure 4* tab | Father-son data, with son age at father death |
| `clark_inds.csv` | Dataset 1 | Individual-level data |
| `clark_rels.csv` | Dataset 2, *Table 2 1910‑97* | Relative pair data for modern status measures |
| `clark_rels_occ.csv` | Dataset 2, *Table 2 Occ Stat 1780-1919* | Occupational status data |
| `clark_rels_ded.csv` | Dataset 3, *Table 2 Ded 1780-1919* | Higher education data |
| `clark_rels_lit.csv` | Dataset 3, *Table 2 Literacy 1754-1889* | Literacy data |

#### Variable dictionary

This table lists the variables referenced in `Clark2023.Rmd` and used from the Clark (2023) supplementary datasets (Datasets 1–4), along with descriptions and units/scales as they are used in the analysis. See Clark (2023) for more details on construction of variables.

| Variable | Description | Units / Scale | Notes |
|----|----|----|----|
| `pid` | Person identifier | ID (integer) |  |
| `pidf`, `pidm`, `pidc` | Person IDs for father (`f`), mother (`m`), child (`c`) | ID |  |
| `pid_grandchild` | Person ID for grandchild | ID |  |
| `pidpgf`, `pidmgf` | Person IDs for paternal / maternal grandfather | ID |  |
| `relationship`, `relationship_rename` | Kinship relation between two records | Categorical | Harmonized across datasets, where relationship names were not standardized |
| `byr0`, `byr1` | Birth year of pair member 0 / 1 | Year (AD) |  |
| `nid` | Numeric lineage identifier | ID (integer) |  |
| `lwealth` | Log wealth | Unitless | Mean-centered log of estimated wealth |
| `lwealthgc`, `lwealthpgf`, `lwealthmgf` | Log wealth of grandchild / paternal GF / maternal GF | Unitless |  |
| `occ` | Occupational status score | Unitless index | Status scale as defined by Clark (2023) |
| `occ0`, `occ1` | Occupational status of pair member 0 / 1 | Unitless index |  |
| `occgc`, `occpgf`, `occmgf` | Occupational status of grandchild / paternal GF / maternal GF | Unitless index |  |
| `ded` | Higher-education attainment indicator | Binary (0/1) | 1 = degree/higher education |
| `ded0`, `ded1` | Higher-education indicator for pair member 0 / 1 | Binary (0/1) |  |
| `dedgc`, `dedpgf`, `dedmgf` | Higher-education indicator (grandchild / paternal GF / maternal GF) | Binary (0/1) |  |
| `lit` | Literacy indicator | Binary (0/1) | 1 = literate |
| `lit0`, `lit1` | Literacy indicator for pair member 0 / 1 | Binary (0/1) |  |
| `litf`, `litm`, `litc` | Literacy of father / mother / child | Binary (0/1) |  |
| `agedeathf` | Individual's age at father's death | Years |  |
| `imd` | Index of multiple deprivation | Index | Weighted average of measures of social deprivation |
| `lhv` | log(house value) | Unitless | Natural log of house value |
| `codir` | Company director indicator | Binary (0/1) | Indicator of whether individual held a company directorship |
| `statmod` | Modern status | Index | PCA index of lhv, imd, and codir |

### Song and Zhang 2024 analyses

All data come from the archived datasets accompanying,

> Song, S., and J. Zhang. 2024. Genetic variants underlying human bisexual behavior are reproductively advantageous. Sci. Adv. 10:eadj6958.

Data can be found here: <https://datadryad.org/dataset/doi:10.5061/dryad.4b8gthtk9>. These datasets are licensed as CC-BY-NC.

#### Variable dictionary

This table lists the variables referenced in `SZ2024.R`, along with descriptions and units/scales as they are used in the analysis. See Song and Zhang (2024) for more details on construction of variables.

| Variable | Description       | Units / Scale | Notes |
|----------|-------------------|---------------|-------|
| `pid`    | Person identifier | ID (integer)  |       |

## Citing this code

If you build upon these scripts, please cite both our paper **and** the appropriate paper for the underlying data:

> Benning, John W., Jedidiah Carlson, Olivia S. Smith, Ruth G. Shaw, and Arbel Harpak. "Confounding Fuels Misinterpretation in Human Genetics." Proceedings of the Royal Society B (2025).

> Clark, Gregory. "The inheritance of social status: England, 1600 to 2022." Proceedings of the National Academy of Sciences 120, no. 27 (2023): e2300926120.

> Song, S., and J. Zhang. 2024. Genetic variants underlying human bisexual behavior are reproductively advantageous. Sci. Adv. 10:eadj6958.

------------------------------------------------------------------------

## Session and package info

```         
 setting  value
 version  R version 4.5.1 (2025-06-13)
 os       macOS Sequoia 15.5
 system   aarch64, darwin20
 ui       RStudio
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/New_York
 date     2025-09-04
 rstudio  2025.05.1+513 Mariposa Orchid (desktop)
 pandoc   3.4 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
 quarto   1.6.42 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto


 package      * version date (UTC) lib source
 backports      1.5.0   2024-05-23 [1] CRAN (R 4.5.0)
 base64enc      0.1-3   2015-07-28 [1] CRAN (R 4.5.0)
 boot           1.3-31  2024-08-28 [2] CRAN (R 4.5.1)
 broom        * 1.0.8   2025-03-28 [1] CRAN (R 4.5.0)
 checkmate      2.3.2   2024-07-29 [1] CRAN (R 4.5.0)
 cli            3.6.5   2025-04-23 [1] CRAN (R 4.5.0)
 cluster        2.1.8.1 2025-03-12 [2] CRAN (R 4.5.1)
 colorspace     2.1-1   2024-07-26 [1] CRAN (R 4.5.0)
 cowplot      * 1.2.0   2025-07-07 [1] CRAN (R 4.5.0)
 data.table     1.17.6  2025-06-17 [1] CRAN (R 4.5.0)
 digest         0.6.37  2024-08-19 [1] CRAN (R 4.5.0)
 dplyr        * 1.1.4   2023-11-17 [1] CRAN (R 4.5.0)
 evaluate       1.0.4   2025-06-18 [1] CRAN (R 4.5.0)
 farver         2.1.2   2024-05-13 [1] CRAN (R 4.5.0)
 fastmap        1.2.0   2024-05-15 [1] CRAN (R 4.5.0)
 forcats      * 1.0.0   2023-01-29 [1] CRAN (R 4.5.0)
 foreign        0.8-90  2025-03-31 [2] CRAN (R 4.5.1)
 Formula        1.2-5   2023-02-24 [1] CRAN (R 4.5.0)
 generics       0.1.4   2025-05-09 [1] CRAN (R 4.5.0)
 ggbrace      * 0.1.2   2025-07-09 [1] CRAN (R 4.5.0)
 ggExtra      * 0.11.0  2025-09-01 [1] CRAN (R 4.5.0)
 ggplot2      * 3.5.2   2025-04-09 [1] CRAN (R 4.5.0)
 ggpmisc      * 0.6.2   2025-07-08 [1] CRAN (R 4.5.0)
 ggpp         * 0.5.9   2025-06-28 [1] CRAN (R 4.5.0)
 ggrepel      * 0.9.6   2024-09-07 [1] CRAN (R 4.5.0)
 glue           1.8.0   2024-09-30 [1] CRAN (R 4.5.0)
 gridExtra      2.3     2017-09-09 [1] CRAN (R 4.5.0)
 gtable         0.3.6   2024-10-25 [1] CRAN (R 4.5.0)
 here         * 1.0.1   2020-12-13 [1] CRAN (R 4.5.0)
 Hmisc        * 5.2-3   2025-03-16 [1] CRAN (R 4.5.0)
 hms            1.1.3   2023-03-21 [1] CRAN (R 4.5.0)
 htmlTable      2.4.3   2024-07-21 [1] CRAN (R 4.5.0)
 htmltools      0.5.8.1 2024-04-04 [1] CRAN (R 4.5.0)
 htmlwidgets    1.6.4   2023-12-06 [1] CRAN (R 4.5.0)
 httpuv         1.6.16  2025-04-16 [1] CRAN (R 4.5.0)
 insight        1.3.1   2025-06-30 [1] CRAN (R 4.5.0)
 knitr          1.50    2025-03-16 [1] CRAN (R 4.5.0)
 labeling       0.4.3   2023-08-29 [1] CRAN (R 4.5.0)
 later          1.4.2   2025-04-08 [1] CRAN (R 4.5.0)
 lattice        0.22-7  2025-04-02 [2] CRAN (R 4.5.1)
 lifecycle      1.0.4   2023-11-07 [1] CRAN (R 4.5.0)
 lme4         * 1.1-37  2025-03-26 [1] CRAN (R 4.5.0)
 lmtest       * 0.9-40  2022-03-21 [1] CRAN (R 4.5.0)
 lubridate    * 1.9.4   2024-12-08 [1] CRAN (R 4.5.0)
 magrittr       2.0.3   2022-03-30 [1] CRAN (R 4.5.0)
 MASS           7.3-65  2025-02-28 [2] CRAN (R 4.5.1)
 Matrix       * 1.7-3   2025-03-11 [2] CRAN (R 4.5.1)
 MatrixModels   0.5-4   2025-03-26 [1] CRAN (R 4.5.0)
 mgcv         * 1.9-3   2025-04-04 [2] CRAN (R 4.5.1)
 mime           0.13    2025-03-17 [1] CRAN (R 4.5.0)
 miniUI         0.1.2   2025-04-17 [1] CRAN (R 4.5.0)
 minqa          1.2.8   2024-08-17 [1] CRAN (R 4.5.0)
 moments      * 0.14.1  2022-05-02 [1] CRAN (R 4.5.0)
 nlme         * 3.1-168 2025-03-31 [2] CRAN (R 4.5.1)
 nloptr         2.2.1   2025-03-17 [1] CRAN (R 4.5.0)
 nnet           7.3-20  2025-01-01 [2] CRAN (R 4.5.1)
 performance  * 0.15.0  2025-07-10 [1] CRAN (R 4.5.0)
 pillar         1.11.0  2025-07-04 [1] CRAN (R 4.5.0)
 pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.5.0)
 polynom        1.4-1   2022-04-11 [1] CRAN (R 4.5.0)
 promises       1.3.3   2025-05-29 [1] CRAN (R 4.5.0)
 purrr        * 1.0.4   2025-02-05 [1] CRAN (R 4.5.0)
 quantreg       6.1     2025-03-10 [1] CRAN (R 4.5.0)
 R6             2.6.1   2025-02-15 [1] CRAN (R 4.5.0)
 rbibutils      2.3     2024-10-04 [1] CRAN (R 4.5.0)
 RColorBrewer * 1.1-3   2022-04-03 [1] CRAN (R 4.5.0)
 Rcpp           1.1.0   2025-07-02 [1] CRAN (R 4.5.0)
 Rdpack         2.6.4   2025-04-09 [1] CRAN (R 4.5.0)
 readr        * 2.1.5   2024-01-10 [1] CRAN (R 4.5.0)
 reformulas     0.4.1   2025-04-30 [1] CRAN (R 4.5.0)
 renv         * 1.1.5   2025-07-24 [1] CRAN (R 4.5.0)
 rlang          1.1.6   2025-04-11 [1] CRAN (R 4.5.0)
 rmarkdown      2.29    2024-11-04 [1] CRAN (R 4.5.0)
 rpart          4.1.24  2025-01-07 [2] CRAN (R 4.5.1)
 rprojroot      2.0.4   2023-11-05 [1] CRAN (R 4.5.0)
 rstudioapi     0.17.1  2024-10-22 [1] CRAN (R 4.5.0)
 sandwich     * 3.1-1   2024-09-15 [1] CRAN (R 4.5.0)
 scales       * 1.4.0   2025-04-24 [1] CRAN (R 4.5.0)
 sessioninfo  * 1.2.3   2025-02-05 [1] CRAN (R 4.5.0)
 shiny          1.11.1  2025-07-03 [1] CRAN (R 4.5.0)
 SparseM        1.84-2  2024-07-17 [1] CRAN (R 4.5.0)
 stringi        1.8.7   2025-03-27 [1] CRAN (R 4.5.0)
 stringr      * 1.5.1   2023-11-14 [1] CRAN (R 4.5.0)
 survival       3.8-3   2024-12-17 [2] CRAN (R 4.5.1)
 tibble       * 3.3.0   2025-06-08 [1] CRAN (R 4.5.0)
 tidyr        * 1.3.1   2024-01-24 [1] CRAN (R 4.5.0)
 tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.5.0)
 tidyverse    * 2.0.0   2023-02-22 [1] CRAN (R 4.5.0)
 timechange     0.3.0   2024-01-18 [1] CRAN (R 4.5.0)
 tzdb           0.5.0   2025-03-15 [1] CRAN (R 4.5.0)
 vctrs          0.6.5   2023-12-01 [1] CRAN (R 4.5.0)
 withr          3.0.2   2024-10-28 [1] CRAN (R 4.5.0)
 xfun           0.52    2025-04-02 [1] CRAN (R 4.5.0)
 xtable         1.8-4   2019-04-21 [1] CRAN (R 4.5.0)
 yaml           2.3.10  2024-07-26 [1] CRAN (R 4.5.0)
 zoo          * 1.8-14  2025-04-10 [1] CRAN (R 4.5.0)
```

------------------------------------------------------------------------

## Contact

Please email [**jbenning\@cornell.edu**](mailto:jbenning@cornell.edu).
