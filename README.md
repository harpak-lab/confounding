# Archived analysis code for: **Benning, John W., Jedidiah Carlson, Olivia S. Smith, Ruth G. Shaw, and Arbel Harpak. "Confounding Fuels Misinterpretation in Human Genetics." bioRxiv (2024): 2023-11.**

This repository preserves the R workflows that accompany our article.  The documents `Clark2023.Rmd` and `SZ2024.R` will reproduce the statistical analyses and figures reported in the paper.

---

> **Note**   No raw data are committed to GitHub because they are covered by the original authors’ license.  Follow the steps in the next section to obtain them.

---

## Data sources

### Clark 2023 analyses

`clark_reportedCors.csv` contains the reported correlations in Clark (2023) Table 2 and is found in this repository.

All other data come from the supplementary datasets accompanying:

> **Clark, Gregory. "The inheritance of social status: England, 1600 to 2022." Proceedings of the National Academy of Sciences 120, no. 27 (2023): e2300926120.**

Supplementary data can be found at https://www.pnas.org/doi/10.1073/pnas.2300926120#supplementary-materials.

The supplementary data contain four Excel workbooks (Datasets 1–4) from which we extracted the following `.csv` files:

| `csv` file                | Source in Clark 2023 Supplementary data   | Contents                                   |
| ------------------------- | ----------------------------- | ------------------------------------------ |
| `clark_wealth.csv`        | Dataset 4, *Figure 3 - ded occ lwealth* tab       | Maternal / paternal wealth and status                   |
| `clark_lit.csv`           | Dataset 4, *Figure 3 - Literacy* tab   | Literacy data      |
| `clark_fatherSon.csv`     | Dataset 4, *Figure 4* tab   | Father-son data, with son age at father death                 |
| `clark_inds.csv`          | Dataset 1  | Individual-level data             |
| `clark_rels.csv`          | Dataset 2, *Table 2 1910‑97*  | Relative pair data for modern status measures                  |
| `clark_rels_occ.csv`      | Dataset 2, *Table 2 Occ Stat 1780-1919* | Occupational status data         |
| `clark_rels_ded.csv`      | Dataset 3, *Table 2 Ded 1780-1919*      | Higher education data             |
| `clark_rels_lit.csv`      | Dataset 3, *Table 2 Literacy 1754-1889* | Literacy data       |

### Song and Zhang 2024 analyses

All data come from the archived datasets found here: https://datadryad.org/dataset/doi:10.5061/dryad.4b8gthtk9

---



## Citing this code

If you build upon these scripts, please cite both our paper **and** the appropriate paper for the underlying data:

> Benning, John W., Jedidiah Carlson, Olivia S. Smith, Ruth G. Shaw, and Arbel Harpak. "Confounding Fuels Misinterpretation in Human Genetics." bioRxiv (2024): 2023-11.

> Clark, Gregory. "The inheritance of social status: England, 1600 to 2022." Proceedings of the National Academy of Sciences 120, no. 27 (2023): e2300926120.

> Song, S., and J. Zhang. 2024. Genetic variants underlying human bisexual behavior are reproductively advantageous. Sci. Adv. 10:eadj6958.


---

## Contact

Please email **[jbenning@cornell.edu](mailto:jbenning@cornell.edu)**.





