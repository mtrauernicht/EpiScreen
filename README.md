# Epigenetic Screening on TRIP Clone
In this repository all code regarding the epigenetic drug screening on the K562 CRISPR-Cas9 DSB clone#5 with 19 integrations is uploaded. 

In short we tested 162 epigenetic drugs for their effect on Cas9 cutting efficiency and DNA repair pathway balance. 

![Figure 1. Overview of the reporter system and screen data.](Figure%201.png)

All relevant processed files can be found here: zenodo.org and in the files_scripts folder.

## Code is structured in the following script units:
**Script #1: Parsing QC**\
General quality checks of the sequencing data. Includes read counts and viability check.

**Scripts #2: Preprocessing drug screen and validations**\
Mutating sequencing data to large data frame. Adding all relevant annotations.

**Script #3: Synergy calculations**\
Scripts related to the calculation of the chromatin context dependencies.

**Scripts #4: Statistics & Plotting**\
All the plots are generated with these scripts, split up by figure.

[![DOI](https://zenodo.org/badge/298230427.svg)](https://zenodo.org/badge/latestdoi/298230427)

