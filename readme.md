---
title: "Phylogenetic Mycorrhizal Nutrient Economic Trait Analysis"
output: html_document
---

## Code to replicate analyses and figures from Averill et al. 2019. (Journal Name)
This repository contains all code to replicate analyses and figures from *Averill et al. 2019. Global imprint of mycorrhizal fungi on whole-plant nutrient economics. (Journal Name).* In order to run this code, you will need three data files: a phylogeny file, an intraspecific trait file and an interspecific trait file. Because of data sharing restrictions imposed by data providers, this data in only available upon request.

The file `paths.r` is sourced at the top of all analysis and figures scripts. This file tracks all file paths within the project. To run this code, the user just needs to change the first line of code in the `paths.r` file to match the location of the data directory on their computer.

Analysis files are numbered 1-4 and need to be run in order. This is because  some analyses depend on the output of other analyses.

This code assumes the working directory is set to the repository containing this R project.

These analyses depend on a number of custom functions, found within the `functions/` directory.

All analysis output is stored in the data directory specified in the `paths.r` script. All figures are stored within the project repository sub directory `figures/`.

This project was built under R version 3.4.4 and depends on the following packages:
`phytools`
`caper`
`picante`
`ggplot2`
`ggalt`
`data.table`
`MCMCglmm`
`lme4`
`MuMIn`